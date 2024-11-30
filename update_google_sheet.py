# update_google_sheet.py
import os
import ast
import json
from google.oauth2 import service_account
from googleapiclient.discovery import build
import docstring_parser  # Import the docstring_parser library

# Load service account credentials from environment variable
SERVICE_ACCOUNT_INFO = os.environ.get('GOOGLE_SERVICE_ACCOUNT_JSON')
if not SERVICE_ACCOUNT_INFO:
    raise ValueError("Service account credentials not found.")

service_account_info = json.loads(SERVICE_ACCOUNT_INFO)
SCOPES = ['https://www.googleapis.com/auth/spreadsheets']
credentials = service_account.Credentials.from_service_account_info(
    service_account_info, scopes=SCOPES)

SPREADSHEET_ID = os.environ.get('SPREADSHEET_ID')  # Set this as a GitHub Secret

def get_functions_from_file(file_path):
    with open(file_path, 'r', encoding='utf-8') as f:
        source = f.read()
    tree = ast.parse(source)
    functions = []

    for node in ast.walk(tree):
        if isinstance(node, ast.FunctionDef):
            func_info = {}
            func_info['module'] = os.path.basename(file_path)
            func_info['name'] = node.name
            func_info['parameters'] = []  # List of dicts with 'name', 'type', and 'description'
            func_info['return_type'] = ''
            func_info['return_description'] = ''
            func_info['description'] = ''
            # Get the function's docstring
            docstring = ast.get_docstring(node)
            if docstring:
                # Parse the docstring
                parsed_docstring = docstring_parser.parse(docstring)
                func_info['description'] = parsed_docstring.short_description or ''
                # Extract parameters
                for param in parsed_docstring.params:
                    param_name = param.arg_name or ''
                    param_type = param.type_name or ''
                    param_description = param.description or ''
                    func_info['parameters'].append({
                        'name': param_name,
                        'type': param_type,
                        'description': param_description
                    })
                # Extract return type and description
                if parsed_docstring.returns:
                    func_info['return_type'] = parsed_docstring.returns.type_name or ''
                    func_info['return_description'] = parsed_docstring.returns.description or ''
            else:
                func_info['description'] = ''
            functions.append(func_info)
    return functions

def update_google_sheet(functions_data):
    service = build('sheets', 'v4', credentials=credentials)
    sheet = service.spreadsheets()

    # Prepare data for the sheet
    data = [['Module', 'Function Name', 'Description',
             'Return Type', 'Return Description',
             'Parameters', 'Parameter Types', 'Parameter Descriptions']]
    for func in functions_data:
        # Separate parameters into names, types, and descriptions
        param_names = '\n'.join([p['name'] for p in func['parameters']])
        param_types = '\n'.join([p['type'] for p in func['parameters']])
        param_descriptions = '\n'.join([p['description'] for p in func['parameters']])

        row = [
            func['module'],
            func['name'],
            func['description'],
            func['return_type'],
            func['return_description'],
            param_names,
            param_types,
            param_descriptions
        ]
        data.append(row)

    # Write data to the Google Sheet
    body = {
        'values': data
    }
    # Clear existing data before updating
    sheet.values().clear(
        spreadsheetId=SPREADSHEET_ID,
        range='Sheet1!A1:Z1000'
    ).execute()

    sheet.values().update(
        spreadsheetId=SPREADSHEET_ID,
        range='Sheet1!A1',
        valueInputOption='RAW',
        body=body
    ).execute()

def main():
    functions_data = []
    for root, dirs, files in os.walk('.'):
        for file in files:
            if file.endswith('.py'):
                path = os.path.join(root, file)
                # Exclude certain directories and files if necessary
                if (
                    'venv' in path or
                    'tests' in path or
                    '.git' in path or
                    os.path.basename(path) == 'update_google_sheet.py'
                ):
                    continue
                functions = get_functions_from_file(path)
                functions_data.extend(functions)
    update_google_sheet(functions_data)

if __name__ == '__main__':
    main()
