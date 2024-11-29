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
            func_info['parameters'] = []
            func_info['returns'] = ''
            func_info['description'] = ''
            # Get the function's docstring
            docstring = ast.get_docstring(node)
            if docstring:
                # Parse the docstring
                parsed_docstring = docstring_parser.parse(docstring)
                func_info['description'] = parsed_docstring.short_description or ''
                # Extract parameters
                for param in parsed_docstring.params:
                    param_str = f"{param.arg_name}: {param.type_name or ''} - {param.description or ''}"
                    func_info['parameters'].append(param_str)
                # Extract return type and description
                if parsed_docstring.returns:
                    func_info['returns'] = f"{parsed_docstring.returns.type_name or ''} - {parsed_docstring.returns.description or ''}"
            else:
                func_info['description'] = ''
            functions.append(func_info)
    return functions

def update_google_sheet(functions_data):
    service = build('sheets', 'v4', credentials=credentials)
    sheet = service.spreadsheets()

    # Prepare data for the sheet
    # Modified the order of 'Parameters' and 'Returns'
    data = [['Module', 'Function Name', 'Description', 'Returns', 'Parameters']]
    for func in functions_data:
        row = [
            func['module'],
            func['name'],
            func['description'],
            func['returns'],
            '\n'.join(func['parameters'])
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
                # Exclude certain directories if necessary
                if 'venv' in path or 'tests' in path or '.git' in path:
                    continue
                functions = get_functions_from_file(path)
                functions_data.extend(functions)
    update_google_sheet(functions_data)

if __name__ == '__main__':
    main()
