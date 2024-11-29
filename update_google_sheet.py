# update_google_sheet.py
import os
import ast
import json
from google.oauth2 import service_account
from googleapiclient.discovery import build

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
            func_info['docstring'] = ast.get_docstring(node)
            func_info['parameters'] = []
            for arg in node.args.args:
                arg_name = arg.arg
                arg_type = ast.get_source_segment(source, arg.annotation) if arg.annotation else ''
                func_info['parameters'].append(f"{arg_name}: {arg_type}")
            func_info['returns'] = ast.get_source_segment(source, node.returns) if node.returns else ''
            functions.append(func_info)
    return functions

def update_google_sheet(functions_data):
    service = build('sheets', 'v4', credentials=credentials)
    sheet = service.spreadsheets()

    # Prepare data for the sheet
    data = [['Module', 'Function Name', 'Description', 'Parameters', 'Returns']]
    for func in functions_data:
        row = [
            func['module'],
            func['name'],
            func['docstring'] or '',
            '\n'.join(func['parameters']),
            func['returns']
        ]
        data.append(row)

    # Write data to the Google Sheet
    body = {
        'values': data
    }
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
                functions = get_functions_from_file(path)
                functions_data.extend(functions)
    update_google_sheet(functions_data)

if __name__ == '__main__':
    main()
