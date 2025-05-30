#!/usr/bin/env python3

import json
import os

project_path = './'
output_file_path = os.path.normpath(os.path.abspath('output.json'))
script_file_path = os.path.normpath(os.path.abspath(__file__))

def load_config(config_path):
    with open(config_path, 'r') as config_file:
        return json.load(config_file)


config = {
    'ignore_extensions': [],
    'ignore_path': []
}

def is_binary_file(file_path):
    try:
        with open(file_path, 'r', encoding='utf-8') as file:
            if '\0' in file.read(1024):
                return True
    except Exception as e:
        print(f"Cannot read file {file_path} due to {e}, assuming binary.")
        return True
    return False

def should_ignore(file_path, config):
    normalized_file_path = os.path.normpath(os.path.abspath(file_path))
    if normalized_file_path == output_file_path or normalized_file_path == script_file_path:
        print('skip', normalized_file_path)
        return True
        
    for ignore_ext in (config['ignore_extensions']):
        if file_path.endswith(ignore_ext):
            # print('skip')
            return True
    return False

def format_file_content(file_path, is_binary):
    if is_binary:
        return file_path, "- BINARY FILE -"
    else:
        try:
            with open(file_path, 'r', encoding='utf-8', errors='ignore') as file:
                content = file.read()
            return file_path, content
        except Exception as e:
            return file_path, f"- ERROR READING FILE: {e} -"

def traverse_and_extract(project_path, config):
    file_contents = {}
    for root, dirs, files in os.walk(project_path):
        dirs[:] = [d for d in dirs if not should_ignore(os.path.join(root, d), config)]
        for filename in files:
            file_path = os.path.join(root, filename)
            if should_ignore(file_path, config):
                continue
            is_binary = is_binary_file(file_path)
            path, content = format_file_content(file_path, is_binary)
            file_contents[path] = content
    return file_contents

extracted_data = traverse_and_extract(project_path, config)


with open(output_file_path, 'w', encoding='utf-8') as output_file:
    json.dump(extracted_data, output_file, ensure_ascii=False, indent=4)
