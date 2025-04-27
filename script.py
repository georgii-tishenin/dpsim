import nbformat
import os

def convert_notebook_to_py_if_assert_found(notebook_path, test_folder):
    # Extract the notebook name without the extension
    notebook_name = os.path.splitext(os.path.basename(notebook_path))[0]

    # Load the notebook content
    with open(notebook_path, 'r') as nb_file:
        notebook_content = nb_file.read()
    notebook = nbformat.reads(notebook_content, as_version=4)

    # Check if any assert statement is found in code cells
    assert_found = False
    for cell in notebook.cells:
        if cell.cell_type == 'code':
            if 'assert' in cell.source:
                assert_found = True
                break

    # If an assert statement is found, convert notebook to .py file
    if assert_found:
        python_script = ''
        for cell in notebook.cells:
            if cell.cell_type == 'code':
                python_script += cell.source + '\n'
        
        # Ensure the 'test' folder exists in the root directory
        os.makedirs(test_folder, exist_ok=True)

        # Define the path for the new Python script in the 'test' folder
        py_file_path = os.path.join(test_folder, f"test_{notebook_name}.py")
        
        # Write the Python script to the file in the 'test' folder
        with open(py_file_path, 'w') as py_file:
            py_file.write(python_script)
        
        print(f"Conversion successful! Saved Python script as: {py_file_path}")
    else:
        print(f"No assert statement found in {notebook_path}. No conversion made.")

def convert_notebooks_in_folder(folder_path, test_folder):
    # Loop through all files in the folder
    for root, _, files in os.walk(folder_path):
        for file in files:
            if file.endswith('.ipynb'):  # Check for notebook files
                notebook_path = os.path.join(root, file)
                convert_notebook_to_py_if_assert_found(notebook_path, test_folder)

# Usage
folder_path = 'examples/Notebooks/Circuits/'  # Replace with the path to your folder
test_folder = 'test'  # Folder to save the converted Python scripts
convert_notebooks_in_folder(folder_path, test_folder)