import os
import json
import re
import nbformat

def extract_code_from_ipynb(file_path):
    """从.ipynb文件中提取代码"""
    with open(file_path, 'r', encoding='utf-8') as f:
        nb = nbformat.read(f, as_version=4)
    
    code_cells = []
    for cell in nb.cells:
        if cell.cell_type == 'code':
            code_cells.append(cell.source)
    return '\n'.join(code_cells)

def extract_code_from_rmd(file_path):
    """从.Rmd文件中提取代码"""
    with open(file_path, 'r', encoding='utf-8') as f:
        content = f.read()
    
    # 使用正则表达式提取代码块
    code_blocks = re.findall(r'```{r.*?}(.*?)```', content, re.DOTALL)
    return '\n'.join(code_blocks)

def extract_code(file_path):
    """根据文件扩展名调用相应的处理函数"""
    if file_path.endswith('.ipynb'):
        return extract_code_from_ipynb(file_path)
    elif file_path.endswith('.Rmd'):
        return extract_code_from_rmd(file_path)
    else:
        raise ValueError("Unsupported file type. Only .ipynb and .Rmd files are supported.")

def main():
    import argparse
    parser = argparse.ArgumentParser(description="Extract code from .ipynb or .Rmd files.")
    parser.add_argument("file_path", type=str, help="Path to the .ipynb or .Rmd file.")
    parser.add_argument("output_path", type=str, help="Path to save the extracted code.")
    
    args = parser.parse_args()
    
    file_path = args.file_path
    output_path = args.output_path
    
    if not os.path.exists(file_path):
        print(f"File not found: {file_path}")
        return
    
    try:
        code = extract_code(file_path)
        with open(output_path, 'w', encoding='utf-8') as f:
            f.write(code)
        print(f"Code extracted and saved to {output_path}")
    except Exception as e:
        print(f"Error processing file: {e}")

if __name__ == "__main__":
    main()