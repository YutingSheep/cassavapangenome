import os
import argparse
from collections import defaultdict

# 读取文件并合并数据
def read_file(file_path):
    data = defaultdict(list)
    file_name = os.path.basename(file_path).split('.')[0]  # 提取文件名（去除后缀）作为物种名
    with open(file_path, 'r') as file:
        for line in file:
            parts = line.strip().split(',')
            key = parts[1]  # 第二列作为键
            data[key].append((file_name, parts))  # 存储文件名和行数据
    return data

# 拆分最后一列
def split_last_column(row):
    last_col = row[-1]
    if last_col.startswith('1'):
        return last_col.split(' ')[-1], 'None'
    elif last_col.startswith('B'):
        genes = last_col.split(' ')[-2:]
        return genes[0], genes[1]
    elif last_col.startswith('2'):
        return 'None', last_col.split(' ')[-1]
    else:
        return 'None', 'None'

# 合并多个文件
def merge_files(file_paths):
    merged_data = defaultdict(dict)
    for file_path in file_paths:
        data = read_file(file_path)
        for key, rows in data.items():
            for file_name, row in rows:
                left_gene, right_gene = split_last_column(row)
                if key not in merged_data:
                    merged_data[key] = {}
                merged_data[key][f"{file_name}_Left"] = left_gene
                merged_data[key][f"{file_name}_Right"] = right_gene
    return merged_data

# 输出合并后的结果
def output_merged_data(merged_data, output_file):
    # 获取所有列标题
    columns = set()
    for key, values in merged_data.items():
        columns.update(values.keys())
    columns = sorted(columns)  # 按字母顺序排序

    # 写入标题行
    with open(output_file, 'w') as file:
        file.write("Key," + ",".join(columns) + "\n")
        for key, values in merged_data.items():
            row = [key]
            for col in columns:
                row.append(values.get(col, 'None'))
            file.write(",".join(row) + "\n")

# 获取文件夹中的所有文件路径
def get_files_in_folder(folder_path):
    file_paths = []
    for root, dirs, files in os.walk(folder_path):
        for file in files:
            file_paths.append(os.path.join(root, file))
    return file_paths

# 主函数
def main():
    # 设置命令行参数
    parser = argparse.ArgumentParser(description="合并指定文件夹下的所有文件，并根据规则拆分最后一列数据。")
    parser.add_argument("folder", help="包含输入文件的文件夹路径")
    parser.add_argument("output", help="输出文件的路径")
    args = parser.parse_args()

    # 获取文件夹中的所有文件路径
    file_paths = get_files_in_folder(args.folder)

    # 合并文件
    merged_data = merge_files(file_paths)

    # 输出合并后的结果
    output_merged_data(merged_data, args.output)
    print(f"合并完成，结果已保存到 {args.output}")

if __name__ == "__main__":
    main()
