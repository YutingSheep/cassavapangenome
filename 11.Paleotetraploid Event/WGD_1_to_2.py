import pandas as pd
import argparse

# 创建解析器
parser = argparse.ArgumentParser(description="Filter rows where genes in column 3 have two unique rows in column 1.")
parser.add_argument("input_file", help="Path to the input file.")
parser.add_argument("output_file", help="Path to the output file.")

# 解析参数
args = parser.parse_args()
input_file = args.input_file
output_file = args.output_file

# 读取输入文件
df = pd.read_csv(input_file, sep=" ", header=None, names=["col1", "col2", "col3", "col4", "col5"])

# 统计第三列基因的对应关系
grouped = df.groupby("col3").filter(lambda x: len(x) == 2 and x["col1"].nunique() == 2)

# 保存结果到文件
grouped.to_csv(output_file, sep="\t", index=False, header=False)

print(f"Filtered data saved to {output_file}")
