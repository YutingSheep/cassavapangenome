import pandas as pd
import argparse

# 解析命令行参数
def parse_arguments():
    parser = argparse.ArgumentParser(description='Expand gene ranges into individual gene IDs.')
    parser.add_argument('input_file', type=str, help='Path to the input CSV file containing the original data.')
    parser.add_argument('output_file', type=str, help='Path to save the expanded gene data CSV file.')
    return parser.parse_args()

# 读取数据并展开基因编号
def expand_genes(input_file, output_file):
    # 读取原始数据
    df = pd.read_csv(input_file)
    
    # 用于保存转换后的新数据
    expanded_rows = []
    
    # 遍历每一行数据
    for index, row in df.iterrows():
        # 获取物种1的起始和终止基因编号
        start1, end1 = int(row['start1']), int(row['end1'])
        # 获取物种2的起始和终止基因编号
        start2, end2 = int(row['start2']), int(row['end2'])
        
        # 物种1的基因编号展开
        if start1 <= end1:
            gene_range1 = range(start1, end1 + 1)  # 正序
        else:
            gene_range1 = range(start1, end1 - 1, -1)  # 倒序
        
        # 物种2的基因编号展开
        if start2 <= end2:
            gene_range2 = range(start2, end2 + 1)  # 正序
        else:
            gene_range2 = range(start2, end2 - 1, -1)  # 倒序

        # 生成物种1的记录
        for gene in gene_range1:
            expanded_rows.append([row['id'], 1, row['chr1'], start1, end1, gene])

        # 生成物种2的记录
        for gene in gene_range2:
            expanded_rows.append([row['id'], 2, row['chr2'], start2, end2, gene])

    # 将展开的数据转换为DataFrame
    expanded_df = pd.DataFrame(expanded_rows, columns=['id', 'species', 'chr', 'start', 'end', 'gene'])

    # 将结果保存为CSV文件
    expanded_df.to_csv(output_file, index=False)
    print(f"Expanded data has been saved to {output_file}")

if __name__ == "__main__":
    args = parse_arguments()
    expand_genes(args.input_file, args.output_file)
