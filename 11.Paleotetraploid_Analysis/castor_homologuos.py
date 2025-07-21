import pandas as pd
import sys

def merge_mdNC_genes(input_file, output_file):
    # 读取输入文件，尝试自动检测分隔符
    df = pd.read_csv(input_file, sep=None, engine='python')

    # 打印列名以确认
    print("Input file columns:", df.columns)

    # 检查是否有 'id' 列
    if 'id' not in df.columns:
        raise KeyError("Column 'id' not found in the input file. Please check the column names.")

    # 创建一个字典，用于存储合并后的结果
    merged_data = []

    # 按 id 分组
    for id_val, group in df.groupby('id'):
        # 获取 species=1 和 species=2 的数据
        species1 = group[group['species'] == 1]
        species2 = group[group['species'] == 2]

        # 遍历 species=1 的行
        for _, row1 in species1.iterrows():
            # 检查是否有 mdNC 基因
            if pd.notna(row1['castor_gene']) and row1['castor_gene'].startswith('mdNC'):
                mdNC = row1['castor_gene']
                # 查找 species=2 中是否有相同的 mdNC 基因
                row2 = species2[species2['castor_gene'] == mdNC]
                if not row2.empty:
                    # 如果找到，合并到一行
                    merged_data.append({**row1, **{f'{col}_2': row2.iloc[0][col] for col in row1.index}})
                else:
                    # 如果没找到，右边为空
                    merged_data.append({**row1, **{f'{col}_2': '' for col in row1.index}})
            else:
                # 如果没有 mdNC 基因，直接保留，右边为空
                merged_data.append({**row1, **{f'{col}_2': '' for col in row1.index}})

        # 遍历 species=2 的行，检查是否有未匹配的 mdNC 基因
        for _, row2 in species2.iterrows():
            if pd.notna(row2['castor_gene']) and row2['castor_gene'].startswith('mdNC'):
                mdNC = row2['castor_gene']
                # 检查是否已经处理过
                if not any(mdNC == row['castor_gene'] for row in merged_data if row['id'] == id_val):
                    # 如果没有处理过，将 species=2 的 mdNC 信息放到 species=1 的 castor_gene 列中
                    merged_data.append({
                        'id': id_val,
                        'species': 1,
                        'chr': '',
                        'start': '',
                        'end': '',
                        'gene': '',
                        'castor_gene': mdNC,
                        **{f'{col}_2': row2[col] for col in row2.index}
                    })

    # 将结果转换为 DataFrame
    merged_df = pd.DataFrame(merged_data)

    # 保存到输出文件
    merged_df.to_csv(output_file, sep=',', index=False)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python merge_mdNC_genes.py <input_file> <output_file>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]
    merge_mdNC_genes(input_file, output_file)
