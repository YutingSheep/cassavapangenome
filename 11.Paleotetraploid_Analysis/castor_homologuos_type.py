import pandas as pd
import sys

def process_files(file_a_path, file_b_path, output_path):
    # 读取文件A
    file_a = pd.read_csv(file_a_path)

    # 读取文件B
    file_b = pd.read_csv(file_b_path, sep='\t', header=None)

    # 创建一个字典来存储文件A中的信息
    castor_gene_info = {}

    # 遍历文件A，填充字典
    for index, row in file_a.iterrows():
        castor_gene = row['castor_gene']
        gene = row['gene']
        gene_2 = row['gene_2']
        id = row['id']
        
        if pd.notna(gene) and pd.isna(gene_2):
            castor_gene_info[castor_gene] = ('1', id, gene)
        elif pd.isna(gene) and pd.notna(gene_2):
            castor_gene_info[castor_gene] = ('2', id, gene_2)
        elif pd.notna(gene) and pd.notna(gene_2):
            castor_gene_info[castor_gene] = ('B', id, gene, gene_2)
        else:
            castor_gene_info[castor_gene] = ('None',)

    # 遍历文件B，根据字典信息补充数据
    for index, row in file_b.iterrows():
        castor_gene = row[1]
        if castor_gene in castor_gene_info:
            info = castor_gene_info[castor_gene]
            file_b.at[index, 7] = ' '.join(map(str, info))
        else:
            file_b.at[index, 7] = 'None'

    # 保存修改后的文件B
    file_b.to_csv(output_path, sep=',', header=False, index=False)

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python process_genes.py <fileA.csv> <fileB.txt> <output.txt>")
        sys.exit(1)

    file_a_path = sys.argv[1]
    file_b_path = sys.argv[2]
    output_path = sys.argv[3]

    process_files(file_a_path, file_b_path, output_path)
