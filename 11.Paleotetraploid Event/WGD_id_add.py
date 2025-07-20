import pandas as pd
import argparse

def parse_arguments():
    parser = argparse.ArgumentParser(description="Add block IDs to genes from two input files.")
    parser.add_argument("-a", "--file_a", required=True, help="Path to file A")
    parser.add_argument("-b", "--file_b", required=True, help="Path to file B")
    parser.add_argument("-o", "--output", required=True, help="Path to the output file")
    return parser.parse_args()

def find_ids(row, file_b):
    cassava_chromosome = row["Cassava_Chromosome"]
    cassava_id = row["Cassava_Gene_Number"]
    castor_chromosome = row["Castor_Chromosome"]
    castor_id = row["Castor_Gene_Number"]
    
    # 匹配木薯基因编号和染色体
    cassava_match = file_b[
        (file_b["Cassava_Chromosome"] == cassava_chromosome) &
        (file_b["Cassava_Gene_Numbers"].apply(lambda x: cassava_id in x))
    ]
    # 匹配蓖麻基因编号和染色体
    castor_match = file_b[
        (file_b["Castor_Chromosome"] == castor_chromosome) &
        (file_b["Castor_Gene_Numbers"].apply(lambda x: castor_id in x))
    ]
    
    # 收集匹配到的ID，用逗号分隔，找不到则为UNID
    cassava_block_id = ",".join(map(str, cassava_match["ID"].tolist())) if not cassava_match.empty else "UNID"
    castor_block_id = ",".join(map(str, castor_match["ID"].tolist())) if not castor_match.empty else "UNID"
    return pd.Series([cassava_block_id, castor_block_id])

def main():
    args = parse_arguments()
    
    # 加载文件A
    file_a = pd.read_csv(args.file_a, sep="\t", header=None)
    file_a.columns = ["Cassava_Gene_Name", "Cassava_Gene_Number", "Castor_Gene_Name", "Castor_Gene_Number"]
    
    # 提取染色体信息
    file_a["Cassava_Chromosome"] = file_a["Cassava_Gene_Name"].str.extract(r"md(.*?)g")
    file_a["Castor_Chromosome"] = file_a["Castor_Gene_Name"].str.extract(r"md(.*?)g")
    
    # 加载文件B
    file_b = pd.read_csv(args.file_b, sep=",", header=None)
    print(f"File B shape: {file_b.shape}")  # 调试信息
    # 动态生成列名
    default_columns = ["ID", "Cassava_Chromosome", "Castor_Chromosome", "Cassava_Gene_Ids", "Castor_Gene_Ids"]
    file_b.columns = default_columns[:file_b.shape[1]]  # 根据列数自动截取

    # 将基因编号字符串分割为列表，并跳过非数值数据
    if "Cassava_Gene_Ids" in file_b.columns:
        file_b["Cassava_Gene_Numbers"] = file_b["Cassava_Gene_Ids"].str.split("_").apply(
            lambda x: [int(i) for i in x if i.isdigit()]
        )
    if "Castor_Gene_Ids" in file_b.columns:
        file_b["Castor_Gene_Numbers"] = file_b["Castor_Gene_Ids"].str.split("_").apply(
            lambda x: [int(i) for i in x if i.isdigit()]
        )
    
    # 匹配ID并添加到文件A
    file_a[["Cassava_Block_ID", "Castor_Block_ID"]] = file_a.apply(find_ids, axis=1, file_b=file_b)
    
    # 保存更新后的文件A
    file_a.to_csv(args.output, sep="\t", index=False)

if __name__ == "__main__":
    main()
