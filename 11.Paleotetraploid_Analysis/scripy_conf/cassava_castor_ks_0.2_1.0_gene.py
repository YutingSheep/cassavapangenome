import sys

# 确保提供了正确数量的命令行参数
if len(sys.argv) != 4:
    print("Usage: python script.py <input_file_A> <input_file_B> <output_file>")
    sys.exit(1)

# 获取输入文件和输出文件路径
file_a_path = sys.argv[1]
file_b_path = sys.argv[2]
output_path = sys.argv[3]

# 构建基因编号到基因名称的映射
gene_mapping = {}
with open(file_a_path, 'r') as file_a:
    for line in file_a:
        parts = line.strip().split()
        species1_gene = parts[1]
        species1_id = parts[2]
        species2_gene = parts[4]
        species2_id = parts[5]
        
        # 构建Species1的基因编号到基因名称的映射
        gene_mapping[(parts[0], species1_id)] = species1_gene
        # 构建Species2的基因编号到基因名称的映射
        gene_mapping[(parts[3], species2_id)] = species2_gene

# 处理文件B，生成结果
with open(file_b_path, 'r') as file_b, open(output_path, 'w') as output:
    for line in file_b:
        parts = line.strip().split(',')
        id_ = parts[0]
        chr1 = parts[1]
        chr2 = parts[2]
        block1 = parts[3].split('_')
        block2 = parts[4].split('_')
        
        # 遍历基因编号对
        for gene1_id, gene2_id in zip(block1, block2):
            # 获取基因名称
            gene1_name = gene_mapping.get((chr1, gene1_id), 'Unknown')
            gene2_name = gene_mapping.get((chr2, gene2_id), 'Unknown')
            
            # 写入结果文件
            output.write(f"{id_}\t{gene1_name}\t{gene2_name}\n")
