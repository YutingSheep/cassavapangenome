import random
import argparse
from collections import defaultdict

def parse_arguments():
    parser = argparse.ArgumentParser(description='自交模拟工具 v3.1')
    parser.add_argument('-i', '--input', required=True, help='输入文件路径')
    parser.add_argument('-o', '--output', required=True, help='输出文件路径')
    parser.add_argument('-s', '--seed', type=int, help='随机数种子')
    return parser.parse_args()

def parse_haplotypes(header):
    """正确解析样本结构：<序号>_<样本名>_hap<数字>"""
    hap_dict = defaultdict(lambda: defaultdict(dict))
    for idx, col in enumerate(header[2:]):  # 跳过CHROM和POS列
        parts = col.split('_')
        if len(parts) < 3: continue
        
        # 解析关键元素
        sample_id = parts[1]  # 第二个元素为样本名
        hap_type = parts[2]   # 第三个元素为单倍型
        
        hap_dict[sample_id][hap_type] = idx + 2  # 存储列索引
    return hap_dict

def simulate_selfing(input_file, output_file, seed=None):
    if seed: random.seed(seed)
    
    # 读取并预处理数据
    with open(input_file, 'r') as f:
        lines = [line.strip() for line in f]
    
    header = lines[0].split('\t')
    data_lines = lines[1:]
    chromosomes = sorted({line.split('\t')[0] for line in data_lines})
    
    # 解析样本结构
    hap_dict = parse_haplotypes(header)
    valid_samples = [s for s in hap_dict if len(hap_dict[s]) >= 2]  # 至少有双单倍型
    
    if not valid_samples:
        raise ValueError("没有符合要求的样本（需至少包含两个单倍型）")
    selected_sample = random.choice(valid_samples)
    
    # 构建染色体选择方案
    chrom_selections = {}
    for chrom in chromosomes:
        # 获取该样本所有可用单倍型（假设单倍型在染色体层面存在）
        available_haps = list(hap_dict[selected_sample].keys())
        hap1 = random.choice(available_haps)
        hap2 = random.choice(available_haps)  # 允许重复选择
        
        chrom_selections[chrom] = (
            hap_dict[selected_sample][hap1],
            hap_dict[selected_sample][hap2]
        )
    
    # 写入输出文件
    with open(output_file, 'w') as f:
        f.write(f"## Selfing Sample: {selected_sample}\n")
        for chrom in chromosomes:
            h1 = header[chrom_selections[chrom][0]].split('_')[-1]
            h2 = header[chrom_selections[chrom][1]].split('_')[-1]
            f.write(f"## {chrom}_selection: {h1}+{h2}\n")
        
        # 写入数据头
        f.write("\t".join(["#CHROM", "POS", "Progeny_hap1", "Progeny_hap2"]) + "\n")
        
        # 处理每个位点
        for line in data_lines:
            fields = line.split('\t')
            chrom = fields[0]
            cols = chrom_selections[chrom]
            
            a1 = fields[cols[0]].split('/')[0]
            a2 = fields[cols[1]].split('/')[0]
            f.write(f"{chrom}\t{fields[1]}\t{a1}/{a2}\n")

if __name__ == "__main__":
    args = parse_arguments()
    try:
        simulate_selfing(args.input, args.output, args.seed)
        print(f"自交模拟完成！结果已保存至 {args.output}")
    except Exception as e:
        print(f"错误: {str(e)}")

