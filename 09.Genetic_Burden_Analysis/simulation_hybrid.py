import random
import argparse
from collections import defaultdict

def parse_arguments():
    parser = argparse.ArgumentParser(description='染色体级遗传模拟工具 v2.1')
    parser.add_argument('-i', '--input', required=True, help='输入文件路径')
    parser.add_argument('-o', '--output', required=True, help='输出文件路径')
    parser.add_argument('-s', '--seed', type=int, help='随机数种子')
    return parser.parse_args()

def parse_haplotypes(header):
    """正确解析样本名称（两个下划线中间部分）"""
    hap_dict = defaultdict(lambda: defaultdict(dict))
    for idx, col in enumerate(header[2:]):
        parts = col.split('_')
        if len(parts) < 3: continue
        sample_id = parts[1]  # 第二个元素为样本名称
        hap_type = parts[2]
        hap_dict[sample_id][hap_type] = idx + 2  # 存储列索引
    return hap_dict

def simulate_progeny(input_file, output_file, seed=None):
    if seed: random.seed(seed)
    
    with open(input_file, 'r') as f:
        lines = [line.strip() for line in f]
    
    header = lines[0].split('\t')
    data_lines = lines[1:]
    
    # 解析染色体结构
    chromosomes = sorted({line.split('\t')[0] for line in data_lines})
    hap_dict = parse_haplotypes(header)
    
    # 选择亲本（基于正确解析的样本名称）
    samples = [s for s in hap_dict.keys() if len(hap_dict[s]) == 2]  # 过滤有效样本
    if len(samples) < 2: 
        raise ValueError("需要至少两个具有双单倍型的样本")
    parentA, parentB = random.sample(samples, 2)
    
    # 为每个染色体独立选择单倍型
    chrom_map = {}
    for chrom in chromosomes:
        # 从每个亲本的两个单倍型中随机选择
        hapA = random.choice(['hap1', 'hap2'])
        hapB = random.choice(['hap1', 'hap2'])
        
        # 验证单倍型存在性
        if hapA not in hap_dict[parentA] or hapB not in hap_dict[parentB]:
            raise KeyError("单倍型数据不完整")
        
        chrom_map[chrom] = {
            'cols': (hap_dict[parentA][hapA], hap_dict[parentB][hapB]),
            'haps': (hapA, hapB)
        }
    
    # 写入输出文件
    with open(output_file, 'w') as f:
        f.write(f"## Global Parents: {parentA} × {parentB}\n")
        for chrom in chromosomes:
            hapA, hapB = chrom_map[chrom]['haps']
            f.write(f"## {chrom}_selection: {parentA}_{hapA} + {parentB}_{hapB}\n")
        
        f.write("\t".join(["#CHROM", "POS", "Progeny_hap1", "Progeny_hap2"]) + "\n")
        
        for line in data_lines:
            fields = line.split('\t')
            chrom = fields[0]
            pos = fields[1]
            cols = chrom_map[chrom]['cols']
            
            a1 = fields[cols[0]].split('/')[0]
            a2 = fields[cols[1]].split('/')[0]
            f.write(f"{chrom}\t{pos}\t{a1}/{a2}\n")

if __name__ == "__main__":
    args = parse_arguments()
    try:
        simulate_progeny(args.input, args.output, args.seed)
        print(f"模拟完成！染色体选择方案已更新")
    except Exception as e:
        print(f"错误: {str(e)}")

