import sys
import bisect
from collections import defaultdict

def parse_gene(gene_id):
    """解析染色体和基因位置"""
    try:
        chrom_part, num_part = gene_id.split('g')
        return chrom_part, int(num_part.lstrip('0') or 0)  # 处理前导零
    except ValueError:
        raise ValueError(f"Invalid gene ID format: {gene_id}")

def build_chrom_index(ref_file):
    """构建染色体位置索引"""
    chrom_index = defaultdict(lambda: {'positions': [], 'genes': []})
    
    with open(ref_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            gene = line.split('\t')[0] if '\t' in line else line
            try:
                chrom, pos = parse_gene(gene)
                chrom_index[chrom]['positions'].append(pos)
                chrom_index[chrom]['genes'].append(gene)
            except ValueError:
                continue
    
    # 对每个染色体的位置排序
    for chrom in chrom_index:
        sorted_pairs = sorted(zip(chrom_index[chrom]['positions'], 
                                chrom_index[chrom]['genes']))
        positions, genes = zip(*sorted_pairs) if sorted_pairs else ([], [])
        chrom_index[chrom]['positions'] = list(positions)
        chrom_index[chrom]['genes'] = list(genes)
    
    return chrom_index

def find_nearest(query_file, chrom_index):
    """执行最近基因查询（仅查找前方基因）"""
    results = []
    
    with open(query_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            query_gene = line.split('\t')[0] if '\t' in line else line
            try:
                q_chrom, q_pos = parse_gene(query_gene)
            except ValueError:
                results.append(f"{query_gene}\tNA\tNA\tNA")
                continue
            
            if q_chrom not in chrom_index:
                results.append(f"{query_gene}\tNA\t{q_chrom}\tNA")
                continue
                
            positions = chrom_index[q_chrom]['positions']
            genes = chrom_index[q_chrom]['genes']
            
            # 使用bisect_right获取前方基因
            idx = bisect.bisect_right(positions, q_pos)
            
            if idx > 0:
                nearest_pos = positions[idx-1]
                nearest_gene = genes[idx-1]
                distance = q_pos - nearest_pos  # 保证非负数
                results.append(f"{query_gene}\t{nearest_gene}\t{q_chrom}\t{distance}")
            else:
                results.append(f"{query_gene}\tNA\t{q_chrom}\tNA")
    
    return results


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python nearest_gene.py query.txt ref.txt")
        sys.exit(1)
    
    # 构建索引
    chrom_index = build_chrom_index(sys.argv[2])
    
    # 执行查询
    output = find_nearest(sys.argv[1], chrom_index)
    
    # 输出结果
    print("QueryGene\tNearestgri\tChromosome\tDistance")
    for line in output:
        print(line)
