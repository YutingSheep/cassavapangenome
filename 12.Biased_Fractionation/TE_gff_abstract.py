import sys

def load_regions(fileB_path):
    """加载文件B的区域信息，返回字典：{chr: [(start1, end1), (start2, end2), ...]}"""
    regions = {}
    with open(fileB_path, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) < 3:
                continue
            chr_b, start_b, end_b = parts[0], int(parts[1]), int(parts[2])
            if chr_b not in regions:
                regions[chr_b] = []
            regions[chr_b].append( (start_b, end_b) )
    return regions

def filter_gff(fileA_path, regions):
    """过滤GFF文件，输出完全包含在区域内的条目"""
    with open(fileA_path, 'r') if fileA_path != '-' else sys.stdin as f:
        for line in f:
            if line.startswith('#'):
                # 保留注释行（可选，按需注释掉）
                # print(line.strip())
                continue
            parts = line.strip().split('\t')
            if len(parts) < 5:
                continue
            chr_a = parts[0]
            start_a = int(parts[3])
            end_a = int(parts[4])
            
            # 检查染色体是否存在区域
            if chr_a in regions:
                for (start_b, end_b) in regions[chr_a]:
                    if start_a > start_b and end_a < end_b:
                        print(line.strip())
                        break  # 匹配任一区域即可

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print(f"Usage: {sys.argv[0]} <fileB.txt> <fileA.gff>")
        print("       Use '-' for fileA.gff to read from stdin")
        sys.exit(1)
    
    fileB_path = sys.argv[1]
    fileA_path = sys.argv[2]
    
    regions = load_regions(fileB_path)
    filter_gff(fileA_path, regions)
