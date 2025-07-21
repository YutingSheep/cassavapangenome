import argparse

def main():
    parser = argparse.ArgumentParser(description='Merge genomic blocks based on VARIANTS/KB threshold.')
    parser.add_argument('-i', '--input', required=True, help='Input file path')
    parser.add_argument('-t', '--threshold', type=float, required=True, help='Threshold value for merging')
    args = parser.parse_args()

    # 读取数据，跳过标题行
    data = []
    with open(args.input, 'r') as f:
        headers = next(f)  # 跳过标题行
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 4:
                continue  # 跳过不完整的行
            chrom = parts[0]
            bin_start = int(parts[1])
            bin_end = int(parts[2])
            vpkb = float(parts[3])  # 第五列是VARIANTS/KB
            data.append({
                'chrom': chrom,
                'start': bin_start,
                'end': bin_end,
                'vpkb': vpkb
            })

    # 按染色体分组
    grouped = {}
    for entry in data:
        chrom = entry['chrom']
        if chrom not in grouped:
            grouped[chrom] = []
        grouped[chrom].append(entry)

    # 处理每个染色体，合并块
    merged = []
    for chrom in grouped:
        entries = grouped[chrom]
        current_block = None
        for entry in entries:
            current_vpkb = entry['vpkb']
            current_condition = current_vpkb > args.threshold
            if current_block is None:
                current_block = {
                    'chrom': chrom,
                    'start': entry['start'],
                    'end': entry['end'],
                    'condition': current_condition
                }
            else:
                if current_block['chrom'] == chrom and current_condition == current_block['condition']:
                    current_block['end'] = entry['end']
                else:
                    merged.append(current_block)
                    current_block = {
                        'chrom': chrom,
                        'start': entry['start'],
                        'end': entry['end'],
                        'condition': current_condition
                    }
        if current_block is not None:
            merged.append(current_block)

    # 输出结果（添加条件状态）
    print("CHROM\tBIN_START\tBIN_END\tCONDITION")
    for block in merged:
        condition = ">%.2f" % args.threshold if block['condition'] else "<=%.2f" % args.threshold
        print(f"{block['chrom']}\t{block['start']}\t{block['end']}\t{condition}")

if __name__ == '__main__':
    main()

