import sys

def process_blocks(input_file):
    blocks = []
    current_block = None
    
    with open(input_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            
            # 解析rate和gene
            parts = line.split('\t')
            if len(parts) < 2:
                continue
            
            try:
                rate = float(parts[0])
                gene = parts[1]
            except:
                continue
            
            # 判断是否进入块条件
            if abs(rate) > 0.05:
                current_sign = 'pos' if rate > 0 else 'neg'
                
                if current_block:
                    # 检查符号连续性
                    if current_block['sign'] == current_sign:
                        current_block['sum'] += rate
                        current_block['count'] += 1
                        current_block['end'] = gene
                    else:
                        # 保存旧块，开始新块
                        blocks.append(current_block)
                        current_block = {
                            'sum': rate,
                            'count': 1,
                            'start': gene,
                            'end': gene,
                            'sign': current_sign
                        }
                else:
                    # 开始新块
                    current_block = {
                        'sum': rate,
                        'count': 1,
                        'start': gene,
                        'end': gene,
                        'sign': current_sign
                    }
            else:
                if current_block:
                    blocks.append(current_block)
                    current_block = None
    
    # 添加最后一个块
    if current_block:
        blocks.append(current_block)
    
    # 格式化为输出结果
    output = []
    for block in blocks:
        avg = block['sum'] / block['count']
        output.append(f"{avg:.4f}\t{block['start']}\t{block['end']}")
    
    return output

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python gene_blocks.py input.txt")
        sys.exit(1)
    
    results = process_blocks(sys.argv[1])
    print("mean_rate\tstart_gene\tend_gene")
    for line in results:
        print(line)
