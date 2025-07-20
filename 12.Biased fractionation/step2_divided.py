import sys

def process_file(input_file):
    current_block = 1
    # 新增第一列的跟踪变量
    prev_col0 = None
    prev_col2 = None
    prev_col5 = None
    prev_col8 = None  
    prev_col3 = None
    prev_col6 = None
    prev_col9 = None
    
    with open(input_file, 'r') as f:
        for line in f:
            line = line.strip()
            cols = line.split('\t')
            # 检查列数是否足够（新增一列后需要至少10列）
            if len(cols) < 10:
                print(f"{line}\tblock{current_block}")
                continue
            
            # 新增首列处理（索引0）
            col0 = cols[0]
            # 调整其他列的索引（原索引+1）
            col2 = cols[2]   # 原第二列现在在索引2
            col5 = cols[5]   # 原第五列现在在索引5
            col8 = cols[8]   # 原第八列现在在索引8
            
            try:
                # 调整数值列的索引
                col3 = int(cols[3])  # 原第三列现在在索引3
                col6 = int(cols[6])  # 原第六列现在在索引6
                col9 = int(cols[9])  # 原第九列现在在索引9
            except (ValueError, IndexError) as e:
                print(f"{line}\tblock{current_block}")
                continue
            
            # 首行初始化（包含新增列）
            if all(v is None for v in [prev_col0, prev_col2, prev_col5, prev_col8]):
                prev_col0, prev_col2, prev_col5, prev_col8 = col0, col2, col5, col8
                prev_col3, prev_col6, prev_col9 = col3, col6, col9
                print(f"{line}\tblock{current_block}")
                continue
            
            # 新增首列必须相同才进行其他判断
            same_columns = (col0 == prev_col0) and (col2 == prev_col2) and (col5 == prev_col5) and (col8 == prev_col8)
            delta3 = abs(col3 - prev_col3)
            delta6 = abs(col6 - prev_col6)
            delta9 = abs(col9 - prev_col9)
            
            # 条件判断逻辑
            if same_columns and delta3 <= 20 and delta6 <= 20 and delta9 <= 20:
                block = current_block
            else:
                current_block += 1
                block = current_block
            
            # 更新所有跟踪变量
            prev_col0, prev_col2, prev_col5, prev_col8 = col0, col2, col5, col8
            prev_col3, prev_col6, prev_col9 = col3, col6, col9
            
            print(f"{line}\tblock{block}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python3 block_annotation.py input.txt")
        sys.exit(1)
    process_file(sys.argv[1])

