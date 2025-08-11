import sys

def main(file_a_path, file_b_path):
    id_to_line = {}
    lines_a = []
    # 解析文件A并建立ID映射
    with open(file_a_path, 'r') as f_a:
        for line_num, line in enumerate(f_a):
            line = line.strip()
            if not line:
                continue
            parts = line.split(',')
            first_id = parts[0]
            id_to_line[first_id] = line_num  # 记录ID与行号的映射
            lines_a.append(line)  # 保存完整行内容 
    
    # 解析文件B并处理序号
    with open(file_b_path, 'r') as f_b:
        for line in f_b:
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) < 3:  # 新增序号后需至少3列 
                continue
            seq_num = parts[0]  # 提取文件B的序号
            start_id, end_id = parts[1], parts[2]  # 调整ID索引位置 
            
            # 获取起止行号
            start_line = id_to_line.get(start_id)
            end_line = id_to_line.get(end_id)
            if start_line is None or end_line is None:
                continue
            
            # 输出带序号的内容
            start = min(start_line, end_line)
            end = max(start_line, end_line)
            for i in range(start, end + 1):
                # 格式：文件B序号 + 文件A原始行（逗号分隔）
                print(f"{seq_num},{lines_a[i]}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("用法：python script.py <文件A路径> <文件B路径>")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])

