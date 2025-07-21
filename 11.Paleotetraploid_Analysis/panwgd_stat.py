import argparse

# 统计 Left 和 Right 基因数目
def count_genes(row):
    left_count = 0
    right_count = 0
    for i, value in enumerate(row):
        if "_Left" in header[i] and value != "None":
            left_count += 1
        if "_Right" in header[i] and value != "None":
            right_count += 1
    return left_count, right_count

# 根据 Left_Count 和 Right_Count 标记
def mark_genes(left_count, right_count):
    if left_count == 117 and right_count == 117:
        return "B"
    elif left_count == 0 and right_count == 117:
        return "2"
    elif left_count == 117 and right_count == 0:
        return "1"
    elif left_count == 0 and right_count == 0:
        return "None"
    else:
        return "V"

# 主函数
def main():
    # 设置命令行参数
    parser = argparse.ArgumentParser(description="统计结果文件中每个键的 Left 和 Right 基因数目，并标记。")
    parser.add_argument("input", help="输入文件路径（结果文件 C）")
    parser.add_argument("output", help="输出文件路径")
    args = parser.parse_args()

    # 读取输入文件
    with open(args.input, 'r') as file:
        lines = file.readlines()

    # 解析标题行
    global header
    header = lines[0].strip().split(',')

    # 统计并写入输出文件
    with open(args.output, 'w') as file:
        file.write("Key,Left_Count,Right_Count,Marker\n")  # 写入标题行
        for line in lines[1:]:
            row = line.strip().split(',')
            key = row[0]  # 第一列是键
            left_count, right_count = count_genes(row)
            marker = mark_genes(left_count, right_count)
            file.write(f"{key},{left_count},{right_count},{marker}\n")

    print(f"统计完成，结果已保存到 {args.output}")

if __name__ == "__main__":
    main()
