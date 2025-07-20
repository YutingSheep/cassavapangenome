def extract_middle_gene(input_file, output_file=None):
    """
    从文件中读取基因对，计算每行的中间基因，并验证染色体一致性。
    """
    try:
        with open(input_file, 'r') as f:
            lines = f.readlines()
    except FileNotFoundError:
        print(f"错误：文件 '{input_file}' 不存在")
        return

    results = []
    for line_num, line in enumerate(lines, 1):
        line = line.strip()
        if not line:
            continue

        try:
            # 分割基因对
            parts = line.split()
            if len(parts) != 2:
                raise ValueError("每行必须包含两个基因编号")

            gene1, gene2 = parts

            # 提取染色体前缀和编号
            prefix1, num_part1 = gene1.split('g', 1)  # 分割一次，防止前缀含'g'
            prefix2, num_part2 = gene2.split('g', 1)
            num1 = int(num_part1)
            num2 = int(num_part2)

            # 检查染色体一致性
            if prefix1 != prefix2:
                raise ValueError(f"染色体前缀不一致: '{prefix1}' vs '{prefix2}'")

            # 检查数字顺序
            if num1 >= num2:
                raise ValueError(f"起始编号 {num1} ≥ 结束编号 {num2}")

            # 计算中间基因
            mid = round((num1 + num2) / 2)
            middle_gene = f"{prefix1}g{mid:05d}"
            results.append(middle_gene)

        except (ValueError, IndexError) as e:
            print(f"第 {line_num} 行错误: {e}，跳过处理")
            continue

    # 输出结果
    if output_file:
        with open(output_file, 'w') as f:
            f.write('\n'.join(results))
        print(f"结果已保存至 {output_file}（共 {len(results)} 个中间基因）")
    else:
        print('\n'.join(results))


# 示例调用
extract_middle_gene(
    input_file='AM560.region',
    output_file='middle_genes.txt'
)

