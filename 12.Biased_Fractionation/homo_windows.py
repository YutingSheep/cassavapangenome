import pandas as pd
import argparse

def calculate_homo_proportions(window, sample_columns):
    """计算每个样本的 homo1、homo2 和 homo_B 的占比"""
    results = {}
    for sample in sample_columns:
        homo1 = ((window[sample] == '1') | (window[sample] == 'B')).sum() / len(window)
        homo2 = ((window[sample] == '2') | (window[sample] == 'B')).sum() / len(window)
        homo_B = (window[sample] == 'B').sum() / len(window)
        results[f"{sample}_homo1"] = homo1
        results[f"{sample}_homo2"] = homo2
        results[f"{sample}_homo_B"] = homo_B
    return results

def process_chromosome(chromosome_data, sample_columns, window_size, step_size):
    """处理单个染色体的数据"""
    results = []
    for start in range(0, len(chromosome_data) - window_size + 1, step_size):
        end = start + window_size
        window = chromosome_data.iloc[start:end]

        # 获取基因窗口范围
        gene_start = window.iloc[0]['Key']
        gene_end = window.iloc[-1]['Key']

        # 计算每个样本的 homo1、homo2 和 homo_B 的占比
        homo_results = calculate_homo_proportions(window, sample_columns)

        # 添加到结果列表
        results.append([chromosome_data.iloc[0]['chr'], f"{gene_start}-{gene_end}"] + list(homo_results.values()))
    return results

def process_file(input_file, output_file, window_size, step_size):
    """处理输入文件并生成输出文件"""
    # 读取数据
    df = pd.read_csv(input_file)

    # 获取样本列（排除固定的前两列）
    sample_columns = df.columns[2:]

    # 初始化结果列表
    all_results = []

    # 按染色体分组处理
    for chromosome, chromosome_data in df.groupby('chr'):
        chromosome_results = process_chromosome(chromosome_data, sample_columns, window_size, step_size)
        all_results.extend(chromosome_results)

    # 创建结果 DataFrame
    result_columns = ['chr', 'gene_window']
    for sample in sample_columns:
        result_columns.extend([f"{sample}_homo1", f"{sample}_homo2", f"{sample}_homo_B"])
    result_df = pd.DataFrame(all_results, columns=result_columns)

    # 保存结果到 CSV 文件
    result_df.to_csv(output_file, index=False)
    print(f"转换完成，结果已保存到 {output_file}")

def main():
    # 设置命令行参数
    parser = argparse.ArgumentParser(description="将基因数据转换为滑动窗口格式")
    parser.add_argument("input_file", help="输入文件路径")
    parser.add_argument("output_file", help="输出文件路径")
    parser.add_argument("--window_size", type=int, default=100, help="窗口大小（默认：100）")
    parser.add_argument("--step_size", type=int, default=10, help="滑动步长（默认：10）")
    args = parser.parse_args()

    # 处理文件
    process_file(args.input_file, args.output_file, args.window_size, args.step_size)

if __name__ == "__main__":
    main()
