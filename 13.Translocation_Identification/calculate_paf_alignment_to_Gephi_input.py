import csv
import argparse

# 创建参数解析器
parser = argparse.ArgumentParser(description="Process alignment data.")
parser.add_argument('-i', '--input', required=True, help="Input file name")
parser.add_argument('-o', '--output', required=True, help="Output file name")

# 解析参数
args = parser.parse_args()
input_file = args.input
output_file = args.output

# 初始化变量
current_ref = None
current_query = None
total_match_length = 0

# 打开输入和输出文件
with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
    reader = csv.reader(infile, delimiter='\t')
    writer = csv.writer(outfile, delimiter='\t')

    # 写入输出文件的表头
    writer.writerow(['ref', 'query', 'ref_length', 'query_length', 'total_match_length'])

    for row in reader:
        # 解析每一行
        ref = row[0]
        ref_length = int(row[1])
        ref_start = int(row[2])
        ref_end = int(row[3])
        query = row[5]
        query_length = int(row[6])

        # 计算当前行的匹配长度
        if float(row[12].split(':')[-1]) > 0.95 and float(row[13].split(':')[-1]) > 0.95:
            match_length = int(row[10]) * float(row[12].split(':')[-1])
        else:
            match_length = 0

        # 如果遇到新的ref或query，输出之前累积的结果
        if current_ref != ref or current_query != query:
            if current_ref is not None and current_query is not None:
                if current_ref_length > 1000000 and current_query_length > 1000000 and (200 * total_match_length) / (current_ref_length + current_query_length) > 10:
                    writer.writerow([current_ref, current_query, current_ref_length, current_query_length, (200 * total_match_length) / (current_ref_length + current_query_length)])

            # 更新当前ref和query
            current_ref = ref
            current_query = query
            current_ref_length = ref_length
            current_query_length = query_length
            total_match_length = 0

        # 累加匹配长度
        total_match_length += match_length

    # 写入最后一个ref-query对的结果
    if current_ref is not None and current_query is not None:
        if current_ref_length > 1000000 and current_query_length > 1000000 and (200 * total_match_length) / (current_ref_length + current_query_length) > 10:
            writer.writerow([current_ref, current_query, current_ref_length, current_query_length, (200 * total_match_length) / (current_ref_length + current_query_length)])

print(f"Analysis complete. Results saved to {output_file}.")

