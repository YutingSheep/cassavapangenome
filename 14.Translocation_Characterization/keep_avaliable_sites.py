import argparse

def filter_vcf(input_file, output_file):
    with open(input_file, 'r') as vcf_in, open(output_file, 'w') as vcf_out:
        for line in vcf_in:
            # 保留 VCF 的头文件部分
            if line.startswith("#"):
                vcf_out.write(line)
                continue

            # 分割每行的字段
            fields = line.strip().split('\t')
            genotypes = fields[9:]  # 样本的基因型从第10列开始

            # 检查基因型，统计 ./., 0/0 和 1/1 的数量
            gt_counts = {"0/0": 0, "1/1": 0, "./.": 0}
            total_gt = 0

            for gt in genotypes:
                gt_value = gt.split(":")[0]  # 取出基因型信息 (GT 位于冒号前)
                if gt_value in gt_counts:
                    gt_counts[gt_value] += 1
                total_gt += 1

            # 如果位点只含有 ./., 0/0 或者只含有 ./., 1/1，则跳过
            if (gt_counts["0/0"] + gt_counts["./."] == total_gt) or (gt_counts["1/1"] + gt_counts["./."] == total_gt):
                continue  # 跳过该行

            # 否则，保留该位点
            vcf_out.write(line)

if __name__ == "__main__":
    # 解析命令行参数
    parser = argparse.ArgumentParser(description="Filter VCF file based on genotype conditions.")
    parser.add_argument("-i", "--input", required=True, help="Input VCF file")
    parser.add_argument("-o", "--output", required=True, help="Output filtered VCF file")
    
    args = parser.parse_args()

    # 调用过滤函数
    filter_vcf(args.input, args.output)

