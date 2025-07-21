import argparse
import pysam
import subprocess

def main(input_vcf, output_file, bed_file):
    # 获取基因组窗口行数并计算len_genome
    command = f"cat {bed_file} | wc -l"
    result = subprocess.getoutput(command)  # 获取wc -l的返回值
    num_lines = int(result.strip())  # 转换为整数
    len_genome = num_lines * 100000  # 基因组长度

    # 打开VCF文件
    vcf_in = pysam.VariantFile(input_vcf)

    # 初始化样本SNP计数的字典
    sample_snp_counts = {sample: 0 for sample in vcf_in.header.samples}

    # 遍历VCF文件中的每个记录（变异位点）
    for record in vcf_in.fetch():
        ref_allele = record.ref
        alt_alleles = record.alts

        # 只考虑SNP（忽略Indel）
        if len(ref_allele) == 1 and all(len(alt) == 1 for alt in alt_alleles):
            # 统计当前位点中各样本的0/0和1/1的数量
            count_00 = 0
            count_11 = 0

            for sample in vcf_in.header.samples:
                sample_call = record.samples[sample]
                genotype = sample_call['GT']  # 获取样本的基因型

                # 跳过./.（缺失数据）
                if genotype is None or genotype == (None, None) or genotype == './.':
                    continue

                # 统计0/0和1/1的数量
                if genotype == (0, 0):
                    count_00 += 1
                elif genotype == (1, 1):
                    count_11 += 1

            # 确定祖先状态：1/1少于或等于0/0，0/0为祖先；否则，1/1为祖先
            ancestral_state = (0, 0) if count_00 >= count_11 else (1, 1)

            # 遍历每个样本，根据祖先状态统计SNP
            for sample in vcf_in.header.samples:
                sample_call = record.samples[sample]
                genotype = sample_call['GT']

                # 跳过./.（缺失数据）
                if genotype is None or genotype == (None, None) or genotype == './.':
                    continue

                # 如果祖先状态是0/0，则直接统计不同于0/0的基因型
                if ancestral_state == (0, 0) and genotype != (0, 0):
                    sample_snp_counts[sample] += 1
                # 如果祖先状态是1/1，则统计不同于1/1的基因型，并且翻转0/0和1/1
                elif ancestral_state == (1, 1):
                    if genotype == (0, 0):
                        sample_snp_counts[sample] += 1  # 0/0被视为突变
                    elif genotype == (1, 1):
                        continue  # 1/1为祖先状态，不计为突变

    # 假设的突变率和每年的繁殖代数
    mutation_rate = 2.7e-8  # 每代的突变率
    generations_per_year = 1  # 每年繁殖一次

    # 输出结果到文件，格式为三列
    with open(output_file, 'w') as out_f:
        out_f.write("Sample\tSNPs\tEstimated_Asexual_Reproduction_Years\n")  # 表头
        for sample, snp_count in sample_snp_counts.items():
            time_estimate = snp_count / (len_genome * mutation_rate * generations_per_year)
            out_f.write(f"{sample}\t{snp_count}\t{time_estimate}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Estimate asexual reproduction time based on VCF.")
    parser.add_argument('-i', '--input', required=True, help='Input VCF file')
    parser.add_argument('-o', '--output', required=True, help='Output file to write results')
    parser.add_argument('-b', '--bed', required=True, help='BED file for genome length calculation')

    args = parser.parse_args()
    
    main(args.input, args.output, args.bed)

