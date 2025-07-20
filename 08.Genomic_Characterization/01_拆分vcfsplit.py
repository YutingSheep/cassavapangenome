#!/usr/bin/env python3
import os
import sys
import argparse
import subprocess

# 预定义的染色体长度信息
CHROM_LENGTHS = {
    'chr01': 43843061,
    'chr02': 39427047,
    'chr03': 33962347,
    'chr04': 36653984,
    'chr05': 34044637,
    'chr06': 31332623,
    'chr07': 36158143,
    'chr08': 42444776,
    'chr09': 39208591,
    'chr10': 32469860,
    'chr11': 33706567,
    'chr12': 37647568,
    'chr13': 37638087,
    'chr14': 29694991,
    'chr15': 33308024,
    'chr16': 34389555,
    'chr17': 37185621,
    'chr18': 33483791
}

def get_vcf_boundaries(vcf_file):
    """获取VCF文件中实际存在的染色体和位置范围"""
    chrom_ranges = {}
    last_chrom = None
    last_pos = 0
    
    with open(vcf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            parts = line.strip().split('\t')
            chrom = parts[0]
            pos = int(parts[1])
            
            if chrom not in chrom_ranges:
                chrom_ranges[chrom] = {'min': pos, 'max': pos}
            else:
                if pos < chrom_ranges[chrom]['min']:
                    chrom_ranges[chrom]['min'] = pos
                if pos > chrom_ranges[chrom]['max']:
                    chrom_ranges[chrom]['max'] = pos
    
    return chrom_ranges

def process_vcf(input_vcf, output_dir, window_size=100000):
    """使用100kb滑窗处理VCF文件"""
    # 确保输出目录存在
    os.makedirs(output_dir, exist_ok=True)
    
    # 获取VCF文件中实际存在的染色体和位置范围
    chrom_ranges = get_vcf_boundaries(input_vcf)
    print("检测到的染色体范围:")
    for chrom, ranges in chrom_ranges.items():
        print(f"{chrom}: {ranges['min']}-{ranges['max']}")
    
    # 对每条染色体进行处理
    for chrom, length in CHROM_LENGTHS.items():
        if chrom not in chrom_ranges:
            print(f"跳过染色体 {chrom} (VCF中不存在)")
            continue
        
        print(f"\n处理染色体: {chrom} (长度: {length})")
        print(f"实际变异范围: {chrom_ranges[chrom]['min']}-{chrom_ranges[chrom]['max']}")
        
        # 计算起始和结束窗口
        first_window = (chrom_ranges[chrom]['min'] - 1) // window_size
        last_window = (chrom_ranges[chrom]['max'] - 1) // window_size
        
        for i in range(first_window, last_window + 1):
            start = i * window_size + 1
            end = min((i + 1) * window_size, length)
            
            # 如果窗口完全在变异范围之外，跳过
            if end < chrom_ranges[chrom]['min'] or start > chrom_ranges[chrom]['max']:
                continue
            
            # 构建输出文件名
            out_file = os.path.join(output_dir, f"{chrom}_{start}-{end}.vcf")
            
            # 构建vcftools命令
            cmd = [
                'vcftools',
                '--vcf', input_vcf,
                '--chr', chrom,
                '--from-bp', str(start),
                '--to-bp', str(end),
                '--recode',
                '--recode-INFO-all',
                '--out', out_file[:-4]  # 去掉.vcf后缀
            ]
            
            # 执行命令
            try:
                subprocess.run(cmd, check=True)
                # 重命名输出文件
                if os.path.exists(out_file[:-4] + ".recode.vcf"):
                    os.rename(out_file[:-4] + ".recode.vcf", out_file)
                    print(f"已创建窗口: {chrom}:{start}-{end}")
                else:
                    print(f"窗口 {chrom}:{start}-{end} 中没有变异，跳过")
            except subprocess.CalledProcessError as e:
                print(f"处理窗口 {chrom}:{start}-{end} 时出错: {e}")

def main():
    parser = argparse.ArgumentParser(description='使用100kb滑窗分割VCF文件')
    parser.add_argument('-i', '--input', required=True, help='输入VCF文件')
    parser.add_argument('-o', '--output', required=True, help='输出目录')
    
    args = parser.parse_args()
    
    # 检查vcftools是否可用
    try:
        subprocess.run(['vcftools', '--version'], 
                      stdout=subprocess.PIPE, 
                      stderr=subprocess.PIPE)
    except FileNotFoundError:
        print("错误: vcftools未安装或不在PATH中")
        sys.exit(1)
    
    # 处理VCF文件
    process_vcf(args.input, args.output)

if __name__ == '__main__':
    main()
