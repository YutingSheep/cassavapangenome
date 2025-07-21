import argparse

def main():
    parser = argparse.ArgumentParser(description='Convert VCF to matrix format')
    parser.add_argument('-i', '--input', required=True, help='Input VCF file')
    parser.add_argument('-o', '--output', required=True, help='Output matrix file')
    args = parser.parse_args()

    with open(args.output, 'w') as out:
        header = ["variant"]
        with open(args.input, 'r') as f:
            for line in f:
                if line.startswith("#"):
                    if line.startswith("#CHROM"):
                        line = line.replace("\n", "").strip().split("\t")
                        for sample in line[9:]:
                            header.append(sample+'_hap1')
                            header.append(sample+'_hap2')
                else:
                    break
        
        out.write(",".join(header) + "\n")
        
        with open(args.input, 'r') as f:
            for line in f:
                if line.startswith("#"):
                    continue
                else:
                    out_gts = []
                    line = line.replace("\n", "").strip().split("\t")
                    variant_id = line[0] + '_' + line[1]
                    out_gts.append(variant_id)
                    gts = line[9:]
                    for gt in gts:
                        haps = gt.split("|")
                        if len(haps) >= 2:  # Make sure we have both haplotypes
                            hap1_gt = haps[0]
                            hap2_gt = haps[1]
                        else:  # Fallback if format is different
                            haps = gt.split("/")  # Try phased format
                            if len(haps) >= 2:
                                hap1_gt = haps[0]
                                hap2_gt = haps[1]
                            else:  # If still no luck, use same value for both
                                hap1_gt = gt
                                hap2_gt = gt
                        out_gts.append(hap1_gt)
                        out_gts.append(hap2_gt)
                    out.write(",".join(out_gts) + "\n")

if __name__ == "__main__":
    main()