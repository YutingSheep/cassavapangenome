import argparse

parser = argparse.ArgumentParser(description='For clusters', add_help=False, 
                               usage='python work.py -i [input.clusters] -o [output.clusters]')
required = parser.add_argument_group('required arguments')
optional = parser.add_argument_group('optional arguments')
required.add_argument('-i', '--input', metavar='INPUT', help='input clusters file', required=True)
required.add_argument('-o', '--output', metavar='OUTPUT', help='output file', required=True)
optional.add_argument('-h', '--help', action='help', help='show this help message and exit')

args = parser.parse_args()

with open(args.output, "w") as out_line:
    list1 = []
    with open(args.input, "r") as f:
        for line in f:
            line = line.strip().replace("\n", " ").split(",")
            list1.append(line[0])
    out_line.write(",".join(list1) + "\n")
    
    list2 = [args.input]
    with open(args.input, "r") as f:
        for line in f:
            line = line.strip().replace("\n", " ").split(",")
            list2.append(line[2])
    list2.remove(list2[1])
    out_line.write(",".join(list2) + "\n")