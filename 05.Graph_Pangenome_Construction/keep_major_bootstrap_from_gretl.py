import numpy as np

for chr in range(1,19):
    if chr<10:
        chr_file="0"+str(chr)+'.bootstrap'
    else:
        chr_file=str(chr)+'.bootstrap'
    outfile=open(chr_file+'.out',"w")
    with open(chr_file,"r") as f:
        for line in f:
            line=line.replace("\n","").split('\t')

            if line[2]=="S":
                int_list = [int(x) for x in line[3:]]
                ll=line[0]+"\t"+str(np.sum(np.array(int_list)))+'\n'
                print(ll)
                outfile.write(ll)
    outfile.close()
