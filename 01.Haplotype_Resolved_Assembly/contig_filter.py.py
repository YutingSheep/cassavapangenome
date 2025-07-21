#@Author : Nan Wang github:wangnan9394 
import argparse
import subprocess
from Bio import SeqIO
import os
import random
import time

#需要提前安装minimap2，需要meryl
#需要下载CP.fa and MT.fa,parameters是optimized for Cassava
#input是Hifiasm产生的contig_hap1.fasta or contig_hap2.fasta
#usage: python work.py -f genome.fa -c CP.fa -m MT.fa -t 10 -o ZZZZZ.txt &

parser = argparse.ArgumentParser(description='Process FASTA and align')
parser.add_argument('-f', '--fasta', required=True, help='Query FASTA file')
parser.add_argument('-c', '--req_cp', required=True, help='CP genome file')
parser.add_argument('-t', '--threads', required=True, help='minimap program file')
parser.add_argument('-m', '--req_mt', required=True, help='MT genome file')
parser.add_argument('-o', '--output', required=True, help='Output statistics')
args = parser.parse_args()
# 2024-1-18

print('Analysis starting...')
name_list = []
for record in SeqIO.parse(args.fasta, "fasta"):
    name_list.append(record.id)
    #print(record.id) 
#print(name_list)
print("Name search. Done.")

# Step 1 - Calculate length
lengths = {}
for record in SeqIO.parse(args.fasta, "fasta"):
    seq_id = record.id
    seq_len = len(record)
    lengths[seq_id] = seq_len
#print(lengths)
print("Length caculate. Done.")

gc_dict = {}
for record in SeqIO.parse(args.fasta, "fasta"):
    seq = record.seq
    gc_count = seq.count("G") + seq.count("C")
    gc_content = gc_count / len(seq) * 100
    gc_dict[record.id] = gc_content  
#print(gc_dict)
print("GC caculate. Done.")

cp_pro={}
cmd_cp='minimap2 -t {} -x asm5 {} {} -o {}_{}.paf >CP.log 2>&1'.format(args.threads,args.req_cp,args.fasta,args.fasta,args.req_cp)
os.system(cmd_cp)
for each in name_list:
    s=0
    with open(args.fasta+'_'+args.req_cp+'.paf') as fcp:
        for line in fcp:
            line=line.replace("\n","").split('\t')
            if line[0]==each:
                s+=abs(int(line[3])-int(line[2]))
    cp_pro[each]=str(s/(lengths[each]))
print("Chloroplast genome caculate. Done.")

mt_pro={}
cmd_mt='minimap2 -t {} -x asm5 {} {} -o {}_{}.paf >MT.log 2>&1'.format(args.threads,args.req_mt,args.fasta,args.fasta,args.req_mt)
os.system(cmd_mt)
for each in name_list:
    s=0
    with open(args.fasta+'_'+args.req_mt+'.paf') as fmt:
        for line in fmt:
            line=line.replace("\n","").split('\t')
            if line[0]==each:
                s+=abs(int(line[3])-int(line[2]))
    mt_pro[each]=str(s/(lengths[each]))
print("Mitochondral genome caculate. Done.")

prepared=0
for each in name_list:
    print(each+'.uniq.meryl')
    if os.path.exists(each+'.uniq.meryl'):
        prepared+=1
if not prepared==len(name_list):
    uniq_k={}
    for each in name_list:
        #seed=random.randint(0,100000000)
        cmd_samtools='samtools faidx {} {}>{}'.format(args.fasta,each,each+'.fa')
        os.system(cmd_samtools)
        cmd_meryl_1='meryl threads=1 count compress k=41 {} output {}>{} 2>&1'.format(each+'.fa',each+'.meryl',each+'.log')
        os.system(cmd_meryl_1)
    #time.sleep(60)#暂停1分钟
    #for each in name_list:
        cmd_meryl_2='meryl statistics {} > {}'.format(each+'.meryl',each+'.uniq.meryl')
        os.system(cmd_meryl_2)
    #time.sleep(1200)#暂停20分钟
    for each in name_list:
        count_uk=0
        with open(each+'.uniq.meryl','r') as f_n:
            for line in f_n:
                line=line.replace('\n','')
                if "distinct" and "non-redundant" in line:
                    if not "missing" in line:
                        count_uk=int(line.replace(' ','').split('(')[0].split('t')[-1])
                        uniq_k[each]=1-count_uk/(lengths[each])
                        #uniq_k[each]=count_uk
        #with open(each+'.uniq.meryl','r') as f_n:
        #    for line in f_n:
        #        line=line.replace('\n','').split('\t')[0]
                        #uniq_kerms=int(line)
                        #uniq_k[each]=count_uk
        #print(uniq_kerms,lengths[each])
    #print(uniq_k)
else:
    print("Meryl files prepared!! Continue...")
    uniq_k={}
    for each in name_list:
        with open(each+'.uniq.meryl','r') as f_n:
            for line in f_n:
                line=line.replace('\n','')
                if "distinct" and "non-redundant" in line:
                    if not "missing" in line:
                        count_uk=int(line.replace(' ','').split('(')[0].split('t')[-1])
                        uniq_k[each]=1-count_uk/(lengths[each])
        #print(uniq_kerms,lengths[each])
    #print(uniq_k)
print("Uniq kmers caculate. Done.")

new_name=[]
out=open(args.output,'w')
for each in name_list:
    ll=each+'\t'+str(lengths[each])+'\t'+str(gc_dict[each])+'\t'+str(cp_pro[each])+'\t'+str(mt_pro[each])+'\t'+str(uniq_k[each])+'\n'
    out.write(ll)
    out.flush()
    if int(lengths[each])>100000 and float(gc_dict[each])< 50 and float(cp_pro[each]) < 0.01 and float(mt_pro[each]) < 0.01 and float(uniq_k[each]) < 0.7:
        new_name.append(each)
print("Calculation finished. Done")
out.close()

print("New FASTA generating")
for new_each in new_name:
    new_cmd_samtools='samtools faidx {} {}>>{}'.format(args.fasta,new_each,args.output+'.CLEAN.fa')
    os.system(new_cmd_samtools)
print("All files finished. Done")


    


                                        
