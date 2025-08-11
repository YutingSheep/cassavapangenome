# -*- coding: utf-8 -*-
"""
Created on March 12 2021

@author: Nan
"""

import vcf
import argparse

parser = argparse.ArgumentParser(description = 'For mandarin', add_help = False, usage = '\npython3 -i [input.vcf] -g [outgroup_list] -o [output.vcf]')
required = parser.add_argument_group()
optional = parser.add_argument_group()
required.add_argument('-i', '--input', metavar = '[input_vcf]', help = 'input_vcf', required = True)
required.add_argument('-g', '--group', metavar = '[outgroup_list]', help = 'outgroup_list', required = True)
required.add_argument('-o', '--output', metavar = '[output]', help = 'output', required = True)
optional.add_argument('-h', '--help', action = 'help', help = 'help')
args = parser.parse_args()
vcf_reader = vcf.Reader(open(args.input, 'r'))

def rever_SNP(gt):
    if gt=='0/0':
        new_gt='1/1'
    if gt=='0/1':
        new_gt='0/1'
    if gt=='1/1':
        new_gt='0/0'
    if gt =='./.':
        new_gt='./.'
    if gt =='./1':
        new_gt='./0'
    if gt =='./0':
        new_gt='./1'
    if gt =='0/.':
        new_gt='1/.'
    if gt =='1/.':
        new_gt='0/.'
    return new_gt

outGrooup_file=[]
with open(args.group,'r') as outgroup_name:
    for each in outgroup_name:
        each=each.replace("\n","").split('\t')[1]
        outGrooup_file.append(each)
outGrooup_file=set(outGrooup_file)
dict_one={}
for one in outGrooup_file:
    list_one=[]
    with open(args.group,'r') as outgroup_name:
        for each in outgroup_name:
            if each.replace("\n","").split('\t')[1]==one:
                list_one.append(each.replace("\n","").split('\t')[0])
    dict_one[one]=list_one
#print(dict_one)             

out_vcf=open(args.output,'w')
with open(args.input,'r') as f:#set3.234.deepvariants.vcf
    for line in f:
        if line.startswith('#'):
            out_vcf.write(line)
        else:
            break

for record in vcf_reader:
    #print(record)
    chr=record.CHROM
    pos=record.POS
    ref=record.REF
    alt=record.ALT[0]
    #print(chr,pos,ref,alt)
    id=str(chr)+'_'+str(pos)+'_'+ref
    qual=record.QUAL
    filter='.'
    info='.'
    format=record.FORMAT
    
    list_freq=[]
    for qw,samp in dict_one.items():
        #print(qw,samp)
        out_group_freq=0
        for i in samp:
            gt=record.genotype(i)['GT'].split('/')
            if '1' in gt:
                out_group_freq+=1
        out_group_freq_prop=out_group_freq/(2*len(samp))
        list_freq.append(out_group_freq_prop)
    tag=0
    #print(list_freq)
    for i in list_freq:
        if i >0:
            tag+=1
    #print(tag)
    if tag==len(outGrooup_file):
        list=[]
        for i in record.samples:
                index=i['GT']
                new_index=rever_SNP(index)
                list.append(new_index)
        l2=''
        for each in list:
            l2+=each+'\t'
    else:
        list=[]
        for i in record.samples:
            index=i['GT']
            list.append(index)
        l2=''
        for each in list:
            l2+=each+'\t'
    ll=str(chr)+'\t'+str(pos)+'\t'+id+'\t'+ref+'\t'+str(alt)+'\t'+str(qual)+'\t'+filter+'\t'+info+'\t'+format+'\t'+l2+'\n'
    out_vcf.write(ll)
