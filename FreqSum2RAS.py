#!/usr/bin/env python3

import sys, argparse, re
from time import strftime


parser = argparse.ArgumentParser(description="Extract the frequency of shared rare variants between each test sample/group and all reference samples/groups from a freqsum file.")
parser.add_argument("-I", "--Input", metavar="<INPUT FILE>", type=argparse.FileType('r'), help="The input freqsum file.", required=False)
parser.add_argument("-M", "--MAF", metavar="<MAX ALLELE COUNT>", type=int, default=10, help="The maximum number of alleles (total) in the reference populations. The minimum allele count is always 2.", required=False)
# parser.add_argument("-mN", "--NormMin", metavar="<NUMBER>", type=int, help="The minimum number of alleles to be taken into account for the normalization factor Zα. Should be equivalent to allele frequency of about 1%% in the reference populations.", required=True)
# parser.add_argument("-MN", "--NormMax", metavar="<NUMBER>", type=int, help="The maximum number of alleles to be taken into account for the normalization factor Zα. Should be equivalent to allele frequency of about 10%% in the reference populations.", required=True)
parser.add_argument("-O", "--Output", metavar="<OUTPUT FILE>", type=argparse.FileType('w'), help="The output file.", required=True)
parser.add_argument("-NT", action='store_true', help="When present, No Transitions are included in the output. Useful for ancient samples with damaged DNA.")
# group = parser.add_mutually_exclusive_group(required=False)
# group.add_argument("-L", "--SampleList", type=argparse.FileType('r'), metavar="<INDIVIDUAL LIST FILE>", required=False, help="A list of column Names from the input FreqSum file to be excluded from the reference group. Ancient individuals should be excluded from the reference group. Can be supplemented with -S.")
# group.add_argument("-S", "--Sample", action="append", metavar="<INDIVIDUAL>", required=False, help="Same as -L, except for flagging columns individually from the command line. Can be called multiple times. Can be called alongside -L to flag additional columns to those in the list.")
args = parser.parse_args()

#If no input file given, read from stdin
if args.Input == None:
    I = sys.stdin
else:
    I = args.Input
M=args.MAF
Transitions = {"A":"G", "G":"A","C":"T","T":"C"}
Samples=[]
Refs=[]
Tests=[]
Sizes={}
Names={}
# if args.SampleList != None:
#     for line in args.SampleList:
#         Samples.append(line.strip())
#
# if args.Sample!=None:
#     for i in args.Sample:
#         Samples.append(i)

for line in args.Input:
    fields=line.strip().split()
    if fields[0][0]=="#":
        PopNames=fields
        for i in range(4,len(fields)):
            if re.split('[(|)]',fields[i])[0] not in Samples:
                Refs.append(i)
                Names[i]=re.split('[(|)]',fields[i])[0]
                Sizes[re.split('[(|)]',fields[i])[0]]=int(re.split('[(|)]',fields[i])[1])
            else:
                Tests.append(i)
                Names[i]=re.split('[(|)]',fields[i])[0]
                Sizes[re.split('[(|)]',fields[i])[0]]=int(re.split('[(|)]',fields[i])[1])
        Matrix=[[[0 for i in range(M+1)] for j in range(len(Names))] for k in range(len(Names))]
    else:
        if fields[2]=="N":
            continue
        if args.NT == True:
            if fields[3]==Transitions[fields[2]]:
                continue
        for i in range(4,len(fields)):
            fields[i]=int(fields[i])
            if fields[i]<0:
                fields[i]=0
        Sum=0
        Sum=sum(fields[4:])
        if Sum>M:
            continue
        elif Sum<2:
            continue
        else:
            for r in Refs:
                c1=fields[r]
                if c1==0:
                    continue
                for x in Refs:
                    c2 = fields[x]
                    if c2==0:
                        continue
                    elif r==x:
                        Matrix[r-4][x-4][Sum]+=(c1*(c2-1)) / (Sizes[Names[r]] * (Sizes[Names[x]]-1))
                    else:
                        Matrix [r-4][x-4][Sum]+=(c1*c2) / (Sizes[Names[r]] * Sizes[Names[x]])

for m in range(2,M+1):
    print(m,*(Names[i] for i in Refs), sep="\t", file=args.Output)
    for i in Refs:
        print (Names[i], *(Matrix[i-4][x-4][m] for x in Refs), sep="\t", file=args.Output)
    print ("", file=args.Output)








