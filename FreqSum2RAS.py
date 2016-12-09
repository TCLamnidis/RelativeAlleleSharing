#!/usr/bin/env python3

import sys, argparse, re


parser = argparse.ArgumentParser(description="Extract the frequency of shared rare variants between each test sample/group and all reference samples/groups from a freqsum file.")
parser.add_argument("-I", "--Input", metavar="<INPUT FILE>", type=argparse.FileType('r'), help="The input freqsum file.", required=False)
parser.add_argument("-M", "--MAF", metavar="<MAX ALLELE COUNT>", type=int, default=10, help="The maximum number of alleles (total) in the reference populations. The minimum allele count is always 2.", required=False)
# parser.add_argument("-mN", "--NormMin", metavar="<NUMBER>", type=int, help="The minimum number of alleles to be taken into account for the normalization factor Zα. Should be equivalent to allele frequency of about 1%% in the reference populations.", required=True)
# parser.add_argument("-MN", "--NormMax", metavar="<NUMBER>", type=int, help="The maximum number of alleles to be taken into account for the normalization factor Zα. Should be equivalent to allele frequency of about 10%% in the reference populations.", required=True)
parser.add_argument("-O", "--Output", metavar="<OUTPUT FILE>", type=argparse.FileType('w'), help="The output file.", required=True)
parser.add_argument("-NT", action='store_true', help="When present, No Transitions are included in the output. Useful for ancient samples with damaged DNA.")
parser.add_argument("-P","--Private", action='store_true', required=False, help="Restrict the RAS calculation to privately shared rare variants only.")
parser.add_argument("-S", "--Sample", type=str, metavar="<POPULATION>", required=True, help="Set the Test population/individual. RAS will be calculated between the Test and all populations in the FreqSum.")
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

#Read -S argumant into sample list
if args.Sample!=None:
    Samples.append(args.Sample)

for line in args.Input:
    fields=line.strip().split()
    #Use FreqSum header to extract Test and Reference Pops
    if fields[0][0]=="#":
        PopNames=fields
        for i in range(4,len(fields)):
            if re.split('[(|)]',fields[i])[0] in Samples:
                Test=int(i)
            Refs.append(i)
            Names[i]=re.split('[(|)]',fields[i])[0]
            Sizes[re.split('[(|)]',fields[i])[0]]=int(re.split('[(|)]',fields[i])[1])
        #Define RAS Matrix
        JackknifeMatrix=[[[0 for i in range(22)] for j in range(M+1)] for j in range(len(Names))]
    else:
        #Ignore sites where hg19 has Ns
        if fields[2]=="N":
            continue
        #Exclude Transitions
        if args.NT == True:
            if fields[3]==Transitions[fields[2]]:
                continue
        #Convert Freqsum input from string to integers
        fields[0]=int(fields[0])
        for i in range(4,len(fields)):
            fields[i]=int(fields[i])
            if fields[i]<0:
                fields[i]=0
        #Only analyse lines where the Test population has variants
        if fields[Test]==0:
            continue
        #Calculate variant counts in all populations and exclude if above MaxAf or 1
        Sum=0
        Sum=sum(fields[4:])
        if Sum>M:
            continue
        elif Sum<2:
            continue
        else:
            Chr=fields[0]
            #Calcualte RAS for each population compared to the test
            c2 = fields[Test]
            for r in Refs:
                c1=fields[r]
                #Exclude non private variants if private flag is given
                if args.Private:
                    if Sum!=c1+c2:
                        continue
                if c1==0:
                    continue
                elif r==Test:
                    JackknifeMatrix[r-4][Sum][Chr]+=(c1*(c2-1)) / (Sizes[Names[r]] * (Sizes[Names[Test]]-1))
                else:
                    JackknifeMatrix [r-4][Sum][Chr]+=(c1*c2) / (Sizes[Names[r]] * Sizes[Names[Test]])

#Print output tables
for m in range(2,M+1):
    print(m,Names[Test], sep="\t", file=args.Output)
    for i in Refs:
        print (Names[i], sum(JackknifeMatrix[i-4][m]), sep="\t", file=args.Output)
    print ("", file=args.Output)








