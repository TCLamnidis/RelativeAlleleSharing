#!/usr/bin/env python3

import sys, argparse, re
from math import sqrt
from time import strftime

parser = argparse.ArgumentParser(description="Extract the frequency of shared rare variants between each test sample/group and all reference samples/groups from a freqsum file.")
parser.add_argument("-I", "--Input", metavar="<INPUT FILE>", type=argparse.FileType('r'), help="The input freqsum file.", required=False)
parser.add_argument("-M", "--MAF", metavar="<MAX ALLELE COUNT>", type=int, default=10, help="The maximum number of alleles (total) in the reference populations. The minimum allele count is always 2.", required=False)
# parser.add_argument("-mN", "--NormMin", metavar="<NUMBER>", type=int, help="The minimum number of alleles to be taken into account for the normalization factor Zα. Should be equivalent to allele frequency of about 1%% in the reference populations.", required=True)
# parser.add_argument("-MN", "--NormMax", metavar="<NUMBER>", type=int, help="The maximum number of alleles to be taken into account for the normalization factor Zα. Should be equivalent to allele frequency of about 10%% in the reference populations.", required=True)
parser.add_argument("-O", "--Output", metavar="<OUTPUT FILE>", type=argparse.FileType('w'), help="The output file.", required=True)
parser.add_argument("-B", "--BedFile", metavar="<BED FILE>", type=argparse.FileType('r'), help="The bed file with the calling mask for the FreqSum.THE FREQSUM SHOULD BE FILTERED THROUGH THE MASK BEFORE INPUT.", required=True)
parser.add_argument("-NT", action='store_true', help="When present, No Transitions are included in the output. Useful for ancient samples with damaged DNA.")
parser.add_argument("-P","--Private", action='store_true', required=False, help="Restrict the RAS calculation to privately shared rare variants only.")
parser.add_argument("-S", "--Sample", type=str, metavar="<POPULATION>", required=True, help="Set the Test population/individual. RAS will be calculated between the Test and all populations in the FreqSum.")
args = parser.parse_args()

print ("Program began running at:", strftime("%D %H:%M:%S"), file=sys.stderr)
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
NumGroups=22

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
        JackknifeMatrix=[[[0 for i in range(NumGroups)] for j in range(M+1)] for k in range(len(Names))]
        ras=[[[0 for i in range(NumGroups)] for j in range(M+1)] for k in range(len(Names))]
        mj=[[[0 for i in range(NumGroups)] for j in range(M+1)] for k in range(len(Names))]
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
                    JackknifeMatrix[r-4][Sum][Chr-1]+=(c1*(c2-1)) / (Sizes[Names[r]] * (Sizes[Names[Test]]-1))
                    mj[r-4][Sum][Chr-1]+=1
                else:
                    JackknifeMatrix [r-4][Sum][Chr-1]+=(c1*c2) / (Sizes[Names[r]] * Sizes[Names[Test]])
                    mj[r-4][Sum][Chr-1]+=1

#Reading chr lengths from bed file
lengths=[0 for x in range(NumGroups)]
for line in args.BedFile:
    fields = line.strip().split()
    Chr=int(fields[0])-1
    start = int(fields[1])
    end = int(fields [2])
    lengths [Chr]+= (end - start)/1000000


#Calculation of ras, RAS per Mb
for j in range(NumGroups):
    for k in range(M+1):
        for l in range(len(Names)):
            ras[l][k][j]=(JackknifeMatrix[l][k][j]/lengths[j])

#Jackknife stimation
Thetahat = [[0 for j in range(M+1)] for k in range(len(Names))]
Thetaminus=[[[0 for c in range(NumGroups)] for j in range(M+1)] for k in range(len(Names))]
for i in range(2,M+1):
    for j in range(len(Names)):
        Thetahat[j][i]=(sum(JackknifeMatrix[j][i])/sum(lengths))
        for c in range(NumGroups):
            Thetaminus[j][i][c]=(sum(ras[j][i])-ras[j][i][c])

ThetaJ=[[0 for j in range(M+1)] for k in range(len(Names))]
for i in range(2,M+1):
    for j in range(len(Names)):
        Sum1=0
        Sum2=0
        for c in range(NumGroups):
            Sum1+=Thetahat[j][i]-Thetaminus[j][i][c]
            Sum2+=((mj[j][i][c] * Thetaminus[j][i][c])/sum(mj[j][i]))
        ThetaJ[j][i]=Sum1+Sum2

Sigma2=[[0 for j in range(M+1)] for k in range(len(Names))]
for i in range(2,M+1):
    for j in range(len(Names)):
        for c in range(NumGroups):
            hj=sum(mj[j][i])/mj[j][i][c]
            pseudovalue=(hj*Thetahat[j][i])-((hj-1) * Thetaminus[j][i][c])
            Sigma2[j][i]+=(((pseudovalue-ThetaJ[j][i])**2)/(hj-1))/NumGroups

#Print output tables
print (*PopNames, file=args.Output, sep="\t", end="\n")
print ("#SAMPLE POPULATION: ", args.Sample, file=args.Output, end="\n\n")
for m in range(2,M+1):
    print(m,"RAS","Jackknife Error", sep="\t", file=args.Output)
    for i in Refs:
        print (Names[i], sum(JackknifeMatrix[i-4][m]), sqrt(Sigma2[i][m]), sep="\t", file=args.Output)
    print ("", file=args.Output)

print ("Program finished running at:", strftime("%D %H:%M:%S"), file=sys.stderr)







