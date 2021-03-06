#!/usr/bin/env python3

import sys, argparse, re
from math import sqrt
from time import strftime

parser = argparse.ArgumentParser(description="Extract the frequency of shared rare variants between each test sample/group and all reference samples/groups from a freqsum file.")
parser.add_argument("-I", "--Input", metavar="<INPUT FILE>", type=argparse.FileType('r'), help="The input freqsum file. Omit to read stom stdin.", required=False)
parser.add_argument("-M", "--MAF", metavar="<MAX ALLELE COUNT>", type=int, default=10, help="The maximum number of alleles (total) in the reference populations. The default maximum allele value is 10.", required=False)
parser.add_argument("-m", "--mAF", metavar="<MIN ALLELE COUNT>", type=int, default=2, help="The minimum number of alleles (total) in the reference populations. The default minimum allele count is 2.", required=False)
# parser.add_argument("-mN", "--NormMin", metavar="<NUMBER>", type=int, help="The minimum number of alleles to be taken into account for the normalization factor Zα. Should be equivalent to allele frequency of about 1%% in the reference populations.", required=True)
# parser.add_argument("-MN", "--NormMax", metavar="<NUMBER>", type=int, help="The maximum number of alleles to be taken into account for the normalization factor Zα. Should be equivalent to allele frequency of about 10%% in the reference populations.", required=True)
group = parser.add_mutually_exclusive_group(required=True)
parser.add_argument("-O", "--Output", metavar="<OUTPUT FILE>", type=argparse.FileType('w'), help="The output file.")
group.add_argument("-C", "--ChromFile", metavar="<FILE>", type=argparse.FileType('r'), help="A file that includes the lengths for each chromosome. The format of this file is Chromosome number    Length. Mutually exclusive with -B.", required=False)
group.add_argument("-B", "--BedFile", metavar="<BED FILE>", type=argparse.FileType('r'), help="The bed file with the calling mask for the FreqSum.THE FREQSUM SHOULD BE FILTERED THROUGH THE MASK BEFORE INPUT. Mutually exclusive with -C.", required=False)
parser.add_argument("-NT", action='store_true', help="When present, No Transitions are included in the output. Useful for ancient samples with damaged DNA.")
parser.add_argument("-P","--Private", action='store_true', required=False, help="Restrict the RAS calculation to privately shared rare variants only.")
parser.add_argument("-S", "--Sample", type=str, metavar="<POPULATION>", required=True, help="Set the Test population/individual. RAS will be calculated between the Test and all populations in the FreqSum.")
parser.add_argument("--restrictAF", type=str, metavar="POP1,POP2,...", required=False, help="Give a list of comma-separated population names that should be considered when computing the allele frequency")
args = parser.parse_args()

if args.ChromFile is True and args.Bedfile is True:
	parser.error("--BedFile and --ChromFile are mutually exclusive.")

if args.mAF<1:
    parser.error("--mAF cannot be lower than 1.")

print ("Program began running at:", strftime("%D %H:%M:%S"), file=sys.stderr)
#If no input file given, read from stdin
if args.Input == None:
    Input = sys.stdin
else:
    Input = args.Input

if args.Output == None:
    args.Output = sys.stdout

mAF=args.mAF
M=args.MAF
Transitions = {"A":"G", "G":"A","C":"T","T":"C"}
Samples=[]
Refs=[]
Tests=[]
Sizes={}
Names={}

restrictPops = args.restrictAF.split(",") if args.restrictAF != None else []
    
def convert_lengths_dict_to_lengths_list(lengths_dict):
    print(lengths_dict)
    nrChroms = len(lengths_dict)
    lengths = []
    for chrom in range(nrChroms):
        if chrom not in lengths_dict:
            print("Error: could not find length information for chromosome",
                chrom, file=sys.stderr)
            sys.exit()
        else:
            lengths.append(lengths_dict[chrom])
    return lengths

def read_chrom_file(filename):
    lengths_dict = {}
    for line in filename:
        fields = line.strip().split()
        Chr=int(fields[0])-1
        lengths_dict[Chr]=int(fields[1])/1000000
    return convert_lengths_dict_to_lengths_list(lengths_dict)

def read_bed_file(filename):
    lengths_dict = {}
    for line in filename:
        fields = line.strip().split()
        Chr=int(fields[0])-1
        start = int(fields[1])
        end = int(fields [2])
        if Chr not in lengths_dict:
            lengths_dict[Chr] = 0
        lengths_dict[Chr] += (end - start)/1000000
    return convert_lengths_dict_to_lengths_list(lengths_dict)

lengths = []
if args.ChromFile:
    lengths = read_chrom_file(args.ChromFile)
else:
    lengths = read_bed_file(args.BedFile)

NumBins=len(lengths)

print("found", NumBins, "chromosomes in bed/chrom file", file=sys.stderr)

def read_Freqsum_Header():
    global Test
    global PopNames
    Pops=fields[4:]
    PopNames=Pops
    for i in range(len(Pops)):
        if re.split('[(|)]',Pops[i])[0] in Samples:
            Test=int(i)
        Refs.append(i)
        Names[i]=re.split('[(|)]',Pops[i])[0]
        Sizes[re.split('[(|)]',Pops[i])[0]]=int(re.split('[(|)]',Pops[i])[1])

#Read -S argument into sample list
if args.Sample!=None:
    Samples.append(args.Sample)

lastChrom = -1
for line in Input:
    fields=line.strip().split()
    Position=fields[1]
    #Use FreqSum header to extract Test and Reference Pops
    if fields[0][0]=="#":
        read_Freqsum_Header()

        #Define RAS Matrix
        RAS=[[[0 for i in range(NumBins)] for j in range(M+1)] for k in range(len(Names))]
        mj=[[[0 for i in range(NumBins)] for j in range(M+1)] for k in range(len(Names))]

    else:
        #Ignore sites where hg19 has Ns
        if fields[2]=="N":
            continue
        #Exclude Transitions
        if args.NT == True:
            if fields[3]==Transitions[fields[2]]:
                continue
        #Convert Freqsum input from string to integers
        if fields[0][0:3]=="chr":
            Chr=int(fields[0][3:])-1
        else:
            Chr=int(fields[0])-1
        if Chr != lastChrom:
            print("processing chromosome", Chr, file=sys.stderr)
        lastChrom = Chr

        Data = [0 if f=="-1" else int(f) for f in fields[4:]] #Assume ref allele when missing data!
        #Only analyse lines where the Test population has variants
        if Data[Test]==0:
            continue
        #Calculate variant counts in all populations and exclude if above MaxAf or 1
        Sum=0
        if restrictPops == []:
            Sum=sum(Data)
        else:
            for i in range(len(Names)):
                if Names[i] in restrictPops:
                    Sum += Data[i] 
        if Sum>M:
            continue
        elif Sum<2:
            continue
        else:
            #Calcualte RAS for each population compared to the test
            c2 = Data[Test]
            for r in Refs:
                c1=Data[r]
                #Exclude non private variants if private flag is given
                if args.Private:
                    if Sum!=c1+c2:
                        continue
                if c1==0:
                    continue
                elif r==Test:
                    RAS[r][Sum][Chr]+=(c1*(c2-1)) / (Sizes[Names[r]] * (Sizes[Names[Test]]-1))
                    mj[r][Sum][Chr]+=1
                    RAS[r][mAF-1][Chr]+=(c1*(c2-1)) / (Sizes[Names[r]] * (Sizes[Names[Test]]-1)) # "mAF-1" stores the sum of RAS across all chromosomes.
                else:
                    RAS [r][Sum][Chr]+=(c1*c2) / (Sizes[Names[r]] * Sizes[Names[Test]])
                    mj[r][Sum][Chr]+=1
                    RAS [r][mAF-1][Chr]+=(c1*c2) / (Sizes[Names[r]] * Sizes[Names[Test]]) # "mAF-1" stores the sum of RAS across all chromosomes.

#Jackknife stimation
Thetahat = [[0 for j in range(M+1)] for k in range(len(Names))]
Thetaminus=[[[0 for c in range(NumBins)] for j in range(M+1)] for k in range(len(Names))]
for i in range(mAF-1,M+1): # M+1 to pick up all chromosomes (0-based to 1-based). mAF-1 to get the Thetas for the Sum of AFs too.
    for j in range(len(Names)):
        Thetahat[j][i]=(sum(RAS[j][i])/sum(lengths))
        for c in range(NumBins):
            Thetaminus[j][i][c]=(sum(RAS[j][i]) - RAS[j][i][c]) / (sum(lengths) - lengths[c])

ThetaJ=[[0 for j in range(M+1)] for k in range(len(Names))]
for i in range(mAF-1,M+1): # M+1 to pick up all chromosomes (0-based to 1-based). mAF-1 to get the Thetas for the Sum of AFs too.
    for j in range(len(Names)):
        Sum1=0
        Sum2=0
        for c in range(NumBins):
            Sum1+=Thetahat[j][i]-Thetaminus[j][i][c]
            Sum2+=((lengths[c] * Thetaminus[j][i][c])/sum(lengths))
        ThetaJ[j][i]=Sum1+Sum2

Sigma2=[[0 for j in range(M+1)] for k in range(len(Names))]
for i in range(mAF-1,M+1): # M+1 to pick up all chromosomes (0-based to 1-based). mAF-1 to get the Thetas for the Sum of AFs too.
    for j in range(len(Names)):
        for c in range(NumBins):
            hj=sum(lengths)/lengths[c]
            pseudovalue=(hj*Thetahat[j][i])-((hj-1) * Thetaminus[j][i][c])
            Sigma2[j][i]+=(((pseudovalue-ThetaJ[j][i])**2)/(hj-1))/NumBins

#Print output tables
print ("#FREQSUM POPS & SIZES:",*PopNames, file=args.Output, sep=" ", end="\n")
print ("#SAMPLE POPULATION: ", Names[Test], file=args.Output, end="\n")
if restrictPops == []:
    print ("#POPULATIONS CONSIDERED FOR ALLELE FREQUENCY CALCULATIONS:", "ALL", file=args.Output, sep="\t", end="\n\n")
else:
    print ("#POPULATIONS CONSIDERED FOR ALLELE FREQUENCY CALCULATIONS:", *restrictPops, file=args.Output, sep="\t", end="\n\n")
print("RefPop","TestPop","RAS","RAS/Mb","Jackknife Estimator", "Jackknife Error ", "Allele Frequency", sep="\t", file=args.Output)
for i in Refs:
    for m in range(mAF,M+1):
        print (Names[i], Names[Test], "{:.5}".format(float(sum(RAS[i][m]))), "{:.15e}".format(Thetahat[i][m]), "{:.15e}".format(ThetaJ[i][m]), "{:.15e}".format(sqrt(Sigma2[i][m])),m, sep="\t", file=args.Output)
    m=mAF-1
    print (Names[i], Names[Test], "{:.5}".format(float(sum(RAS[i][m]))), "{:.15e}".format(Thetahat[i][m]), "{:.15e}".format(ThetaJ[i][m]), "{:.15e}".format(sqrt(Sigma2[i][m])),"Total [{},{}]".format(mAF,M), sep="\t", file=args.Output)
    print ("", file=args.Output)

print ("Program finished running at:", strftime("%D %H:%M:%S"), file=sys.stderr)
