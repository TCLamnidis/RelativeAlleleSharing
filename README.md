# RelativeAlleleSharing
A script to calculate Relative Allele Sharing (RAS) between a Sample population and all other populations in a FreqSum file. The script also calculates RAS per Mb (θ-hat) and preforms weighted Jackknife to give a Jackknife estimate of θ-hat as well as a Jackknife Error estimate for an Error bar.

`Freqsum2RAS.py` requires a FreqSum file containing all the Populations for use in the analysis as an input file. The FreqSum format and various tools for working with this format are documented in https://github.com/stschiff/rarecoal-tools.

# Filtering the FreqSum file
The input FrqSum file should first be filtered through the Mappability, which can be found here: (link to be added here). This filtering needs to be done for each chromosome separately (if using the script provided here), and for that a separate mask for each chromosome will be needed. Such a mask can be easily acquired with:

    for i in {1..22}; do awk -v X=$i '$1==X' /path/to/mask > /path/to/out.chr$i.bed ; done

The resulting `.bed` files can then be used for filtering of each FreqSum. The awk script `filterThroughMask.awk` included in this package, and can be used to filter a FreqSum file using a `.bed` file. 
   
    ...| awk -f filterThroughMask.awk -v maskFile=/path/to/mask.bed | ...

At the moment, each chromosome has to be filtered separately, before merging all of the with `concatFreqSum`.

# Using `FreqSum2RAS.py`
* The input FreqSum file can be input using the option `-I`, or as stdin, by directly piping the FreqSum to `FreqSum2RAS.py` and omitting this option.
* The Bed file used to mask the data is input with the `-B` flag. The length of the chromosomes of the bed file is needed for jackknifing.
* I instead of a bed file, a chromosome length file can be given with the `-C` option. The chromosome file should contain one line per chromosome, containing the chromosome number and the chromosome length in bp. This option is mutually exclusive with `-B`. The Chromosome file format is defined as follows:

```
1   2490000000
2   2200000000
    ...
```

* The output file is defined using the `-O` option.
* The maximum allele frequency for the non-Reference allele in can be changed with the `-M` option (default is 10).
* By default, the minimum allele frequency of the non-reference allele is 2, but it can be changed with option `-m` to analyse other slices of RAS.
* By default, all populations are used to compute allele counts. When using the `--restrictAF` option, you can give a set of comma-separated population names to restrict computing the allele count to those populations.
* The script reads the FreqSum header to catalogue the populations in the FreqSum file. The sample, for which all the RAS will be calculated can be given with `-S` or `--Sample`. Only one population can be given as the Sample population, and the name given should match the name in the FreqSum header.
* With the `-P` flag it is possible to only output private shared variants between the Sample population/individual and another population. 
  * _(It should be noted that at the moment variants found multiple times in the Sample population and multiple times in a Reference population will be counted as private but will also be added to the self-sharing RAS of the Sample population. Hence when finding variants in Allele Frequencies above the number of chromosomes in your Populations, this is the cause. When looking at individuals as Sample population this should not be happening as there shouldn't be any rare variants that are homozygous.)_
* Finally, the `-NT` flag will exclude variants that are transitions from the RAS calculation. 
* The above information can also be found in the command line with the option `-h`.
```
    $FreqSum2RAS.py -h 
usage: FreqSum2RAS.py [-h] [-I <INPUT FILE>] [-M <MAX ALLELE COUNT>]
                      [-m <MIN ALLELE COUNT>] [-O <OUTPUT FILE>]
                      (-C <FILE> | -B <BED FILE>) [-NT] [-P] -S <POPULATION>
                      [--restrictAF POP1,POP2,...]

Extract the frequency of shared rare variants between each test sample/group
and all reference samples/groups from a freqsum file.

optional arguments:
  -h, --help            show this help message and exit
  -I <INPUT FILE>, --Input <INPUT FILE>
                        The input freqsum file. Omit to read stom stdin.
  -M <MAX ALLELE COUNT>, --MAF <MAX ALLELE COUNT>
                        The maximum number of alleles (total) in the reference
                        populations. The default maximum allele value is 10.
  -m <MIN ALLELE COUNT>, --mAF <MIN ALLELE COUNT>
                        The minimum number of alleles (total) in the reference
                        populations. The default minimum allele count is 2.
  -O <OUTPUT FILE>, --Output <OUTPUT FILE>
                        The output file.
  -C <FILE>, --ChromFile <FILE>
                        A file that includes the lengths for each chromosome.
                        The format of this file is Chromosome number Length.
                        Mutually exclusive with -B.
  -B <BED FILE>, --BedFile <BED FILE>
                        The bed file with the calling mask for the FreqSum.THE
                        FREQSUM SHOULD BE FILTERED THROUGH THE MASK BEFORE
                        INPUT. Mutually exclusive with -C.
  -NT                   When present, No Transitions are included in the
                        output. Useful for ancient samples with damaged DNA.
  -P, --Private         Restrict the RAS calculation to privately shared rare
                        variants only.
  -S <POPULATION>, --Sample <POPULATION>
                        Set the Test population/individual. RAS will be
                        calculated between the Test and all populations in the
                        FreqSum.
  --restrictAF POP1,POP2,...
                        Give a list of comma-separated population names that
                        should be considered when computing the allele
                        frequency

```
# RAS Output Format

The first two lines of a RAS output start with a `#` and contain the populations and sizes of the populations in the input FreqSum, and the Sample Population given for that run. 

    #FREQSUM POPS & SIZES: PopA(8) PopB(14) PopC(14) PopD(28) PopE(20) PopF(24) PopG(24)       Individual1(2) Individual2(2) Individual3(2) Individual4(2) Individual5(2)
    #SAMPLE POPULATION:  Individual1

This information is followed by the result tables.
The first line of the result table contains the labels for each result column. 
After this line the script will output a table for each population, containing the above information per allele frequency (2-10 by default) as well as one line for the sum across all allele frequencies.

    RefPop	TestPop			RAS	    RAS/Mb	              Jackknife Estimator     Jackknife Error   Allele Frequency
    PopA    Individual1 	0	0.000000000000000e+00	0.000000000000000e+00	0.000000000000000e+00	2
    PopA    Individual1 	0	0.000000000000000e+00	0.000000000000000e+00	0.000000000000000e+00	3
    PopA    Individual1 	0	0.000000000000000e+00	0.000000000000000e+00	0.000000000000000e+00	4
    PopA    Individual1 	0	0.000000000000000e+00	0.000000000000000e+00	0.000000000000000e+00	5
    PopA    Individual1 	0	0.000000000000000e+00	0.000000000000000e+00	0.000000000000000e+00	6
    PopA    Individual1 	0	0.000000000000000e+00	0.000000000000000e+00	0.000000000000000e+00	7
    PopA    Individual1 	0	0.000000000000000e+00	0.000000000000000e+00	0.000000000000000e+00	8
    PopA    Individual1 	0	0.000000000000000e+00	0.000000000000000e+00	0.000000000000000e+00	9
    PopA    Individual1 	0	0.000000000000000e+00	0.000000000000000e+00	0.000000000000000e+00	10
	PopA    Individual1 	0	0.000000000000000e+00	0.000000000000000e+00	0.000000000000000e+00	Total [2,10]
    
    PopB    Individual1 	0	0.000000000000000e+00	0.000000000000000e+00	0.000000000000000e+00	2
    PopB    Individual1 	0	0.000000000000000e+00	0.000000000000000e+00	0.000000000000000e+00	3
                        ...

The Reference Population is given in the first column, followed by the Sample population in the second. The Sample population will be constant across the whole output file. The third output column shows the RAS statistic as defined in (paper in prep). The fourth output column gives a measure of RAS/Mb, based on the Bed file/Chromosome file provided, and is followed by the Jackknife Estimator of RAS/Mb and the Jackknife Error. The final column shows the frequency of non-reference variants and is between -m (default=2) and -M (default=10). 
_The above output has been formatted slightly to align column labels with the columns._

_Please note that we are still exploring the optimal way to present the output for plotting, and the above format might change in future renditions of the script._
