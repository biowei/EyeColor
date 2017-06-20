'''
Created on 19 Apr 2016

@author: weilu1
'''

from __future__ import division
import collections,sys
import pysam

#Check the command line arguments
if len(sys.argv) < 4:
    print "Usage: <Aligned short reads (bam)> <chromosome position (string)> <Quality (integer)>"
    sys.exit(0)

#Read the reference sequence and initiate the aligner
try:
    samfile = pysam.AlignmentFile(sys.argv[1], "rb")
except IOError as e:
    print "Could not read reference sequence file (see below)!"
    print e
    sys.exit(1)

Quality = int(sys.argv[3])
SNP_poistion = sys.argv[2]

for i, pileupcolumn in enumerate(samfile.pileup(region=SNP_poistion)):
    if pileupcolumn.reference_pos == int(SNP_poistion.split(':')[1])-1:
        total_reads = pileupcolumn.nsegments
        print "There are %d reads overlapping the coordinate at %s." % (total_reads,sys.argv[2])

        PileupBase = collections.namedtuple('PileupBase', ['Pos', 'Base','Qual','Strand'])
        PileupBases = []
        for pileupread in pileupcolumn.pileups:
            k = pileupread.query_position
            p = PileupBase(k, pileupread.alignment.seq[k], pileupread.alignment.query_qualities[k], not pileupread.alignment.is_reverse)
            PileupBases.append(p)

Bases = [p.Base for p in PileupBases]
Bases_f = [p.Base for p in PileupBases if p.Strand == True]

Bases_filter = [p.Base for p in PileupBases if p.Qual >=Quality]
Bases_filter_f = [p.Base for p in PileupBases if p.Strand == False and p.Qual>=Quality]

print "The reads support the following nucleotides (total/with minimum base quality of %d):"% Quality
print "\tA: %d/%d" % (Bases.count('A'),Bases_filter.count('A')),
if Bases.count('A') != 0 and Bases_filter.count('A') != 0:
    print "(of these %d%%/%d%% map to the forward strand)" % (
        Bases_f.count('A')/Bases.count('A')*100,Bases_filter_f.count('A')/Bases_filter.count('A')*100)
else: print
print "\tC: %d/%d" % (Bases.count('C'),Bases_filter.count('C')),
if Bases.count('C') != 0 and Bases_filter.count('C') != 0:
    print "(of these %d%%/%d%% map to the forward strand)" % (
        Bases_f.count('C')/Bases.count('C')*100,Bases_filter_f.count('C')/Bases_filter.count('C')*100)
else: print
print "\tG: %d/%d" % (Bases.count('G'),Bases_filter.count('G')),
if Bases.count('G') != 0 and Bases_filter.count('G') != 0:
    print "(of these %d%%/%d%% map to the forward strand)" % (
        Bases_f.count('G')/Bases.count('G')*100,Bases_filter_f.count('G')/Bases_filter.count('G')*100)
else: print
print "\tT: %d/%d" % (Bases.count('T'),Bases_filter.count('T')),
if Bases.count('T') != 0 and Bases_filter.count('T') != 0:
    print "(of these %d%%/%d%% map to the forward strand)" % (
        Bases_f.count('T')/Bases.count('T')*100,Bases_filter_f.count('T')/Bases_filter.count('T')*100)
else: print
print "\tN: %d/%d" % (Bases.count('N'),Bases_filter.count('N')),
if Bases.count('N') != 0 and Bases_filter.count('N') != 0:
    print "(of these %d%%/%d%% map to the forward strand)" % (
        Bases_f.count('N')/Bases.count('N')*100,Bases_filter_f.count('N')/Bases_filter.count('N')*100)
else: print


"""Computational Genomics Assignment 2

Task 2 Genotyping and phenotyping
	a)	and b) Short summary of sample filtered SNP counting and inferred genotypes.
Results of SNP counting:
WEIs-MBP:Assignment2 weilu1$ python snp_programe.py sample1.bam chr15:28365618 20
There are 36 reads overlapping the coordinate at chr15:28365618.
The reads support the following nucleotides (total/with minimum base quality of 20):
	A: 21/19 (of these 38%/57% map to the forward strand)
	C: 0/0
	G: 15/15 (of these 40%/60% map to the forward strand)
	T: 0/0
	N: 0/0
WEIs-MBP:Assignment2 weilu1$ python snp_programe.py sample2.bam chr15:28365618 20
There are 25 reads overlapping the coordinate at chr15:28365618.
The reads support the following nucleotides (total/with minimum base quality of 20):
	A: 25/25 (of these 48%/52% map to the forward strand)
	C: 0/0
	G: 0/0
	T: 0/0
	N: 0/0
WEIs-MBP:Assignment2 weilu1$ python snp_programe.py sample3.bam chr15:28365618 20
There are 13 reads overlapping the coordinate at chr15:28365618.
The reads support the following nucleotides (total/with minimum base quality of 20):
	A: 0/0
	C: 0/0
	G: 13/13 (of these 23%/76% map to the forward strand)
	T: 0/0
	N: 0/0

WEIs-MBP:Assignment2 weilu1$ python snp_programe.py sample1.bam chr15:28356859 20
There are 38 reads overlapping the coordinate at chr15:28356859.
The reads support the following nucleotides (total/with minimum base quality of 20):
	A: 1/0
	C: 16/15 (of these 37%/60% map to the forward strand)
	G: 0/0
	T: 21/18 (of these 47%/50% map to the forward strand)
	N: 0/0
WEIs-MBP:Assignment2 weilu1$ python snp_programe.py sample2.bam chr15:28356859 20
There are 20 reads overlapping the coordinate at chr15:28356859.
The reads support the following nucleotides (total/with minimum base quality of 20):
	A: 0/0
	C: 19/18 (of these 63%/38% map to the forward strand)
	G: 1/0
	T: 0/0
	N: 0/0
WEIs-MBP:Assignment2 weilu1$ python snp_programe.py sample3.bam chr15:28356859 20
There are 9 reads overlapping the coordinate at chr15:28356859.
The reads support the following nucleotides (total/with minimum base quality of 20):
	A: 0/0
	C: 0/0
	G: 0/0
	T: 9/9 (of these 22%/77% map to the forward strand)
	N: 0/0
"""

"""Summary as below:
Sample 1
Genotype
rs12913832
A/G heterozygous
rs1129038
C/T heterozygous

Sample 2
Genotype
rs12913832
A/A homozygous
rs1129038
C/C homozygous

Sample 3
Genotype
rs12913832
G/G homozygous
rs1129038
T/T homozygous
"""

'''

Genotypes (rs1129038, rs12913832)
Eye Color
Sample 1
A/G, C/T
Brown

Sample 2
A/A, C/C
Blue

Sample 3
G/G, T/T
Brown

c) In this question we only consider two SNP alleles. 
Similarly we can perform Fisher's exact test on our data, the p-value is 2.2e-16. 
This test p-value suggests that the four genotypes and eye color are strongly associated.
By comparing sample's genotypes with the above table we can confidently infer the eye color of each sample.

d) Probability of eye color phenotyping
Inferring eye color from genotypes is essentially a conditional probability calculate, 
namely probability of sample's eye color given that certain genotype happens. 
By definition of conditional probability we can calculate by following 
P(eyecolor|genotype)=P(genotype|eyecolor)/P(genotype).

In this example we estimate P(genotype) based on limited sample size. 
So we can only infer the eye color with certain probability. Further, 
multiple genes may control the eye color. In this inference we only consider two loci of these genes. 
So we can only infer the eye color with certain probabilities.  

'''


