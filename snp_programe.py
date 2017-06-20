'''
Created on 19 Apr 2016
@author: weilu1
This script is used to analysis genotypes of eye color associated SNPs and infer eye color 
using a simple probabilty model.
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