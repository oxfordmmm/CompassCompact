#!/usr/bin/python

# ============================================================================ #
# g4_merge.py - Reimplementation of gorm_mergevcfs.py
# Main changes: Data sctructures, OOP
#
# Execution of same sample in Bespin [AMD Opteron(TM) Processor 6274 - 32 Core - 256 Gb Ram]
#   gorm_mergevcfs.py: 8m11s - 5.3Gb Ram
#         g4_merge.py: 5m37s - 171Mb Ram   Improvements: 30% and 96% respectively
#
# ============================================================================ #
# Camilla Ip
# camilla.ip@stats.ox.ac.uk
# August 2012
#
# Carlos del Ojo Elias
# carlos.delojoelias@ndm.ox.ac.uk
# Juny 2013
# ============================================================================ #

import argparse

"""
	Merge VCF files containing different information into a single VCF
	
	The Gorm Mapcall pipeline creates a set of VCF files with the following
	information with the prefix comid_refix.mapcall.SUFFIX:
	
	  4.varannotDP.vcf : the genotype (FORMAT:GT) for sites with some depth, but
						 no high-quality depth (INFO:DP4)
	  4.varannot.vcf   : the high-quality depth and all of the values for variant sites
	  5.indels.vcf	 : the indel records
	  2.mpileupDP.txt  : the number of spanning reads (INFO:DP)
	
	Before now, the Mapcall pipeline has been run to obtain the number of spanning
	reads (INFO:DP) as well as the statistics based on high-quality depth (INFO:DP4)
	and the reference GC (INFO:GC).
	
	So this program starts with the 4.varannot.vcf file containing the high-quality
	depth values then:
	  - relabel MQ with MQ4
	  - replace DP with DP from 4.varannotDP.vcf
	  - replace BaseCounts with BaseCounts from 4.varannotDP.vcf
	  - add MQ from 4.varannotDP.vcf
	
	Note that the 4.varannotDP.vcf file has more records than the 4.varannot.vcf
	file because these files only contain the sites with non-zero depth, and the
	there are 1.5% of the reference sites with zero high-depth coverage but non-zero
	number of spanning reads. For these cases, the VCF record should contain sensible
	INFO and FORMAT values inferred from the 4.varannotDP.vcf file.
	
	The output.bam file is called
	
	  7.merged.vcf	 : a single VCF file containing all statistics on the mapped
						 reads required to assess the callability of each site
	
"""

from extendibleArray import ExtendibleArray
from fastafile import FastaFile
from gormVCF import GormVcf
import array
import itertools
import sys
from genomemap import GenomeMap


class GormVcfMerger:
    '''		The categories are defined as:
            nilspanning                    allsites - inDP
            yesspanning_nilmapped            DPsites_with_DP=0
            yesspanning_yesmapped_yeshiqual   DP4sites_with_DP4
            yesyesspanning_yesmapped_nilhiqual   yesspanning_yesmapped - yesspanning_yesmapped_yeshiqual
    '''

    CAT_yesspanning_yesmapped_yeshiqual = 0
    CAT_yesspanning_nilmapped = 1
    CAT_yesspanning_yesmapped_nilhiqual = 2
    CAT_nilspanning = 3

    def __init__(self, vcfDP4Path, vcfDPPath, vcfIndelPath, referencePath, mpileupDP, extraHeaders):
        # Open all the input files
        self.vcfDP4 = GormVcf(vcfDP4Path)
        self.vcfDP = GormVcf(vcfDPPath)
        self.vcfIndel = GormVcf(vcfIndelPath)
        self.fasta = FastaFile(referencePath)
        self.mpileupDP = {}
        self.extraHeaders = extraHeaders

        # Merging headers in self.vcfDP:
        for i in itertools.chain(self.vcfDP4._header(), self.vcfIndel._header()):
            self.vcfDP.addHeader(i)

        # indexing mpileup information
        # Actually we are getting the coverage depth and storing it into the
        # array self.mpileupDP (Long)
        a = open(mpileupDP)
        for i in a:
            i = i.split()
            if i[0] not in self.mpileupDP:
                self.mpileupDP[i[0]] = ExtendibleArray("L")
            self.mpileupDP[i[0]].insert(int(i[1]), int(i[3]))

        # Creating a genome map for every category
        # each category has a genome map which indicates
        # wether a specific belongs to this category
        self.yesspanning_nilmapped = GenomeMap()
        self.yesspanning_yesmapped_yeshiqual = GenomeMap()
        self.yesspanning_yesmapped_nilhiqual = GenomeMap()

        # Filling genome maps according to the flags found in records
        # we use indexRecord which is an interator that indexes every
        # record at the same time it iterates over them
        for i in self.vcfDP4.indexRecords():
            if "DP4=" in i[7]:
                self.yesspanning_yesmapped_yeshiqual.addSite(i[0], int(i[1]))

        for i in self.vcfDP.indexRecords():
            pos = int(i[1])
            if "DP=0" in i[7]:
                self.yesspanning_nilmapped.addSite(i[0], pos)
            elif not self.yesspanning_yesmapped_yeshiqual.getSite(i[0], pos):
                self.yesspanning_yesmapped_nilhiqual.addSite(i[0], pos)

        for i in self.vcfIndel.indexRecords():
            # Here we only index records
            pass

    def getCategory(self, contig, pos):
        # Given a site it reports to which category does it belong
        if self.yesspanning_yesmapped_yeshiqual.getSite(contig,
                                                        pos):
                                                            return GormVcfMerger.CAT_yesspanning_yesmapped_yeshiqual
        if self.yesspanning_nilmapped.getSite(contig, pos):
            return GormVcfMerger.CAT_yesspanning_nilmapped
        if self.yesspanning_yesmapped_nilhiqual.getSite(contig,
                                                        pos):
                                                            return GormVcfMerger.CAT_yesspanning_yesmapped_nilhiqual
        return GormVcfMerger.CAT_nilspanning

    def mergeHeaders(self):
        # Generator that returns merged headers
        for i in self.vcfDP._header():
            yield i
        for i in range(len(self.extraHeaders)):
            yield "##mpileupCommand{0}={1}".format(i + 1, self.extraHeaders[i])

    def mergeRecords(self):
        # Generator that returns merged records

        contigs = self.fasta.getChromosomes().items()
        contigs.sort()

        for contig, seqlen in contigs:
            for pos in xrange(1, seqlen + 1):
                category = self.getCategory(contig, pos)
                line = self.newVcfLine(category, contig, pos)
                yield line
                try:
                    indelR = self.vcfIndel.getRecord(contig, pos, indel=True)
                    line = self.newVcfIndelLine(indelR)
                    yield line
                except GormVcf.RecordNotExistent as e:
                    pass

    # newVcfIndelLine and newVcfLine create normalized record without missing
    # fields, assigning default values if needed

    def newVcfIndelLine(self, record):
        output = []
        output += record[:7]
        info = dict([i.split("=") for i in record[7].split(";") if "=" in i])
        info.setdefault("DP", "0")
        info.setdefault("DP4", "0,0,0,0")
        info.setdefault("MQ", "0")
        # ?? don't know apparenty should be prints MQ4 instead of MQ
        info.setdefault("MQ4", info["MQ"])
        output.append("INDEL;DP={0};DP4={1};MQ={2}".format(
            info["DP"], info["DP4"], info["MQ"]))
        output.append("GT:DP")
        formatD = dict(zip(record[8].split(':'), record[9].split(':')))
        formatD.setdefault("GT", "./.")
        formatD.setdefault("DP", "0")
        output.append(formatD["GT"] + ":" + formatD["DP"])

        return "\t".join(output)

    def newVcfLine(self, cat, contig, pos):
        output = []

        if cat == GormVcfMerger.CAT_nilspanning:
            output = [contig, str(pos), '.', self.fasta.getChunk(contig, pos - 1, 1), ".", '999', '.',
                      'BaseCounts=0,0,0,0;DP=0;DP4=0,0,0,0;GC=0;MQ=0;MQ4=0', 'GT:DP', "./.:0"]

        elif cat == GormVcfMerger.CAT_yesspanning_nilmapped:
            dpRecord = self.vcfDP.getRecord(contig, pos)
            output += dpRecord[:7]
            dpInfo = dict([i.split("=")
                           for i in dpRecord[7].split(";") if "=" in i])
            dpInfo.setdefault("GC", "0")
            output.append("BaseCounts=0,0,0,0;DP={0};DP4=0,0,0,0;GC={1};MQ=0;MQ4=0".format(self.mpileupDP[contig][pos],
                                                                                           dpInfo["GC"]))
            output.append("GT:DP")
            formatD = dict(zip(dpRecord[8].split(':'), dpRecord[9].split(':')))
            formatD.setdefault("GT", "./.")
            formatD["DP"] = "0"
            output.append(formatD["GT"] + ":" + formatD["DP"])

        elif cat == GormVcfMerger.CAT_yesspanning_yesmapped_nilhiqual:
            dpRecord = self.vcfDP.getRecord(contig, pos)
            output += dpRecord[:7]
            dpInfo = dict([i.split("=")
                           for i in dpRecord[7].split(";") if "=" in i])
            dpInfo.setdefault("BaseCounts", "0,0,0,0")
            dpInfo.setdefault("MQ", "0")
            dpInfo.setdefault("GC", "0")
            output.append("BaseCounts={0};DP={1};DP4=0,0,0,0;GC={2};MQ={3};MQ4=0".format(dpInfo["BaseCounts"],
                                                                                         self.mpileupDP[contig][pos],
                                                                                         dpInfo["GC"], dpInfo["MQ"]))
            output.append("GT:DP")
            formatD = dict(zip(dpRecord[8].split(':'), dpRecord[9].split(':')))
            formatD.setdefault("GT", "./.")
            formatD["DP"] = "0"
            output.append(formatD["GT"] + ":" + formatD["DP"])

        elif cat == GormVcfMerger.CAT_yesspanning_yesmapped_yeshiqual:
            dpRecord = self.vcfDP.getRecord(contig, pos)
            dp4Record = self.vcfDP4.getRecord(contig, pos)
            output += dp4Record[:7]
            dp4Info = dict([i.split("=")
                            for i in dp4Record[7].split(";") if "=" in i])
            dpInfo = dict([i.split("=")
                           for i in dpRecord[7].split(";") if "=" in i])
            dp4Info.setdefault("BaseCounts", "0,0,0,0")
            dp4Info.setdefault("GC", "0")
            dp4Info.setdefault("DP4", "0,0,0,0")
            dp4Info.setdefault("MQ", "0")
            dpInfo.setdefault("MQ", "0")
            output.append("BaseCounts={0};DP={1};DP4={2};GC={3};MQ={4};MQ4={5}".format(dp4Info["BaseCounts"],
                                                                                       self.mpileupDP[contig][pos],
                                                                                       dp4Info["DP4"], dp4Info["GC"],
                                                                                       dpInfo["MQ"], dp4Info["MQ"]))
            output.append("GT:DP")
            formatD = dict(
                zip(dp4Record[8].split(':'), dp4Record[9].split(':')))
            formatD.setdefault("GT", "./.")
            formatD.setdefault("DP", "0")
            output.append(formatD["GT"] + ":" + formatD["DP"])

        return "\t".join(output)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Merge VCF files (HighQ depth + some depth + Indels)')
    parser.add_argument("-vcfhq", dest="vcfhq", required=True,
                        help="Vcf for sites containing high quality depth (varannot.vcf)")
    parser.add_argument("-vcf", dest="vcf", required=True,
                        help="Vcf for sites containing some depth (varannotDP.vcf)")
    parser.add_argument("-vcfindel", dest="vcfindel",
                        required=True, help="Vcf Containig indels (indels.vcf)")
    parser.add_argument("-mpileup", dest="mpileup", required=True,
                        help="output.bam of samtools mpileup (mpileupDP.txt)")
    parser.add_argument("-ref", dest="ref", required=True,
                        help="Reference used for the alignment (fasta)")
    parser.add_argument("-o", dest="output.bam", required=True,
                        help="Output Vcf file")
    args = parser.parse_args()

    a = GormVcfMerger(args.vcfhq, args.vcf, args.vcfindel,
                      args.ref, args.mpileup)
    out = open(args.output, "w")
    for i in a.mergeHeaders():
        out.write(i + "\n")
    for i in a.mergeRecords():
        out.write(i + "\n")
