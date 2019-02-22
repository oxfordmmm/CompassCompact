#!/usr/bin/env python

# ============================================================================ #
# g4_annotvcf.py - Reimplementation of gorm_annotvcf.py
# Main changes: Data sctructures, OOP
#
# Execution of same sample in Bespin [AMD Opteron(TM) Processor 6274 - 32 Core - 256 Gb Ram]
#   gorm_mergevcfs.py: 29m - 15Gb Ram
#         g4_merge.py: 26m - 692Mb Ram   Improvements: 11% and 96% respectively
#
# ============================================================================ #
# Camilla Ip
# camilla.ip@stats.ox.ac.uk
# April 2012
#
# Carlos del Ojo Elias
# carlos.delojoelias@ndm.ox.ac.uk
# July 2013
# ============================================================================ #


VERSION = "1.0.1"

#
#  TODO variants near variant/indel
#


"""
Add extra INFO to VCF file

Add position-specific INFO fields needed to decide if a site is callable.
These are statistics that were not available from third-party software.

The new fields added are listed alphabetically below.

ID           Number  Type    Description

ABQ4         1       Float   Average PHRED-scaled quality of the BASES4 high-quality bases from the pileup
BASES4       1       String  The high-quality bases in the pileup that contribute to high-quality depth in DP4
BaseCounts4  4       String  Number of ACGTs contributing to high-quality depth in DP4 (or NA if could not be computed from the mpileup string)
BCALL        1       String  Base call that would be made by standard Gorm vN.N pipeline
DZ           1       Float   Z-score (number of stddevs from mean) of absolute depth DP (where mean=FLOAT and stddev=FLOAT for sample)
DZ4          1       Float   Z-score (number of stddevs(=FLOAT) from mean(=FLOAT) of high-quality depth)
							 Z-score (number of stddevs from mean) of high-quality depth DP4 (where mean=FLOAT and stddev=FLOAT for sample)
DZ4L         1       Float   Studentized residuals after removing GC bias
DM           1       Float   Robust Z-score (number of scaled MAD=FLOAT from median=FLOAT of absolute depth)
DM4          1       Float   Robust Z-score (number of scaled MAD=FLOAT from median=FLOAT of high-quality depth)
DM4L         1       Float   MAD (Median Absolute Deviation, number of medians from mean) of GC-corrected high-quality depth DPT4L (where mean=FLOAT and median=FLOAT)
DPT4L        1       Float   Raw residuals after removing GC bias: residual=sum(DP4)-(slope*GC)-intercept
IqFLOATwINT  1       String  Number of indels with QUAL>=FLOAT in a window of +/- INT bp
PCALL        1       Float   Proportion of the all BaseCounts bases that support the base call that would be made by the standard Gorm pipeline (zero if not called)
PCONS        1       Float   Proportion of the all BaseCounts bases that support the GT genotype base(s)
PCALL4       1       Float   Proportion of high-quality DP4 bases that support the Gorm pipeline base call (1.0 if not called)
PCONS4       1       Float   Proportion of high-quality DP4 bases that support the GT genotype base(s)
QUALS4       1       String  The PHRED-scaled qualities for bases in the pileup that contribute to high-quality depth in DP4
SBR          1       Integer Self-self BLAST repeat-region with 1 meaning in a repeat region and 0 otherwise
VqFLOATwINT  1       String  Number of variants with QUAL>=FLOAT in a window of +/- INT bp

Notes:

Note that it not possible to routinely print the QUALS4 field because characters
'"' (decimal 34) and ';' (decimal 59) are valid base quality figures in the
33-offset system and VCF format cannot cope with either of these characters inside
the INFO var=val tokens. The QUALS4 field should thus only be printed when running
the program to debug some other problem.

"""

import array
import itertools
import math
import numpy
import os
import struct
import sys
import argparse
from lib.extendibleArray import ExtendibleArray
from lib.fastafile import FastaFile
from lib.genomemap import GenomeMap
from lib.gormVCF import GormVcf
from lib.logerror import LE, dump_exc
from lib.mpileup import PileupReader

class VcfAnnotator:
    ######################################################
    #
    #  STATIC METHODS FOR CALCULATIONS
    #
    # 0
    # We do I implement Mean and Median? easy, when numpy deals with lists of millions of
    # numbers the memory taken by the software increases in a few hundreds...
    # we don't want that, no nO NO!

    @staticmethod
    def Mean(x):

        k = len(x)
        if not k:
            return 0.0
        return float(sum(x)) / k

    @staticmethod
    def Median(x):
        return x[len(x) / 2 + 1] if len(x) % 2 else (x[len(x) / 2 + 1] + x[len(x) / 2]) / 2

    @staticmethod
    def linreg_coeff(x_list, y_list):
        """
        Use numpy to compute and return the b and a coefficients for the linear regression.
        """
        if not len(x_list) or not len(y_list):
            return 0.0, 0.0
        b, a = list(numpy.linalg.lstsq(numpy.vstack(
            [x_list, numpy.ones(len(x_list))]).T, y_list)[0])
        # Slope, Intercept
        return b, a

    @staticmethod
    def stdev(a):
        mean = VcfAnnotator.Mean(a)
        res = float(0)
        for i in a:
            res += (i - mean) ** 2

        return math.sqrt(res / len(a))

    @staticmethod
    def Compute_MAD(a, c=0.6745, axis=None):
        """
        Median Absolute Deviation along given axis of an array:

        median(abs(a - median(a))) / c

        c = 0.6745 is the constant to convert from MAD to std; it is used by
        default

        Copied from http://code.google.com/p/agpy/source/browse/trunk/agpy/mad.py
        Downloaded 7-Dec-2012.
        """

        LE.debug("Computing MAD for {0} #elements ({1}bytes)".format(
            len(a), sys.getsizeof(a)))

        d = VcfAnnotator.Median(a)
        nary = array.array('d')
        for i in a:
            nary.append(float(abs(i - d)) / c)
        nary = list(nary)
        nary.sort()
        mad = VcfAnnotator.Median(nary)
        del nary

        return mad

    ############################################
    # END CALCULATION METHODS
    ############################################

    def __init__(self, vcf, pileup, repmask, annot_basecall, annot_selfblastregions, annot_hqdepthinfo,
                 annot_lgcdepthinfo, annot_absdepthinfo, annot_bases4, annot_quals4):
        self.ANNOT_BASECALL = annot_basecall
        self.ANNOT_SELFBLASTREGIONS = annot_selfblastregions
        self.ANNOT_HQDEPTHINFO = annot_hqdepthinfo
        self.ANNOT_LGCDEPTHINFO = annot_lgcdepthinfo

        self.ANNOT_ABSDEPTHINFO = annot_absdepthinfo
        self.ANNOT_BASES4 = annot_bases4
        self.ANNOT_QUALS4 = annot_quals4
        
        self.vcf = GormVcf(vcf)
        self.pileup = PileupReader(pileup, index=True)
        self.repmask = GenomeMap()

        fastarepmask = FastaFile(repmask)
        # fastarepmask=FastaFile(mask)
        for ch, lgth in fastarepmask.getChromosomes().items():
            for i in xrange(lgth):
                if fastarepmask.getChunk(ch, i, 1) == "1":
                    # We add 1 because the genome itself is 1-coordinate
                    self.repmask.addSite(ch, i + 1)

        self.DPavg = self.DP4avg = self.DP4Lavg = 0.0
        self.DPstd = self.DP4std = self.DP4Lstd = 1.0
        self.DPmedian = self.DP4median = self.DP4Lmedian = 0.0
        self.DPmad = self.DP4mad = self.DP4Lmad = 1

        self.lm_slope = self.lm_intercept = None

        self.DP4s = self.DP4Ls = None

        self.calculateStats()

        self.addNewHeaders()

    def addNewHeaders(self):
        self.vcf.addHeader(
            '##INFO=<ID=BCALL,Number=1,Type=String,Description="Base call or indel call that would be made by standard g4_annotvcf {versionnumber} pipeline">'.format(
                versionnumber=VERSION))
        self.vcf.addHeader(
            '##INFO=<ID=PCONS,Number=1,Type=Float,Description="Proportion of the all BaseCounts bases that support the GT genotype base(s)">')
        self.vcf.addHeader(
            '##INFO=<ID=PCALL,Number=1,Type=Float,Description="Proportion of the all BaseCounts bases that support the base call that would be made by the standard Gorm pipeline (zero if not called)">')
        self.vcf.addHeader(
            '##INFO=<ID=DZ,Number=1,Type=Float,Description="Z-score (number of stddevs from mean) of absolute depth DP (where mean={0:.3f} and stddev={1:.3f} for sample)">'.format(
                self.DPavg, self.DPstd))
        self.vcf.addHeader(
            '##INFO=<ID=DM,Number=1,Type=Float,Description="Robust Z-score (number of scaled MAD={0:.3f} from median={1:.3f} of absolute depth)">'.format(
                self.DPmad, self.DPmedian))
        self.vcf.addHeader(
            '##INFO=<ID=BaseCounts4,Number=4,Type=String,Description="Number of ACGTs contributing to high-quality depth in DP4 or or NA if could not be computed from the mpileup string">')
        self.vcf.addHeader(
            '##INFO=<ID=ABQ4,Number=1,Type=Float,Description="Average PHRED-scaled quality of the BASES4 high-quality bases from the pileup">')
        self.vcf.addHeader(
            '##INFO=<ID=PCONS4,Number=1,Type=Float,Description="Proportion of high-quality DP4 bases that support the GT genotype base(s)">')
        self.vcf.addHeader(
            '##INFO=<ID=PCALL4,Number=1,Type=Float,Description="Proportion of high-quality DP4 bases that support the Gorm pipeline base call (1.0 if not called)">')
        self.vcf.addHeader(
            '##INFO=<ID=DZ4,Number=1,Type=Float,Description="Z-score (number of stddevs={0:.3f} from mean={1:.3f} of high-quality depth)">'.format(
                self.DP4std, self.DP4avg))
        self.vcf.addHeader(
            '##INFO=<ID=DM4,Number=1,Type=Float,Description="Robust Z-score (number of scaled MAD={0:.3f} from median={1:.3f} of high-quality depth)">'.format(
                self.DP4mad, self.DP4median))
        self.vcf.addHeader(
            '##INFO=<ID=BASES4,Number=1,Type=String,Description="The high-quality bases in the pileup that contribute to high-quality depth in DP4">')
        self.vcf.addHeader(
            '##INFO=<ID=QUALS4,Number=1,Type=String,Description="The PHRED-scaled qualities for bases in the pileup that contribute to high-quality depth in DP4">')
        self.vcf.addHeader(
            '##INFO=<ID=DPT4L,Number=1,Type=Float,Description="Raw residuals after removing GC bias: residual=sum(DP4)-({0:.3f}*GC)-{1:.3f}">'.format(
                self.lm_slope, self.lm_intercept))
        self.vcf.addHeader(
            '##INFO=<ID=DZ4L,Number=1,Type=Float,Description="Studentized residuals after removing GC bias (mean={0:.3f} and stddev={1:.3f})">'.format(
                self.DP4Lstd, self.DP4Lavg))
        self.vcf.addHeader(
            '##INFO=<ID=DM4L,Number=1,Type=Float,Description="Robust Z-score (number of scaled MAD={0:.3f} from median={1:.3f} of GC-corrected high-quality depth)">'.format(
                self.DP4Lmad, self.DP4Lmedian))
        self.vcf.addHeader(
            '##INFO=<ID=SBR,Number=1,Type=Integer,Description="Self-self BLAST repeat-region with 1 meaning in a repeat region and 0 otherwise">')

    def calculateStats(self):
        '''
            Calculates stats (stdev,mean,median and MAD) for:
            DP: Raw read depth, the absolute depth that includes all covering reads,
                including those with a gap in the read at this site
            DP4: # high-quality ref-forward bases, ref-reverse, alt-forward and alt-reverse bases
            slope and intercept: Linear regression coefficients for DP4 vs GC
        '''

        # We extract all GC values from the VCF
        GCs = self.vcf.getValue("GC=([0-9\.]+).*", 'd')

        # We extract all DP values from the VCF amd filter them (only positive
        # ones)
        DPs = self.vcf.getValue("DP=([0-9]+).*", 'i')
        filtered_DPs = []
        for i in DPs.values():
            filtered_DPs += [j for j in i if j > 0]
        filtered_DPs.sort()

        # Calculate all the stats on DP values
        if filtered_DPs:
            self.DPavg = self.Mean(filtered_DPs)
            self.DPmedian = VcfAnnotator.Median(filtered_DPs)
            self.DPstd = VcfAnnotator.stdev(filtered_DPs)
            self.DPmad = VcfAnnotator.Compute_MAD(filtered_DPs)

        # We extract all DP4 values from the VCF, they are 4 integers comma
        # separated (we sum them up)
        DP4s = self.vcf.getValue(
            "DP4=([0-9,]+).*", 'i', mod=lambda x: sum([int(i) for i in x.split(",")]))

        # We create two arrays (integers for DP4 and floats for GC) to store those
        # values in positions where DP4 is positive

        filtered_DP4s = array.array("i")
        filtered_GCs = array.array('d')
        contigs = DP4s.keys()
        for i in contigs:
            for j in range(len(DP4s[i])):
                if DP4s[i][j] > 0:
                    filtered_DP4s.append(DP4s[i][j])
                    filtered_GCs.append(GCs[i][j])

        # We calculate the linear regression coeficients for these points
        self.lm_slope, self.lm_intercept = VcfAnnotator.linreg_coeff(
            filtered_GCs, filtered_DP4s)

        # Then we calculate same statistics for DP4
        filtered_DP4s = list(filtered_DP4s)
        filtered_DP4s.sort()

        if filtered_DP4s:
            self.DP4avg = VcfAnnotator.Mean(filtered_DP4s)
            self.DP4median = VcfAnnotator.Median(filtered_DP4s)
            self.DP4std = VcfAnnotator.stdev(filtered_DP4s)
            self.DP4mad = VcfAnnotator.Compute_MAD(filtered_DP4s)

        # DP4L contains Raw residuals after removing GC bias
        DP4Ls = {}

        # we fill DP4L with new calculations
        for i in contigs:
            DP4Ls[i] = ExtendibleArray('d')
            for j in range(len(DP4s[i])):
                if not DP4s[i][j]:
                    DP4Ls[i].insert(j, 0.0)
                else:
                    DP4Ls[i].insert(
                        j, DP4s[i][j] - (self.lm_slope * GCs[i][j]) - self.lm_intercept)

        # We filter all values inside DP4L and filter out the positve ones only
        filtered_DP4Ls = []
        for i in DP4Ls.values():
            filtered_DP4Ls += [j for j in i if j > 0]
        filtered_DP4Ls.sort()

        # And calculate same stats
        if filtered_DP4Ls:
            self.DP4Lavg = VcfAnnotator.Mean(filtered_DP4Ls)
            self.DP4Lmedian = VcfAnnotator.Median(filtered_DP4Ls)
            self.DP4Lstd = VcfAnnotator.stdev(filtered_DP4Ls)
            self.DP4Lmad = VcfAnnotator.Compute_MAD(filtered_DP4Ls)

        self.DP4s = DP4s
        self.DP4Ls = DP4Ls


    ################ TEMPORAL METHODS TO STORE AND LOAD STATS WITHOUT CALCULAT

    def storeStats(self, fil):
        import pickle
        a = open(fil, "w")
        pickle.dump(((self.DPavg, self.DP4avg, self.DP4Lavg), (self.DPstd, self.DP4std, self.DP4Lstd),
                     (self.DPmedian, self.DP4median,
                      self.DP4Lmedian), (self.DPmad, self.DP4mad, self.DP4Lmad),
                     (self.lm_slope, self.lm_intercept), (self.DP4s, self.DP4Ls)), a)
        a.close()

    def loadStats(self, fil):
        import cPickle
        a = open(fil)
        ((self.DPavg, self.DP4avg, self.DP4Lavg), (self.DPstd, self.DP4std, self.DP4Lstd),
         (self.DPmedian, self.DP4median,
          self.DP4Lmedian), (self.DPmad, self.DP4mad, self.DP4Lmad),
         (self.lm_slope, self.lm_intercept), (self.DP4s, self.DP4Ls)) = cPickle.load(a)
        a.close()

    def Headers(self):
        for i in self.vcf._header():
            yield i
        yield "##g4_vcfannot={0}".format(" ".join(sys.argv))

    def Annotate(self):
        '''
            Methods that returns a generator to iterate over all the anottated records
            Every anottated record will be a string
        '''

        for i in self.vcf.indexRecords():
            pass

        for contig, pos in self.vcf.iterContigsPos():
            (snp, indel, baseCall, indelCall, bases4, quals4, basecounts4, pileupinfo) = self.vcf.getAdditionalInfo(
                contig, pos, self.pileup)

            if snp:
                info = snp[GormVcf.fld_PARSED_INFO]
                # TODO variants near variant/indel
                # -> if ANNOT_VARNEARINDSPEC:
                # -> if ANNOT_VARNEARVARSPEC:
                ####
                if self.ANNOT_ABSDEPTHINFO:
                    self.annot_absdepthinfo(snp, baseCall, pileupinfo)

                if self.ANNOT_HQDEPTHINFO:
                    self.annot_hqdepthinfo(
                        snp, contig, pos, baseCall, basecounts4, bases4, quals4)

                if self.ANNOT_BASECALL:
                    info["BCALL"] = baseCall

                if self.ANNOT_LGCDEPTHINFO:
                    self.annot_lgcdepthinfo(snp, contig, pos, pileupinfo)

                if self.ANNOT_SELFBLASTREGIONS:
                    masked = self.repmask.getSite(contig, pos)
                    if masked:
                        info["SBR"] = 1
                    else:
                        info["SBR"] = 0

                yield self.dump(snp)

            if indel:
                info = indel[GormVcf.fld_PARSED_INFO]
                if self.ANNOT_BASECALL:
                    if indelCall:
                        info["BCALL"] = indelCall
                    else:
                        info["BCALL"] = GormVcf.APPROXCOVZDEPTHCHAR

                if self.ANNOT_SELFBLASTREGIONS:
                    masked = self.repmask.getSite(contig, pos)
                    if masked:
                        info["SBR"] = 1
                    else:
                        info["SBR"] = 0

                yield self.dump(indel, indel=True)

    ######################## DIFFERENT ANNOTATION MODULES ####################

    def annot_absdepthinfo(self, snp, baseCall, pileupinfo):
        info = snp[GormVcf.fld_PARSED_INFO]

        bcsum = sum(info["BaseCounts"])
        info["DZ"] = info["DM"] = 0.0

        if not snp[GormVcf.fld_PARSED_DATA]["GTBASES"] or not bcsum:
            info["PCONS"] = 0.0
        elif snp[GormVcf.fld_PARSED_DATA]["GTBASES"] in [GormVcf.UNCALLABLECHAR, GormVcf.APPROXCOVZDEPTHCHAR]:
            info["PCONS"] = 1.0
        else:
            info["PCONS"] = snp[GormVcf.fld_PARSED_DATA]["GTDEPTH"] / \
                float(bcsum)

        if not baseCall or not bcsum or not pileupinfo:
            info["PCALL"] = 0.0
        elif baseCall in [GormVcf.UNCALLABLECHAR, GormVcf.APPROXCOVZDEPTHCHAR]:
            info["PCALL"] = 1.0
        else:
            info["PCALL"] = snp[GormVcf.fld_PARSED_DATA]["BCALLDEPTH"] / \
                float(bcsum)

        if self.DPstd:
            info["DZ"] = (info["DP"] - self.DPavg) / self.DPstd

        if self.DPmad:
            info["DM"] = (info["DP"] - self.DPmedian) / self.DPmad

    def annot_hqdepthinfo(self, snp, contig, pos, baseCall, basecounts4, bases4, quals4):
        info = snp[GormVcf.fld_PARSED_INFO]
        bc4sum = sum(basecounts4)

        info["BaseCounts4"] = ",".join([str(i) for i in basecounts4])
        info["ABQ4"] = info["PCALL4"] = info["DZ4"] = info["DM4"] = 0.0

        if quals4:
            if self.ANNOT_QUALS4:
                info["QUALS4"] = struct.pack(
                    "b" * len(quals4), *[i + 33 for i in quals4])
            info["ABQ4"] = self.Mean(quals4)

        if bases4:
            if self.ANNOT_BASES4:
                info["BASES4"] = "".join([str(i) for i in bases4])
            info["PCALL4"] = float(
                len([i for i in bases4 if i.upper() == baseCall])) / len(bases4)

        if not snp[GormVcf.fld_PARSED_DATA]["GTBASES"] or not bc4sum:
            info["PCONS4"] = 0.0
        elif snp[GormVcf.fld_PARSED_DATA]["GTBASES"] in [GormVcf.UNCALLABLECHAR, GormVcf.APPROXCOVZDEPTHCHAR]:
            info["PCONS4"] = 1.0
        else:
            info["PCONS4"] = snp[GormVcf.fld_PARSED_DATA]["GTHQDEPTH"] / \
                float(bc4sum)

        if self.DP4std:
            info["DZ4"] = float(self.DP4s[contig][pos] -
                                self.DP4avg) / self.DP4std

        if self.DP4mad:
            info["DM4"] = float(self.DP4s[contig][pos] -
                                self.DP4median) / self.DP4mad

    def annot_lgcdepthinfo(self, snp, contig, pos, pileupinfo):
        info = snp[GormVcf.fld_PARSED_INFO]

        info["DPT4L"] = info["DZ4L"] = info["DM4L"] = 0.0
        if self.DP4Ls[contig].is_set(pos) and pileupinfo:
            DP4L = self.DP4Ls[contig][pos]
            info["DPT4L"] = DP4L

            if self.DP4Lstd:
                info["DZ4L"] = (DP4L - self.DP4Lavg) / self.DP4Lstd
            if self.DP4Lmad:
                info["DM4L"] = (DP4L - self.DP4Lmedian) / self.DP4Lmad

    ############### METHOD THAT RETURNS RECORDS IN STRING FORMAT ###########

    def dump(self, record, indel=False):
        newinfo = record[GormVcf.fld_PARSED_INFO]
        Q4 = None
        if "QUALS4" in newinfo:
            Q4 = newinfo["QUALS4"]
            del newinfo["QUALS4"]

        newinfo2 = []

        for i, j in newinfo.items():
            if type(j) == float:
                j = "{0:.3f}".format(j)
                newinfo2.append((i, j))

            elif type(j) == list:
                newinfo2.append((i, ",".join([str(k) for k in j])))
            else:
                newinfo2.append((i, str(j)))

        newinfo2.sort()
        if Q4:
            newinfo2 += (("QUALS4", Q4),)

        newinfo2 = ";".join(["=".join(i) for i in newinfo2])

        if indel:
            newinfo2 = "INDEL;" + newinfo2

        record[GormVcf.fld_INFO] = newinfo2
        record[5] = "{0:.2f}".format(float(record[5]))

        data = record[GormVcf.fld_PARSED_DATA]
        fmt = record[GormVcf.fld_FORMAT].split(":")
        record[GormVcf.fld_DATA] = ":".join([data[i] for i in fmt])

        return ("\t".join(record[:10]))


################### END VcfAnnotator #########


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Calculate extra stats for a VCF file (previously merged with extra information)')
    parser.add_argument("-vcf", dest="vcf", required=True,
                        help="Vcf input file")
    parser.add_argument("-mpileup", dest="mpileup",
                        required=True, help="Mpileup File")
    parser.add_argument("-refmask", dest="ref_id", required=True,
                        help="Masked reference (for repeated regions) used for the alignment (fasta)")  # type=lambda x: parser.error("File does not exist [{0}]".format(x)) if not os.path.isfile(x) else x) #dest="refmask"
    # parser.add_argument("-u", dest="guuid", required=True)
    parser.add_argument("-basecall", dest="basecall", help="Annotate basecall [BCALL]", default=False,
                        action='store_true')
    parser.add_argument("-selfblastR", dest="selfblastR", help="Annotate Self-self BLAST repeat-regions [SBR]",
                        default=False, action='store_true')
    parser.add_argument("-hqdepthinfo", dest="hqdepthinfo",
                        help="Annotate HQ depth info fields [PCONS4,PCALL4,DM4,DZ4]", default=False,
                        action='store_true')
    parser.add_argument("-lgcdepthinfo", dest="lgcdepthinfo",
                        help="Annotate raw residuals after removing GC bias [DPT4L,DZ4L,DM4L]", default=False,
                        action='store_true')
    parser.add_argument("-absdepthinfo", dest="absdepthinfo",
                        help="Annotate absolute depth info fields [PCONS,PCALL,DM,DZ]", default=False,
                        action='store_true')
    parser.add_argument("-bases4", dest="bases4",
                        help="Annotate high-quality bases in the pileup that contribute to high-quality depth in DP4 [BASES4]",
                        default=False, action='store_true')
    parser.add_argument("-quals4", dest="quals4",
                        help="Annotate PHRED-scaled qualities for bases in the pileup that contribute to high-quality depth in DP4 [QUALS4]",
                        default=False, action='store_true')

    parser.add_argument("-o", dest="output", required=True,
                        help="Output Vcf file")
    args = parser.parse_args()

    try:
        a = VcfAnnotator(args.vcf, args.mpileup, args.ref_id, annot_basecall=args.basecall,
                         annot_selfblastregions=args.selfblastR, annot_hqdepthinfo=args.hqdepthinfo,
                         annot_lgcdepthinfo=args.lgcdepthinfo, annot_absdepthinfo=args.absdepthinfo,
                         annot_bases4=args.bases4, annot_quals4=args.quals4)
        out = open(args.output, "w")
        for i in itertools.chain(a.Headers(), a.Annotate()):
            out.write(i + "\n")
        out.close()
    except:
        dump_exc()
