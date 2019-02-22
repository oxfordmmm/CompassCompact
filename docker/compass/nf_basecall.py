#!/usr/bin/env python
import array
from collections import Counter
from lib.gormVCF import GormVcf
import argparse
import sys
import gzip
from lib.extendibleArray import ExtendibleArray
import re
from lib.logerror import LE, dump_exc
import os

VERSION = "1.0.1"
APPROXCOVZDEPTHCOUNT = 0
UNCALLABLECHAR = 'N'
APPROXCOVZDEPTHCHAR = "-"


class MappingStats:
    class VectorData(object):
        def __init__(self, typ, data):
            self.typ = typ
            if type(data) == str:
                self.data = array.array(self.typ)
                self.data.fromstring(data)
            elif type(data) == list:
                self.data = array.array(self.typ, data)
            else:
                self.data = data

            self.ceil = 0

        def ceilMode(self):
            self.ceil = Counter(self.data).most_common()
            if self.ceil[0][0]:
                self.ceil = self.ceil[0][0]
            else:
                self.ceil.pop(0)
                if not self.ceil:
                    self.ceil = 1
                else:
                    self.ceil = self.ceil[0][0]

        def ceilAvg(self):
            if not self.data:
                return 0
            self.ceil = sum(self.data) / len(self.data)

        def ceilMax(self):
            if not self.data:
                return 0
            self.ceil = max(self.data)

        def vectorSummary(self, l, normalized=True, binary=False):
            """Normalized will transform data to values from 0 to 1
               binary will interpret original positions as 0 or 1, 1 meaning different than 0 (so any value)
               the binary flag will convert the vector to 1's and 0's and then averages will be performed 
               leading to values between 0 and 1 """

            if binary:
                normalized = False

            res = [i for i in (float(sum(i)) / len(i)
                               for i in self.chunks(self.data, l, binary))]
            if not normalized:
                return res
            return [min(float(i) / self.ceil, 1) for i in res]

        def toList(self):
            return [self.typ, self.data.tostring(), self.ceil]

        @staticmethod
        def chunks(l, n, binary):
            """Yield successive n chunks from l."""
            chunksize = float(len(l)) / n
            intchunk = int(chunksize)
            if not chunksize.is_integer():
                intchunk += 1

            if binary:
                for i in range(n):
                    pos = int(chunksize * i)
                    yield [1 if j else 0 for j in l[pos:pos + intchunk]]
            else:
                for i in range(n):
                    pos = int(chunksize * i)
                    yield l[pos:pos + intchunk]

        @staticmethod
        def fromList(l):
            typ, data, ceil = l
            v = MappingStats.VectorData(typ, data)
            v.ceil = ceil
            return v

        def mean(self):
            return sum(self.data) / float(len(self.data))

        def median(self):
            lst = self.data
            lst = sorted(lst)
            if len(lst) < 1:
                return None
            if len(lst) % 2 == 1:
                return lst[((len(lst) + 1) / 2) - 1]
            else:
                return float(sum(lst[(len(lst) / 2) - 1:(len(lst) / 2) + 1])) / 2.0

        def mode(self):
            return Counter(self.data).most_common(1)[0][0]

        def __gt__(self, val):
            assert type(val) in [int, float]
            lst = self.data
            vals = [i for i in self.data if i > val]
            return float(len(vals)) / len(lst) * 100

        def __lt__(self, val):
            assert type(val) in [int, float]
            lst = self.data
            vals = [i for i in self.data if i < val]
            return float(len(vals)) / len(lst) * 100

    def __init__(self, sampleid="NA"):
        self.data = {'coverages': {}, 'mapq': {}, 'highcov': {},
                     'called': {}, 'mutated': {}, 'hz': {}}
        self.ID = sampleid

    def covStats(self):
        # ofile.write("ID\tREF\tmean\tmedian\tmode\tcov1\tcov5\tcov10\n")
        res = {}
        for ch, data in self.data["coverages"].items():
            res[ch] = [data.mean(), data.median(), data.mode(),
                       data > 0, data > 4, data > 9]
        # ofile.close()
        return res

    def varStats(self, ofile):
        pass

    def processVcf(self, vcf):
        if self.ID == "NA":
            self.ID = os.path.basename(vcf)
        self.vcf = gzip.GzipFile(vcf)
        self.readVcf()

    def vSum(self, dtype, chromosome, size, norm=True, binary=False):
        """Normalized will transform data to values from 0 to 1
           binary will interpret original positions as 0 or 1, 1 meaning different than 0 (so any value)
           the binary flag will convert the vector to 1's and 0's and then averages will be performed 
           leading to values between 0 and 1, this is useful for when you are not interested in values but 
           just PRESENCE of value"""

        assert dtype in self.data and chromosome in self.data[dtype]
        return self.data[dtype][chromosome].vectorSummary(size, norm, binary)

    def switchData(self):
        for datatype, data in self.data.items():
            for chromosome in data.keys():
                if type(data[chromosome]) != list:
                    data[chromosome] = data[chromosome].toList()
                else:
                    data[chromosome] = MappingStats.VectorData.fromList(
                        data[chromosome])

    def chromosomes(self):
        return set(self.data["coverages"].keys() + self.data["mapq"].keys() + self.data["highcov"].keys() + self.data[
            "called"].keys())

    @staticmethod
    def zeros(n):
        for i in xrange(n):
            yield 0

    @staticmethod
    def mapq2array(d, length):
        buckets = [(pos, sum(l) / len(l)) for pos, l in d.items()]
        data = array.array('b', MappingStats.zeros(length))
        for pos, val in buckets:
            for i in xrange(pos, min(pos + 100, length)):
                data[i] = val

        return data

    def extractField(self, f, string, sep):
        fpos = string.index(f) + len(f)
        endpos = string.find(sep, fpos)
        return string[fpos:endpos]

    def readVcf(self):
        # return
        # set(self.data["coverages"].keys()+self.data["mapq"].keys()+self.data["highcov"].keys()+self.data["called"].keys())
        for i in self.vcf:
            if i[0] == "#":
                continue
            ch = i[:i.index('\t')]

            if ch not in self.data["coverages"]:
                self.data["coverages"][ch] = MappingStats.VectorData('H', [])
                self.data["mapq"][ch] = MappingStats.VectorData('b', [])
                self.data["highcov"][ch] = MappingStats.VectorData('H', [])
                self.data["called"][ch] = MappingStats.VectorData('b', [])
                self.data["mutated"][ch] = MappingStats.VectorData('b', [])
                self.data["hz"][ch] = MappingStats.VectorData('b', [])

            bc = self.extractField("BaseCounts=", i, ';')
            bc4 = self.extractField("BaseCounts4=", i, ';')
            bc4 = sum([int(j) for j in bc4.split(",")])
            self.data["highcov"][ch].data.append(bc4)
            bc = sum([int(j) for j in bc.split(",")])
            self.data["coverages"][ch].data.append(bc)
            mq = int(self.extractField("MQ=", i, ';'))
            self.data["mapq"][ch].data.append(mq)
            bcall = self.extractField("BCALL=", i, ';')
            if "PASS" in i and bcall in "ACTG":
                self.data["called"][ch].data.append(1)
                called = True
            else:
                self.data["called"][ch].data.append(0)
                called = False
            if 'z\t' in i:
                self.data["hz"][ch].data.append(1)
            else:
                self.data["hz"][ch].data.append(0)

        for i in self.data["coverages"].values():
            i.ceilMode()
        for i in self.data["mapq"].values():
            i.ceilMax()
        for i in self.data["highcov"].values():
            i.ceilMode()
        for i in self.data["mapq"].values():
            i.ceilMax()
        for i in self.data["called"].values():
            i.ceilMax()
        for i in self.data["mutated"].values():
            i.ceilMax()


class ParameterError(Exception):
    pass


class VcfFilter:
    def __init__(self):
        self.AVAIL_FILTERS = []
        self.FILTER_ON = []

    def addFilter(self,parameter_flag,parameter_name,par_type,vcf_par_id,par_description,vcf_par_description,filter_function):
        self.AVAIL_FILTERS.append([parameter_flag,parameter_name,par_type,vcf_par_id,par_description,vcf_par_description,filter_function])

    def setUpFilter(self, parameter_name, value):
        for i in self.AVAIL_FILTERS:
            if i[1] == parameter_name:
                newfilter = i[:]

                if newfilter[2] == bool and value == True:
                    newfilter[2] = True
                    newfilter[3] = newfilter[3].format(value)
                    self.FILTER_ON.append(newfilter)
                elif i[2] != bool:
                    try:
                        newfilter[3] = newfilter[3].format(value)
                        newfilter[5] = newfilter[5].format(value)
                        newfilter[2] = newfilter[2](value)
                    except:
                        raise ParameterError("Error, parameter {0} should be {1} [ WRONG: {2} ]".format(i[1],i[2],value))
                    self.FILTER_ON.append(newfilter)
                break

    def availableFilters(self):
        for i in self.AVAIL_FILTERS:
            yield i

    def writeHeaders(self, outvcf):
        for i in self.FILTER_ON:
            vcfFile.addHeader(i[5])

        vcfFile.addHeader(
            '''##FILTER=<ID=PASS,Number=0,Type=None,Description="Variant or indel has passed all filters">''')
        vcfFile.addHeader('''##g4_basecall={0}'''.format(" ".join(sys.argv)))

        for i in vcfFile._header():
            outvcf.write(i + "\n")

        outvcf.write(
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tEXTRA\n")

    # def filterVcf(self, vcfFile, outvcf_path, outvcfIndel=None, outfasta=None, stats=None, guuid=None, refid=None):
    def filterVcf(self, vcfFile, outvcf_path, outvcfIndel=None, outfasta=None, stats=None, guuid = None, refid = None):
        
        for i in vcfFile.indexRecords():
            pass
        sampleID = guuid
       
        outvcf = gzip.GzipFile(outvcf_path, "w")
        outvcfIndelS = gzip.GzipFile(outvcfIndel, "w")

        self.writeHeaders(outvcf)
        self.writeHeaders(outvcfIndelS)

        fasta = {}
        statistics = {}

        contigs = vcfFile.getChromosomes().items()
        
        # print contigs
        
        contigs.sort(
            key=lambda x: 0 if not "-" in x[0] else int(x[0].split("-")[1]))

        for ch, pos in contigs:
            fasta[ch] = ExtendibleArray("c", "-")
            statistics[ch] = [0, 0, 0, 0, 0]

            for pos in xrange(1, pos):
                indel = None
                snp = vcfFile.getParsedRecord(ch, pos)
                nopass = []
                if snp[GormVcf.fld_PARSED_INFO]["DP"] > APPROXCOVZDEPTHCOUNT:
                    for i in self.FILTER_ON:
                        if i[6](False, snp, i[2]):
                            nopass.append(i[3])
                if not nopass:
                    nopass.append("PASS")
                nopass.sort()
                snp[GormVcf.fld_FILTER] = ";".join(nopass)

                call = UNCALLABLECHAR
                if snp[GormVcf.fld_FILTER] in ["PASS", ".", ""]:
                    if not snp[GormVcf.fld_PARSED_INFO]["DP"]:
                        call = APPROXCOVZDEPTHCHAR
                    call = snp[GormVcf.fld_PARSED_INFO]["BCALL"]
                    if call == snp[GormVcf.fld_REF]:
                        statistics[ch][1] += 1
                if 'z' in snp[GormVcf.fld_FILTER]:
                    statistics[ch][4] += 1

                fasta[ch].insert(pos, call)

                outvcf.write("\t".join(snp[:10]) + "\n")

                if vcfFile.isindel(ch, pos):
                    nopass = []
                    indel = vcfFile.getParsedRecord(ch, pos, True)
                    if indel[GormVcf.fld_PARSED_INFO]["DP"] > APPROXCOVZDEPTHCOUNT:
                        for i in self.FILTER_ON:
                            if i[6](True, indel, i[2]):
                                nopass.append(i[3])
                    if not nopass:
                        nopass.append("PASS")
                    nopass.sort()
                    indel[GormVcf.fld_FILTER] = ";".join(nopass)
                    outvcfIndelS.write("\t".join(indel[:10]) + "\n")
        outvcf.close()

        ms = MappingStats(guuid)
        ms.processVcf(outvcf_path)

        if outfasta:
            fil = outfasta
            fil = gzip.GzipFile(outfasta, "w")

            fitems = fasta.items()
            fitems.sort(
                key=lambda x: 0 if not "-" in x[0] else int(x[0].split("-")[1]))

            for ch, fa_array in fitems:
                # HP: added sample ID
                fil.write(">{0}_{1}_{2} Software:{3} Version:{4}\n".format(
                    sampleID, refid, ch, __file__, VERSION))

                fa_string = fa_array[1:len(fa_array)].tostring()
                for i in xrange(0, len(fa_string), 100):
                    sub_fa = fa_string[i:i + 100]
                    fil.write(sub_fa + "\n")

            fil.close()

        if stats:
            covstats = ms.covStats()
            TOTAL = [0, 0, 0, 0, 0]  # numACGT numsameasref numN numZ
            for ch, calls in fasta.items():
                statistics[ch][0] = calls[1:len(calls)].count("A") + calls.count("C") + calls.count("T") + calls.count(
                    "G")
                statistics[ch][2] = calls[1:len(calls)].count("N")
                statistics[ch][3] = calls[1:len(calls)].count("-")
                statistics[ch] += covstats[ch]
                for i in range(len(TOTAL)):
                    TOTAL[i] += statistics[ch][i]
            fil = open(stats, "w")

            TOTAL += [0, 0, 0, 0, 0, 0]

            fil.write(
                "guuid\trefid\tcontigid\tnumACGT\tnumsameasref\tnumN\tnumZ\tnumHz\tmean\tmedian\tmode\tcov1\tcov5\tcov10\n")
            for ch, info in statistics.items():
                fil.write("\t".join([guuid, refid, ch] +
                                    [str(i) for i in info]) + "\n")
            fil.write("\t".join([guuid, refid, "ALL"] +
                                [str(i) for i in TOTAL]))

            fil.close()


############ START FILTERS ###############################################
############ START FILTERS ###############################################
############ START FILTERS ###############################################
############ START FILTERS ###############################################

FILTERS = VcfFilter()

FILTERS.addFilter('-Q', 'MINSNPMAPQ', float, 'Q{0}',
                  'RMS mapping quality of variant site >= threshold [not is_indel and INFO:MQ >= FLOAT]',
                  '##FILTER=<ID=Q{0},Number=0,Type=Float,Description="Variant base not called because mapping quality less than {0}">',
                  lambda indel, record, filter_val: not indel and record[
                      GormVcf.fld_PARSED_INFO]["MQ"] < filter_val
                  )
FILTERS.addFilter('-q', 'MINGAPMAPQ', float, 'q{0}',
                  'RMS mapping quality of indel >= threshold [is_indel and INFO:MQ >= FLOAT]',
                  '##FILTER=<ID=q{0},Number=0,Type=Float,Description="Indel rejected because mapping quality less than {0}">',
                  lambda indel, record, filter_val: indel and record[
                      GormVcf.fld_PARSED_INFO]["MQ"] < filter_val
                  )
FILTERS.addFilter('-m', 'MINCNSMAPQ', float, 'm{0}',
                  'RMS mapping quality for site same as reference >= threshold [not is_indel and sum(BaseCounts)>0 and INFO:MQ >= FLOAT]',
                  '##FILTER=<ID=m{0},Number=0,Type=Float,Description="Sites same as reference not called because mapping quality less than {0}">',
                  lambda indel, record, filter_val: not indel and not (
                      sum(record[GormVcf.fld_PARSED_INFO]["BaseCounts"]) and record[GormVcf.fld_PARSED_INFO][
                          "MQ"] >= filter_val)
                  )
FILTERS.addFilter('-n', 'MINDEPTHHARD', int, 'n{0}',
                  'Number of high quality bases >= threshold [total(DP4) >= INT-equivalent]',
                  '##FILTER=<ID=n{0},Number=0,Type=Integer,Description="Base not called because number of high quality bases at this position < {0}">',
                  lambda indel, record, filter_val: sum(
                      record[GormVcf.fld_PARSED_INFO]['DP4']) < filter_val
                  )
FILTERS.addFilter('-S', 'MINSNPQ', float, 'S{0}',
                  'Phred-scaled quality score for the assertion made in ALT for a variant site >= threshold [not is_indel and QUAL>=FLOAT]',
                  '##FILTER=<ID=S{0},Number=0,Type=Float,Description="Base not called because min SNP quality less than {0}">',
                  lambda indel, record, filter_val: not indel and float(
                      record[GormVcf.fld_QUAL]) < filter_val
                  )
FILTERS.addFilter('-I', 'MININDELQ', float, 'I{0}',
                  'Phred-scaled quality score for the assertion made in ALT for start of an indel >= threshold [is_indel and QUAL>=FLOAT]',
                  '##FILTER=<ID=I{0},Number=0,Type=Float,Description="Indel not recognised because min INDEL quality less than {0}">',
                  lambda indel, record, filter_val: indel and record[GormVcf.fld_ALT] != '.' and float(
                      record[GormVcf.fld_QUAL]) < filter_val
                  )

# [
# 'VARNEARIND',
# 'Variant with <= n callable indels in a window of w bp on either side [not is_indel and >= N passed indels within +/ w bp] (note: the value must be encased on double-quotes to parse at the command-line)',
# [
# 'VARNEARVAR',
# 'Variant with <= n callable variants in a window of w bp on either side [not is_indel and >= N passed variants within +/- w bp] (note: the value must be encased on double-quotes to parse at the command-line)',

FILTERS.addFilter('-i', 'REMOVEINDELS', bool, 'i',
                  'Site is not the start of an indel [not is_indel]',
                  '##FILTER=<ID=i,Number=0,Type=None, Description="Position cannot be called because it marks the start of an indel">',
                  lambda indel, record, filter_val: indel
                  )
FILTERS.addFilter('-z', 'REMOVEHETS', bool, 'z',
                  'Site has only a homozygous genotype under a diploid model [VALUES:GT of form X/Y where X!=Y]',
                  '##FILTER=<ID=z,Number=0,Type=None,Description="Base not called because SAMtools made an initial heterozygote call (genotype X/Y where X!=Y)">',
                  lambda indel, record, filter_val: len(
                      set(record[GormVcf.fld_PARSED_DATA].setdefault("GT", "./.").split("/"))) > 1
                  )
FILTERS.addFilter('-B', 'BOTHSTRANDS', int, 'B{0}',
                  'Reads in each direction >= threshold [not is_covz and INFO:DP4fwd >= threshold and INFO:DP4rev >= threshold]',
                  '##FILTER=<ID=B{0},Number=0,Type=Integer,Description="Base not called because not supported by at least {0} high quality read(s) in forward and {0} in reverse directions">',
                  lambda indel, record, filter_val: not (
                      record[GormVcf.fld_PARSED_INFO]["DP4"][0] + record[GormVcf.fld_PARSED_INFO]["DP4"][
                          2] >= filter_val and record[GormVcf.fld_PARSED_INFO]["DP4"][1] +
                      record[GormVcf.fld_PARSED_INFO]["DP4"][3] >= filter_val)
                  )
FILTERS.addFilter('-p', 'REMOVESBRS', bool, 'p',
                  'Site is not in a repeated region of the genome as identified by a self-self BLAST of the reference [not is_indel and INFO:SBR == 0]',
                  '##FILTER=<ID=p,Number=0,Type=None,Description="Base not called because in a duplicated region of the reference (as estimated by self-self blast)">',
                  lambda indel, record, filter_val: record[GormVcf.fld_PARSED_INFO]["SBR"] == 1
                  )
FILTERS.addFilter('-N', 'MARKNONGENUINEINDELS', bool, 'N',
                  'INDEL VCF record is not a genuine INDEL',
                  '##FILTER=<ID=N,Number=0,Type=None,Description="INDEL not called because the BCALL could not be made or the BCALL is the same as the reference">',
                  lambda indel, record, filter_val: indel and (
                      record[GormVcf.fld_PARSED_INFO]['BCALL'] == UNCALLABLECHAR or record[GormVcf.fld_PARSED_INFO][
                          'BCALL'] == record[GormVcf.fld_REF])
                  )
FILTERS.addFilter('-f', 'MINPROPHQBASES', float, 'f{0}',
                  'Base not called because proportion of the high-quality bases sum(BaseCounts4) < threshold of total spanning reads DP',
                  '##FILTER=<ID=f{0},Number=0,Type=None,Description="Base not called because proportion of the high-quality mapped bases sum(BaseCounts4) < {0} of total spanning reads DP">',
                  lambda indel, record, filter_val: not indel and (not record[GormVcf.fld_PARSED_INFO]["DP"] or sum(
                      record[GormVcf.fld_PARSED_INFO]["BaseCounts4"]) / float(
                      record[GormVcf.fld_PARSED_INFO]["DP"]) < filter_val)
                  )
FILTERS.addFilter('-A', 'MINABQ4', float, 'A{0}',
                  'Mean of the high-quality base qualities >= threshold [INFO:ABQ4 >= threshold]',
                  '##FILTER=<ID=A{0},Number=0,Type=None,Description="Base not called because mean of the high-quality (DP4) base qualities < {0}">',
                  lambda indel, record, filter_val: not indel and record[
                      GormVcf.fld_PARSED_INFO]['ABQ4'] < filter_val
                  )
# FILTERS.addFilter( '-c','MAXDZ',float,'c{0}',
# 'Z-score of non-gap depth <= threshold [INFO:DZ <= threshold]',
# '##FILTER=<ID=c{0},Number=0,Type=None,Description="Base not called because positive Z-score of non-zero absolute depth (DZ) > {0}">',
# lambda indel,record,filter_val: not indel and record[GormVcf.fld_PARSED_INFO]["DZ"] > filter_val
# )
# FILTERS.addFilter( '-C','MAXDM',float,'C{0}',
# 'MAD of non-gap depth <= threshold [INFO:DM <= threshold]',
# '##FILTER=<ID=C{0},Number=0,Type=None,Description="Base not called because positive mean absolute deviation (MAD) of non-zero absolute depth (DM) > {0}">',
# lambda indel,record,filter_val: not indel and record[GormVcf.fld_PARSED_INFO]["DM"] > filter_val
# )
FILTERS.addFilter('-e', 'MAXDZ4', float, 'e{0}',
                  'Z-score of high-quality depth <= threshold [INFO:DZ4 <= threshold]',
                  '##FILTER=<ID=e{0},Number=0,Type=None,Description="Base not called because positive Z-score of high-quality depth (DZ4) > {0}">',
                  lambda indel, record, filter_val: not indel and record[
                      GormVcf.fld_PARSED_INFO]["DZ4"] > filter_val
                  )
FILTERS.addFilter('-E', 'MAXDM4', float, 'E{0}',
                  'MAD of high-quality depth <= threshold [INFO:DM4 <= threshold]',
                  '##FILTER=<ID=E{0},Number=0,Type=None,Description="Base not called because positive mean absolute deviation (MAD) of high-quality depth (DM4) > {0}">',
                  lambda indel, record, filter_val: not indel and record[
                      GormVcf.fld_PARSED_INFO]["DM4"] > filter_val
                  )
FILTERS.addFilter('-g', 'MAXDZ4L', float, 'g{0}',
                  'Z-score of GC-corrected high-quality depth <= threshold [INFO:DZ4L <= threshold]',
                  '##FILTER=<ID=g{0},Number=0,Type=None,Description="Base not called because positive Z-score of GC-corrected high-quality depth (DZ4L) > {0}">',
                  lambda indel, record, filter_val: not indel and record[
                      GormVcf.fld_PARSED_INFO]["DZ4L"] > filter_val
                  )
FILTERS.addFilter('-G', 'MAXDM4L', float, 'G{0}',
                  'MAD of GC-corrected high-quality depth <= threshold [INFO:DM4L <= threshold]',
                  '##FILTER=<ID=G{0},Number=0,Type=None,Description="Base not called because positive mean absolute deviation (MAD) of the GC-corrected high-quality depth (DM4L) > {0}">',
                  lambda indel, record, filter_val: not indel and record[
                      GormVcf.fld_PARSED_INFO]["DM4L"] > filter_val
                  )
# FILTERS.addFilter( '-k','MINPCALL',float,'k{0}',
# 'Proportion of the non-gap depth that supports the BCALL base call >= threshold [INFO:PCALL >= threshold]',
# '##FILTER=<ID=k{0},Number=0,Type=None,Description="Base not called because proportion of the non-zero absolute depth that supports the BCALL base call < {0}">',
# lambda indel,record,filter_val: not indel and record[GormVcf.fld_PARSED_INFO]['PCALL'] < filter_val
# )
FILTERS.addFilter('-K', 'MINPCALL4', float, 'K{0}',
                  'Proportion of the high-quality DP4 bases that supports the BCALL base call >= threshold [INFO:PCALL4 >= threshold]',
                  '##FILTER=<ID=K{0},Number=0,Type=None,Description="Base not called because proportion of the high-quality bases that supports the BCALL base call < {0}">',
                  lambda indel, record, filter_val: not indel and record[
                      GormVcf.fld_PARSED_INFO]['PCALL4'] < filter_val
                  )

# FILTERS.addFilter( '-j','MINPCONS',float,'j{0}',
# 'Proportion of the non-gap depth (sum of BaseCounts) that supports the GT genotype calls(s) >= threshold [INFO:PCONS >= threshold]',
# '##FILTER=<ID=j{0},Number=0,Type=None,Description="Base not called because proportion of the non-zero absolute depth that supports the GT genotype call(s) < {0}">',
# lambda indel,record,filter_val: not indel and record[GormVcf.fld_PARSED_INFO]['PCONS'] < filter_val
# )

FILTERS.addFilter('-J', 'MINPCONS4', float, 'J{0}',
                  'Proportion of the high-quality depth (DP4) that supports the GT genotype call(s) >= threshold',
                  '##FILTER=<ID=J{0},Number=0,Type=None,Description="Base not called because proportion of the high-quality bases that supports the GT genotype call(s) < {0}">',
                  lambda indel, record, filter_val: not indel and record[
                      GormVcf.fld_PARSED_INFO]['PCONS4'] < filter_val
                  )

##################### END FILTERS ##################################
##################### END FILTERS ##################################
##################### END FILTERS ##################################
##################### END FILTERS ##################################

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Filter VCF file')

    parser.add_argument("-invcf", dest="invcf", required=True)
    parser.add_argument("-outvcf", dest="outvcf", required=True)
    parser.add_argument("-outvcfIndel", dest="outvcfIndel", required=True)
    parser.add_argument("-outfasta", dest="outfasta",
                        required=False, default=False)
    parser.add_argument("-outstats", dest="outstats",
                        required=False, default=False)
    parser.add_argument("-u", dest="guuid", required=True)
    parser.add_argument("-refid", dest="ref_id", required=True)

    for i in FILTERS.availableFilters():
        if i[2] == bool:
            parser.add_argument(i[0], dest=i[1], help=i[4],
                                default=False, action='store_true')
        else:
            parser.add_argument(i[0], dest=i[1], help=i[4], default="DISABLED")
    args = parser.parse_args()

    try:
        for i in FILTERS.availableFilters():
            value = getattr(args, i[1])
            if value and value != "DISABLED":
                FILTERS.setUpFilter(i[1], value)
    except ParameterError as e:
        LE.critical("Parameter setup Error: {0}".format(e.message))
        dump_exc()

    try:
        vcfFile = GormVcf(args.invcf)
        FILTERS.filterVcf(vcfFile, outvcf_path=args.outvcf, outvcfIndel=args.outvcfIndel, outfasta=args.outfasta,
                          stats=args.outstats, guuid = args.guuid, refid=args.ref_id)
        print ("Done")
    except:
        dump_exc()
