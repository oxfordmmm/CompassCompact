#!/usr/bin/env python

# ============================================================================ #
# gorm_fixsamtoolsvcf.py													   #
# ============================================================================ #
"""
Make a legal VCF file

Ideally, we would like to call bases in the sample wrt the reference for
every site in the reference that has non-zero coverage in the sample by
omitting the bcftools option '-v output.bam potential variant sites only'.

However, since using this option from the bcftools program in version
0.1.12-10 (r896) of the samtools distribution only outputs "PL:DP:SP" of the
FORMAT information instead of "GT:PL:DP:SP:GQ", this script adds tags and
values for genotype GT = "0/0" and genotype quality GQ = "99" so that
subsequent programs have a consistent format to work with.

Has been updated to work on either Gorm v2 files or Gorm v3 (and upwards),
where the version is obtained from the current gorm.ini file.

Input  : VCF file path (or stdin if path specified as '-')

Output : stdout (VCF format)
"""
# ============================================================================
# Madeleine Cule									 Camilla Ip
# madeleine.cule@stats.ox.ac.uk					  camilla.ip@stats.ox.ac.uk
# December 2010									  January 2011
#
# Optimized and restructured by Carlos del Ojo (June 2013)
# carlos.delojoelias@ndm.ox.ac.uk
# Juny 2013
# ============================================================================

import os
import sys
import itertools
import argparse
import re
from logerror import LE, dump_exc


class filter_VCF():
    """Filter the VCF file produced by Gorm v3+."""

    def __init__(self, fp):
        self.lenCR = None
        self.fp = fp

    def getCleanFile(self):
        INFOTAGS = set()
        lenCR = None
        buff = []
        LE.info("Filtering headers, adding additional headers.")
        for l in self.fp:
            if not l.startswith("#"):
                buff.append(l)
                break

            match = re.findall("##INFO=<ID=([^,]+),", l)
            if match:
                INFOTAGS.add(match[0])

            if not lenCR:
                lenCR = len(l) - len(l.strip())

            l = l[:-lenCR]
            if l.startswith('##FORMAT=<ID=PL,Number=-1,'):
                # Don't do anything - there's a good line and an invalid line -
                # just remove the invalid one
                pass
            elif l.startswith('##INFO=<ID=LowMQ,Number=3,Type=Integer'):
                l = l.replace('##INFO=<ID=LowMQ,Number=3,Type=Integer',
                              '##INFO=<ID=LowMQ,Number=3,Type=Float')
                yield l, None
            elif l.startswith('##INFO=<ID=GC,Number=1,Type=Integer,'):
                l = l.replace('##INFO=<ID=GC,Number=1,Type=Integer,',
                              '##INFO=<ID=GC,Number=1,Type=Float,')
                yield l, None
            elif l.startswith('##INFO=<ID=MQ,Number=1,Type=Float,Description="RMS Mapping Quality">'):
                # There's two of these MQ lines - remove the less informative
                # one.
                pass
            elif l.startswith('##INFO=<ID=MQ,Number=1,Type=Integer,'):
                l = l.replace('##INFO=<ID=MQ,Number=1,Type=Integer,',
                              '##INFO=<ID=MQ,Number=1,Type=Float,')
                yield l, None
            elif l.startswith('##INFO=<ID=LowMQ,Number=3,Type=Integer,'):
                l = l.replace('##INFO=<ID=LowMQ,Number=3,Type=Integer,',
                              '##INFO=<ID=LowMQ,Number=3,Type=Float,')
                yield l, None
            else:
                yield l, None

        LE.info("Filtering VCF records...")
        for l in itertools.chain(buff, self.fp):
            line_split = l[:-lenCR].split("\t")
            indel = "INDEL" in line_split[7]

            fields = line_split[7].split(";")
            line_split[7] = ";".join(
                i for i in fields if i.split("=")[0] in INFOTAGS)

            if not indel:
                # CI 23-Apr-2012: Some of the VCF files produced either by samtools v0.1.18, bcftools v0.1.18
                # or GATK VariantAnnotator v1.5-31 produce INFO fields that have either a ';;' or start or end with a
                # trailing ';' (i.e., starts with a blank INFO field). These
                # need to be removed.

                if line_split[4] == 'X':
                    line_split[4] = "."

                if ",X" in line_split[4]:
                    line_split[4] = line_split[4].replace(",X", "")

                if line_split[7][0] == ';' or ';;' in line_split[7] or line_split[7][-1] == ';':
                    line_split[7] = ';'.join(
                        [x for x in line_split[7].split(';') if x])

                # CI 24-Jan-2011: This was written by Madeleine when she was
                # testing.
                if (line_split[8] == "PL"):
                    line_split[8] = "GT:" + line_split[8]
                    line_split[9:] = ["0/0:" + x for x in line_split[9:]]
                # CI 24-Jan-2011: This was written by Madeleine when she was
                # testing.
                elif (line_split[8] == "GT:GQ:PL"):
                    pl, gt, gq = line_split[9].split(":")
                    line_split[9] = ":".join(
                        line_split[9].split(":")[:3].rotate(-1))
                # CI 24-Jan-2011: This is the section that is used by the new VCF-based filtering method.
                # For every 'same as reference' position called by samtools mpileup, the 'FORMAT' fields are 'PL:DP:SP' instead of
                # 'GT:PL:DP:SP:GQ' so this section just replaces the first FORMAT field with 'GT:PL:DP:SP:GQ' and every subsequent
                # format value with '0/0:plvalue:dpvalue:spvalue:99'. The first transformation line operates on everything after
                # the FORMAT specifier, then the second line operates on every set of values for all readgroups on the rest of the
                # line.
                elif (line_split[8] == "PL:DP:SP"):
                    line_split[8] = "GT:" + line_split[8] + ":GQ"
                    for i in range(9, len(line_split)):
                        line_split[i] = "0/0:" + line_split[i] + ":99"

            yield "%s" % "\t".join(line_split), indel

    def store(self, out, out_indel, noindels):
        out = open(out, "w")

        if out_indel:
            out_indel = open(out_indel, "w")

        gen = self.getCleanFile()

        if noindels:
            for line, typ in gen:
                if typ != True:
                    out.write(line + "\n")

        elif not out_indel:
            for line, j in gen:
                out.write(line + "\n")
        else:
            for line, typ in gen:
                if typ != None:
                    break
                out.write(line + "\n")
                out_indel.write(line + "\n")

            for line, typ in itertools.chain([(line, typ)], gen):
                if typ:
                    out_indel.write(line + "\n")
                else:
                    out.write(line + "\n")


if __name__ == "__main__":
    '''bcftools view C00001770_R00000003.noa.Q.bcf | gorm_fixsamtoolsvcf.py | bcftools view -b -S - > C00001770_R00000003.noa.Q.fixd.bcf'''
    parser = argparse.ArgumentParser(
        description='Merge VCF files (HighQ depth + some depth + Indels)')
    parser.add_argument("-noindels", dest="noindels",
                        action="store_true", default=False)
    parser.add_argument("-o", dest="output.bam", required=True,
                        help="Output Vcf file")
    parser.add_argument("-oi", dest="indelout", default=None,
                        help="Output Vcf file for indels [Disabled]")
    args = parser.parse_args()

    if args.noindels and args.indelout:
        print "Error: -noindels and -oi are mutual exclusive"
        sys.exit(-1)

    clean_file = filter_VCF(sys.stdin)
    clean_file.store(args.output, args.indelout, args.noindels)
