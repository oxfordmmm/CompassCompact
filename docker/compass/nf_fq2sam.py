#!/usr/bin/env python
from __future__ import print_function
from __future__ import print_function
import argparse
import re
from lib.logerror import LE, dump_exc
import subprocess
import sys
import StringIO
from lib.fqreader import *
import itertools
from threading import Thread
from lib.compassconfig import COMPASSCFG


class Fq2Sam:
    EMPTYHEADER = '''@HD	VN:1.0	SO:coordinate
@SQ	SN:unmapped	LN:0
@RG	ID:{rg}	PL:{platform}	PU:unknown	LB:{lib}	SM:{samplename}	CN:{centre}	PI:0
'''

    def __init__(self, header, readgroup, fq1, fq2, out, shuffled=False):
        self.header = StringIO.StringIO(header)
        self.rg = readgroup

        fqfiles = [fq1, fq2]

        self.fqfiles = []
        for i in fqfiles:
            if not i:
                continue
            if type(i) == FastQReader:
                self.fqfiles.append(i)
            else:
                self.fqfiles.append(FastQReader(i))
        self.singleshuffled = shuffled
        self.out = out
        self.thread = None

    def run(self):
        self.thread = Thread(target=self.worker, kwargs={})
        self.thread.start()

    def worker(self):
        for i in self:
            self.out.write(i)

    def join(self):
        self.thread.join()

    def SingleEndDump(self, f1):
        for i in f1:
            yield i[0] + "\t4\t*\t0\t0\t*\t*\t0\t0\t" + i[2] + "\t" + i[4] + "\tRG:Z:" + self.rg + "\n"

    def ShuffledPairedEndDump(self, f1):
        reads = {}
        for i in f1:
            name = i[0]
            l = reads.setdefault(name, [])
            l.append(i)
            if len(l) == 2:
                del reads[name]
                i = l[0]
                yield name + "\t77\t*\t0\t0\t*\t*\t0\t0\t" + i[2] + "\t" + i[4] + "\tRG:Z:" + self.rg + "\n"
                i = l[1]
                yield name + "\t141\t*\t0\t0\t*\t*\t0\t0\t" + i[2] + "\t" + i[4] + "\tRG:Z:" + self.rg + "\n"

        for name, i in reads.items():
            if len(i) == 1:
                i = i[0]
                yield name + "\t77\t*\t0\t0\t*\t*\t0\t0\t" + i[2] + "\t" + i[4] + "\tRG:Z:" + self.rg + "\n"
                yield name + "\t141\t*\t0\t0\t*\t*\t0\t0\tN\tA\tRG:Z:" + self.rg + "\n"

    def PairedEndDump(self, f1, f2):
        reads = {}
        for i, j in itertools.izip_longest(f1, f2, fillvalue=("", "", "", "", "")):
            name1 = i[0]
            name2 = j[0]
            l1 = reads.setdefault(name1, [None, None])
            l2 = reads.setdefault(name2, [None, None])
            l1[0] = i
            l2[1] = j
            if l1[0] and l1[1]:
                k, l = l1
                yield k[0] + "\t77\t*\t0\t0\t*\t*\t0\t0\t" + k[2] + "\t" + k[4] + "\tRG:Z:" + self.rg + "\n"
                yield l[0] + "\t141\t*\t0\t0\t*\t*\t0\t0\t" + l[2] + "\t" + l[4] + "\tRG:Z:" + self.rg + "\n"
                l1[0] = l1[1] = None
                del reads[k[0]]

            if l2[0] and l2[1]:
                k, l = l2
                yield k[0] + "\t77\t*\t0\t0\t*\t*\t0\t0\t" + k[2] + "\t" + k[4] + "\tRG:Z:" + self.rg + "\n"
                yield l[0] + "\t141\t*\t0\t0\t*\t*\t0\t0\t" + l[2] + "\t" + l[4] + "\tRG:Z:" + self.rg + "\n"
                l2[0] = l2[1] = None
                del reads[name1]

    def __iter__(self):
        self.header.seek(0)

        for i in self.header:
            yield i
        if self.singleshuffled:
            for i in self.ShuffledPairedEndDump(self.fqfiles[0]):
                yield i
        elif len(self.fqfiles) == 1:
            for i in self.SingleEndDump(self.fqfiles[0]):
                yield i
        else:
            for i in self.PairedEndDump(self.fqfiles[0], self.fqfiles[1]):
                yield i


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Fix FastqFiles')
    parser.add_argument('-s', action='store_true', dest="shuffled",
                        help="Fastq file with shuffled fastQ reads (paired end data) (- for stdin)", default=False)
    parser.add_argument('-fq1', dest="fq1", help="First Fastq file (can be single end data) (- for stdin)",
                        default=None)
    parser.add_argument('-fq2', dest="fq2",
                        help="Second Fastq file (second set or reads, fq1 must be specified as well)", default=None)
    parser.add_argument('-r', dest="rgroup", help="Read group", required=True)
    parser.add_argument('-H', dest="header", help="You must specifiy a file with the SAM header in text format",
                        default=None)
    parser.add_argument('-dh', dest="headerinfo",
                        help="Default ummaped header, you must specify [readgroup,platform,lib,sample,SeqCentre] ex: -dh RG0045,ILLUMINA,LIB03,SN123,Sanger",
                        default=None)
    parser.add_argument('-o', dest="output.bam", help="You must specifiy a file with the SAM header in text format",
                        default="-")
    args = parser.parse_args()

    # print " ".join(sys.argv)

    LE.info("Input fastq files in {0} , {1}".format(args.fq1, args.fq2))
    LE.info("Output bam in {0}".format(args.output))

    if args.headerinfo == "None":
        args.headerinfo = None

    if (not args.header and not args.headerinfo) or (args.header and args.headerinfo):
        print("You must specify either a header file for the SAM header or default header information (-H/-dh)")
        sys.exit(-1)

    '''if args.output.bam=="-":
        output.bam=sys.stdout
    else:
        output.bam=open(args.output.bam,"w")'''

    if not args.fq1 or (args.fq2 and args.shuffled):
        print("You must specify at least one FastQ File(fq1). If you specify two -s parameter is not compatible")
        sys.exit(0)

    if args.header:
        LE.info("Reading header from {0}".format(args.header))
        if args.rgroup == "UNKNOWN":
            head = open(args.header).read()
            if '@RG' in head:
                args.rgroup = re.findall(
                    "[0-9a-z]{8}-[0-9a-z]{4}-[0-9a-z]{4}-[0-9a-z]{4}-[0-9a-z]{12}", head, re.I)[0]
            else:
                LE.error("Can not find read group on header file")
        md = COMPASSCFG['tools']['samtools'].popen(append="view -h -S -b -o {0} -".format(args.output),
                                                   stderr=subprocess.PIPE, stdin=subprocess.PIPE)
        SamFile = Fq2Sam(open(args.header).read(), args.rgroup,
                         args.fq1, args.fq2, md.stdin, shuffled=args.shuffled)
    else:
        try:
            rg, pl, lib, sn, cn = args.headerinfo.split(",")
        except:
            LE.error(
                "Wrong format for header information! ex: -dh RG0045,ILLUMINA,LIB03,SN123,Sanger")

        header = Fq2Sam.EMPTYHEADER.format(
            rg=rg, platform=pl, lib=lib, samplename=sn, centre=cn)
        md = COMPASSCFG['tools']['samtools'].popen(append="view -h -S -b -o {0} -".format(args.output),
                                                   stderr=subprocess.PIPE, stdin=subprocess.PIPE)
        SamFile = Fq2Sam(header, args.rgroup, args.fq1,
                         args.fq2, md.stdin, shuffled=args.shuffled)

    SamFile.run()
    SamFile.join()
    md.stdin.close()
    md.stderr.read()
