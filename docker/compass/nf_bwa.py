#!/usr/bin/env python
import pysam
import sys
import tempfile
import argparse
import os
import subprocess
import shlex
from lib.fixfastq import PairFqNormalizer
from lib.fqreader import *
from lib.fqsort import fqSort
import re
import uuid
import shutil
import ntpath
from lib.logerror import LE, dump_exc
import numpy as np
from lib.compassconfig import COMPASSCFG


def compute_mad(data, med, axis=None):
    return np.median(np.abs(data - med))


class BwaMap:
    def __init__( self, input, output, ref, execute, seqstats, flagstats, guuid="-"):
        # bam="/tmp/"+str(uuid.uuid4())+".bam"
        # os.symlink(input,bam)
        self.input = input
        self.output = output
        self.ref = ref
        self.origref = ref
        self.execute = execute
        self.sampleguuid = guuid
        self.seqstatsfile = seqstats
        self.flagstatsfile = flagstats
        self.commandsHistory = []
        bamf = pysam.Samfile(input, check_sq=True, mode="rb") # Open bam file using Samfile method, which is only available for pysam < 0.8.0
        guuid = bamf.header["RG"][0]["ID"] #Extract guuid
        bamf.close() # Then close it
        self.tmpdir = tempfile.mkdtemp("tmp")
        self.path = ntpath.split(os.path.realpath(input))[0] #Path of the input
        self.insertsizes = {}
        
        cmd, stdout, stderr, errcode = COMPASSCFG["tools"]["bwa"].execute(
            append="index {0}".format(self.ref)) # Build index for ref when initialise a BwaMap instance.
        if errcode:
            raise Exception("Error indexing reference")

    def prepareCommand(self, cad):
        self.commandsHistory.append(cad) # Append command to command history
        return shlex.split(cad) #Also split the command into human readable

    def run(self):
        if self.execute:
            self.readgroups = self.get_readgroups()
            self.createFQs()
            self.map()
            self.merge()
            self.markDuplicates()
            self.generateStats()

        else:
            # no bwa, return same input
            shutil.copyfile(self.input, self.output)

    def createFQs(self):
        fileList = {}
        for i in self.readgroups:
            fileList[i] = [open(self.tmpdir + "/" + i + "-1.fq", "w"),
                           open(self.tmpdir + "/" + i + "-2.fq", "w")]

        for i in pysam.Samfile(self.input):
            rg = dict(i.tags)["RG"]
            if i.flag & 64:
                fileList[rg][0].write(
                    "@{0}\n{1}\n+\n{2}\n".format(i.qname, i.seq, i.qual))
            else:
                fileList[rg][1].write(
                    "@{0}\n{1}\n+\n{2}\n".format(i.qname, i.seq, i.qual))

        names = os.listdir(self.path)
        k = 2
        for name in names:
            pattern = re.compile(
                'output.bam[0-9]_[0-9a-z]{8}-[0-9a-z]{4}-[0-9a-z]{4}-[0-9a-z]{4}-[0-9a-z]{12}_bam.bam')
            match = re.search(pattern, name)
            if match:
                bam = os.path.join(self.path, name)
                for i in pysam.Samfile(bam):
                    rg = dict(i.tags)["RG"]
                    if i.flag & 64:
                        fileList[rg][0].write(
                            "@{0}\n{1}\n+\n{2}\n".format(i.qname, i.seq, i.qual))
                    else:
                        fileList[rg][1].write(
                            "@{0}\n{1}\n+\n{2}\n".format(i.qname, i.seq, i.qual))
                k += 1

        for i in fileList.values():
            i[0].close()
            i[1].close()

        for i in self.readgroups:
            fqSort(FastQReader(self.tmpdir + "/" + i + "-1.fq"),
                   self.tmpdir + "/" + i + "-1.sort")
            fqSort(FastQReader(self.tmpdir + "/" + i + "-2.fq"),
                   self.tmpdir + "/" + i + "-2.sort")
            os.unlink(self.tmpdir + "/" + i + "-1.fq")
            os.unlink(self.tmpdir + "/" + i + "-2.fq")

        for i in self.readgroups:
            pn = PairFqNormalizer(self.tmpdir + "/" + i + "-1.sort", self.tmpdir + "/" + i + "-2.sort",
                                  FastQWriter(self.tmpdir + "/" + i + "-A.fq"),
                                  FastQWriter(self.tmpdir + "/" + i + "-B.fq"), True, 1)
            pn.normalize()
            os.unlink(self.tmpdir + "/" + i + "-1.sort")
            os.unlink(self.tmpdir + "/" + i + "-2.sort")

        self.fqfiles = [(self.tmpdir + "/" + i + "-A.fq",
                         self.tmpdir + "/" + i + "-B.fq") for i in self.readgroups]
        LE.debug("FQs created")

    def clean(self):
        try:
            shutil.rmtree(self.tmpdir, ignore_error=True)
        except:
            pass

    def map(self):
        for i in self.readgroups:
            fq1 = self.tmpdir + "/" + i + "-A.fq"
            fq2 = self.tmpdir + "/" + i + "-B.fq"

            output = open(self.tmpdir + "/" + i + ".sam", "w")
            p = COMPASSCFG["tools"]["bwa"].popen(
                append="mem -R '@RG\\tID:{0}' {1} {2} {3} -L 20 -B 3 -O 6 -T 20".format(
                    i, self.ref, fq1, fq2),
                stderr=subprocess.PIPE, stdout=output)
            self.commandsHistory.append(p.cmd)

            LE.debug(p.stderr, "stderr")
            errcode = p.wait()

            if errcode:
                LE.error("BWA tool failed")
                raise Exception("BWA tool failed")
            os.unlink(fq1)
            os.unlink(fq2)
            output = self.tmpdir + "/" + i + ".sam"
            self.insertsizes[i] = self.generateSeqStats(output)

    def merge(self):
        LE.info("Merging SAMfiles from different readgroup mappings")
        input = pysam.Samfile(self.input)
        newheaders = dict(input.header.items())
        samheader = pysam.Samfile(
            self.tmpdir + "/" + self.readgroups[0] + ".sam")
        newheaders["SQ"] = samheader.header["SQ"]
        samheader.close()
        input.close()

        newheaders["PG"] = [
            {"PN": "bwa", "VN": "0.7.10", "CL": self.commandsHistory[0]}]
        if "CO" in newheaders.keys():
            newheaders["CO"] = list(set(newheaders["CO"]))
        else:
            newheaders["CO"] = []
        newheaders["CO"].append("CMD:{0}".format(" ".join(sys.argv)))
        for i in self.commandsHistory:
            newheaders["CO"].append("CMD:{0}".format(i))

        LE.debug("Doing merge, writing in " + self.output)
        unsortedBamName = self.tmpdir + "/" + self.output + "_unsorted.bam"
        output = pysam.Samfile(unsortedBamName, "wb", header=newheaders)

        for i in self.readgroups:
            with pysam.Samfile(self.tmpdir + "/" + i + ".sam") as source:
                for j in source:
                    if not j.flag & 2048:
                        output.write(j)
            os.unlink(self.tmpdir + "/" + i + ".sam")
        output.close()
        outputNames = self.output.split(".")
        outputPrefix = ".".join(outputNames[0: len(outputNames) - 1])
        pysam.sort(unsortedBamName, self.output)
        os.unlink(unsortedBamName)
        shutil.move(self.output + ".bam", self.output)

    def get_readgroups(self):
        groups = []
        input = pysam.Samfile(self.input)
        groups.append(input.header["RG"][0]["ID"])
        names = os.listdir(self.path)
        for i in input.header["RG"]:
            groups.append(i["ID"])
        for name in names:
            pattern = re.compile(
                'output.bam[0-9]_[0-9a-z]{8}-[0-9a-z]{4}-[0-9a-z]{4}-[0-9a-z]{4}-[0-9a-z]{12}_bam.bam')
            match = re.search(pattern, name)
            if match:
                bam = os.path.join(self.path, name)
                input = pysam.Samfile(bam)
                groups.append(input.header["RG"][0]["ID"])
        return list(set(groups))

    def markDuplicates(self):
        cmd, stdout, stderr, errcode = COMPASSCFG["tools"]["picard"].execute(source="path", prepend="java -jar",
                                                                             file="MarkDuplicates.jar",
                                                                             append="I={0} O={1}.dedup METRICS_FILE=a.txt ASSUME_SORTED=true VERBOSITY=DEBUG VALIDATION_STRINGENCY=SILENT".format(
                                                                                 self.output, self.output))
        LE.debug(StringIO(stdout), "stdout")
        LE.debug(StringIO(stderr), "stderr")

        if errcode:
            LE.critical(
                "MarkDuplicates execution failed {0}".format(md.returncode))
            raise Exception(
                "MarkDuplicates execution failed {0}".format(md.returncode))

        shutil.move(self.output + ".dedup", self.output)

    def generateSeqStats(self, inpbam):
        a = pysam.Samfile(inpbam)

        isize = {}
        sizes = []

        for i in a:
            if i.mate_is_unmapped:
                continue
            data = isize.setdefault(i.qname, [None, None])
            if i.flag & 64:
                data[0] = i
            if i.flag & 128:
                data[1] = i
            if data[0] and data[1]:

                if data[0].is_unmapped or data[1].is_unmapped:
                    continue
                if data[0].pos > data[1].pos:
                    sizes.append(data[0].pos + data[0].alen - data[1].pos)
                else:
                    sizes.append(data[1].pos + data[1].alen - data[0].pos)
                del isize[i.qname]

        try:
            npSizes = np.array(sizes)
            median = np.median(npSizes)
            mad = compute_mad(npSizes, median)
            newsizes = [i for i in sizes if i < median + 10 * mad]
            std = np.std(np.array(newsizes))

        except:
            return ["0", "0"]

        return [str(int(median)), str(int(std))]

    def getFlagstats(self, bam):
        a = pysam.Samfile(bam, "rb")

        tot = qcfail = dup = map = pair = r1 = r2 = proppair = bothmapped = singleton = mapdiffchrs = mapdiffchrsmqg5 = 0
        mates = {}

        for i in a:
            tot += 1
            if i.flag & 512:
                qcfail += 1
            if i.flag & 1024:
                dup += 1
            if i.flag & 1:
                pair += 1
            if i.flag & 2:
                proppair += 1
            if not i.flag & 4:
                map += 1
            if i.flag & 64:
                r1 += 1
            if i.flag & 128:
                r2 += 1
            if not i.flag & 12:
                bothmapped += 1
            if i.flag & 12 == 8:
                singleton += 1
            if i.flag & 2048:
                continue
            rds = mates.setdefault(i.qname, [])
            rds.append(i)
            if len(rds) == 2:
                rr1, rr2 = rds
                if not rr1.is_unmapped and not rr2.is_unmapped and rr1.tid != rr2.tid:
                    mapdiffchrs += 1
                    if rr1.mapq > 5 and rr2.mapq > 5:
                        mapdiffchrsmqg5 += 1
                del mates[i.qname]

        return tot, qcfail, dup, map, "{0:3.2f}".format(
            float(map * 100) / tot), pair, r1, r2, proppair, "{0:3.2f}".format(
            float(proppair * 100) / tot), bothmapped, singleton, "{0:3.2f}".format(
            float(singleton * 100) / tot), mapdiffchrs, mapdiffchrsmqg5

    def generateStats(self):
        flagstatValues = [self.sampleguuid, self.origref,
                          "compass"] + list(self.getFlagstats(self.output))

        outfile = open(self.flagstatsfile, "w")

        outfile.write("\t".join(
            ["comid", "refid", "datasrc", "totalreadsnum", "qcfailurenum", "duplicatesnum", "mappednum", "mappedpct",
             "pairedinseqnum", "read1num", "read2num", "properlypairednum", "properlypairedpct", "itselfandmatenum",
             "singletonsnum", "singletonspct", "todiffchrnum", "todiffchrmapQgeq5num"]) + "\n")
        outfile.write("\t".join([str(i) for i in flagstatValues]))
        outfile.close()

        outfile = open(self.seqstatsfile, "w")
        outfile.write(
            "\t".join(["seqid", "medianinsertsize", "sdinsertsize"]) + "\n")
        for rg, stats in self.insertsizes.items():
            outfile.write("\t".join([rg] + stats) + "\n")

        outfile.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Bwa wrapper')

    parser.add_argument("-b", dest="bam", required=True, help="input bam file")
    parser.add_argument("-r", dest="ref_id", required=True,
                        help="reference path", default="unknown")
    parser.add_argument("-o", dest="output.bam", help="Output file", required=True)
    parser.add_argument("-execute", dest="execute",
                        help="Run Bwa", action="store_true", default=False)
    parser.add_argument("-ss", dest="seqstats",
                        help="Seqstats Output file", required=True)
    parser.add_argument("-fs", dest="flagstats",
                        help="Flagstats Output file", required=True)
    parser.add_argument("-g", dest="guuid", help="Sample guuid", default="-")
    args = parser.parse_args()
    mapping = BwaMap(args.bam, args.output, args.ref_id,
                     args.execute, args.seqstats, args.flagstats, args.guuid)
    try:
        mapping.run()
        mapping.clean()
    except:
        mapping.clean()
        dump_exc()
