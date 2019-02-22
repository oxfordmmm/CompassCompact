#!/usr/bin/env python

import os
import argparse
import subprocess
import shlex
import uuid
import pysam
import StringIO
from lib.g4_fixsamtoolsvcf import filter_VCF
from lib.g4_merge import *
import sys
from lib.logerror import LE, dump_exc
from lib.compassconfig import COMPASSCFG

# Variant Calling
# Wrapper for tools: bcftools, samtools, gatk


tmp = '../tmp/'




class Compass_Pileup:
    def __init__(self, gapopeningerrprob=None,
                 gapexterrprob=None,
                 coeffhomopolerrs=None,
                 mingappedreadsforindels=None,
                 minfracgappedreads=None,
                 outputDPinBCF=None,
                 outputstrandbiasPval=None,
                 maxmappingquality=None,
                 skipalnswithmapqltint=None,
                 useanomalousreadpairs=None,
                 minbaseq=None,
                 baqcomp=None,
                 callsnps=None,
                 callgenotypes=None,
                 keepallaltalleles=None,
                 calcLDforadjacentsites=None,
                 scaledsubstmutrate=None,
                 indeltosubstratio=None,
                 variantifprobltint=None,
                 typeofprior=None,
                 inbam=None,
                 inref=None,
                 pileup_out="None"):

        self.mpileup_gapopeningerrprob = "-o{0}".format(
            gapopeningerrprob) if gapopeningerrprob != "DISABLED" else ""
        self.mpileup_gapexterrprob = "-e{0}".format(
            gapexterrprob) if gapexterrprob != "DISABLED" else ""
        self.mpileup_coeffhomopolerrs = "-h{0}".format(
            coeffhomopolerrs) if coeffhomopolerrs != "DISABLED" else ""
        self.mpileup_mingappedreadsforindels = "-m{0}".format(
            mingappedreadsforindels) if mingappedreadsforindels != "DISABLED" else ""
        self.mpileup_minfracgappedreads = "-F{0}".format(
            minfracgappedreads) if minfracgappedreads != "DISABLED" else ""
        self.mpileup_outputDPinBCF = "-D" if outputDPinBCF else ""
        self.mpileup_outputstrandbiasPval = "-S" if outputstrandbiasPval else ""
        self.mpileup_maxmappingquality = "-M{0}".format(
            maxmappingquality) if maxmappingquality != "DISABLED" else ""
        self.mpileup_skipalnswithmapqltint = "-q{0}".format(
            skipalnswithmapqltint) if skipalnswithmapqltint != "DISABLED" else ""
        self.mpileup_useanomalousreadpairs = "-A" if useanomalousreadpairs else ""
        self.mpileup_minbaseq = "-Q{0}".format(
            minbaseq) if minbaseq != "DISABLED" else ""
        self.mpileup_baqcomp = "-E" if baqcomp else ""

        self.bcftools_callsnps = "-c" if callsnps else ""
        self.bcftools_callgenotypes = "-g" if callgenotypes else ""
        self.bcftools_keepallaltalleles = "-A" if keepallaltalleles else ""
        self.bcftools_calcLDforadjacentsites = "-L" if calcLDforadjacentsites else ""
        self.bcftools_scaledsubstmutrate = "-t{0}".format(
            scaledsubstmutrate) if scaledsubstmutrate != "DISABLED" else ""
        self.bcftools_indeltosubstratio = "-i{0}".format(
            indeltosubstratio) if indeltosubstratio != "DISABLED" else ""
        self.bcftools_variantifprobltint = "-p{0}".format(
            variantifprobltint) if variantifprobltint != "DISABLED" else ""
        self.bcftools_typeofprior = "-P{0}".format(typeofprior)
        
        #Check and creat folder
        if not os.path.exists(tmp):
            os.mkdir(tmp)
        # Temp folder
        bam = tmp + str(uuid.uuid4()) + ".bam"
        
        os.symlink(os.path.join(os.getcwd(), inbam), bam)
        self.opts_inbam = bam
        self.opts_inref = inref

        self.opts_out_vars = self.tempPath()
        self.opts_out_varsDP = self.tempPath()
        self.opts_out_pileup = pileup_out
        self.opts_out_pileupDP = self.tempPath()
        self.opts_out_indels = self.tempPath()
        self.opts_out_varsAnnot = self.tempPath()
        self.opts_out_varsDPAnnot = self.tempPath()
        self.stderror = open(self.tempPath(), "w+")

        self.executedCommands = []

    def tempPath(self):
        return tmp + str(uuid.uuid4())

    def substitutePars(self, cad):
        vardict = dict([(i, getattr(self, i)) for i in dir(self) if
                        i.startswith("mpileup") or i.startswith("bcftools") or i.startswith("opts")])
        LE.debug("Running command: [{0}]".format(cad.format(**vardict)))
        cmdStr = cad.format(**vardict)
        return cmdStr

    def runPileup(self):
        p1 = COMPASSCFG['tools']['samtools'].popen(append=self.substitutePars(
            "mpileup -f {opts_inref} -B -M0 -Q0 -q0 {mpileup_gapopeningerrprob} {mpileup_gapexterrprob} {mpileup_coeffhomopolerrs} {mpileup_mingappedreadsforindels} {mpileup_minfracgappedreads} {mpileup_outputDPinBCF} -g {mpileup_outputstrandbiasPval} -A {opts_inbam}"),
            stdout=subprocess.PIPE, stderr=self.stderror)
        p2 = COMPASSCFG['tools']['bcftools'].popen(append=self.substitutePars(
            "view {bcftools_callsnps} {bcftools_callgenotypes} -b {bcftools_keepallaltalleles} {bcftools_calcLDforadjacentsites} {bcftools_scaledsubstmutrate} {bcftools_indeltosubstratio} {bcftools_variantifprobltint} {bcftools_typeofprior} -"),
            stdin=p1.stdout, stdout=subprocess.PIPE, stderr=self.stderror)
        p3 = COMPASSCFG['tools']['bcftools'].popen(append=self.substitutePars("view -"), stdin=p2.stdout,
                                                   stdout=subprocess.PIPE, stderr=self.stderror)
        clean_file = filter_VCF(p3.stdout)
        clean_file.store(self.opts_out_varsDP, None, True)

        self.executedCommands.append(p1.cmd)
        self.executedCommands.append(p2.cmd)
        self.executedCommands.append(p3.cmd)

        p3, p2, p1 = p3.wait(), p2.wait(), p1.wait()
        
        if p1 or p2 or p3:
            raise Exception("Error in Pileup")

        p1 = COMPASSCFG['tools']['samtools'].popen(append=self.substitutePars(
            "mpileup -f {opts_inref} {mpileup_baqcomp} {mpileup_maxmappingquality} {mpileup_minbaseq} {mpileup_skipalnswithmapqltint} {mpileup_gapopeningerrprob} {mpileup_gapexterrprob} {mpileup_coeffhomopolerrs} {mpileup_mingappedreadsforindels} {mpileup_minfracgappedreads} {mpileup_outputDPinBCF} -g {mpileup_outputstrandbiasPval} {mpileup_useanomalousreadpairs} {opts_inbam}"),
            stdout=subprocess.PIPE, stderr=self.stderror)
        p2 = COMPASSCFG['tools']['bcftools'].popen(append=self.substitutePars(
            "view {bcftools_callsnps} {bcftools_callgenotypes} -b {bcftools_keepallaltalleles} {bcftools_calcLDforadjacentsites} {bcftools_scaledsubstmutrate} {bcftools_indeltosubstratio} {bcftools_variantifprobltint} {bcftools_typeofprior} -"),
            stdin=p1.stdout, stdout=subprocess.PIPE, stderr=self.stderror)
        p3 = COMPASSCFG['tools']['bcftools'].popen(append=self.substitutePars("view -"), stdin=p2.stdout,
                                                   stdout=subprocess.PIPE, stderr=self.stderror)
        
        # print p3
        
        clean_file = filter_VCF(p3.stdout)
        clean_file.store(self.opts_out_vars, self.opts_out_indels, False)

        self.executedCommands.append(p1.cmd)
        self.executedCommands.append(p2.cmd)
        self.executedCommands.append(p3.cmd)

        p3, p2, p1 = p3.wait(), p2.wait(), p1.wait()
        if p1 or p2 or p3:
            raise Exception("Error in Pileup")

        pysam.index(self.opts_inbam)

        out = open(self.opts_out_pileupDP, "w")
        p1 = COMPASSCFG['tools']['samtools'].popen(append=self.substitutePars(
            "mpileup -f {opts_inref} -B -M0 -Q0 -q0 {mpileup_gapopeningerrprob} {mpileup_gapexterrprob} {mpileup_coeffhomopolerrs} {mpileup_mingappedreadsforindels} {mpileup_minfracgappedreads} {mpileup_outputDPinBCF} -s -O -A {opts_inbam}"),
            stdout=out, stderr=self.stderror)
        if p1.wait():
            raise Exception("Error en pileup")
        self.executedCommands.append(p1.cmd)
        out.close()
        out = open(self.opts_out_pileup, "w")
        p1 = COMPASSCFG['tools']['samtools'].popen(append=self.substitutePars(
            "mpileup -f {opts_inref} {mpileup_baqcomp} {mpileup_maxmappingquality} {mpileup_minbaseq} {mpileup_skipalnswithmapqltint} {mpileup_gapopeningerrprob} {mpileup_gapexterrprob} {mpileup_coeffhomopolerrs} {mpileup_mingappedreadsforindels} {mpileup_minfracgappedreads} {mpileup_outputDPinBCF} -s -O {mpileup_useanomalousreadpairs} {opts_inbam}"),
            stdout=out, stderr=self.stderror)
        if p1.wait():
            raise Exception("Error en pileup")
        self.executedCommands.append(p1.cmd)
        out.close()

    def annotate(self):
        retcode = 0
        # TODO: Need to create a dict for ref using gatk since the latest veriosn requires (thanhlv)
        # picard CreateSequenceDictionary R=ref_cholera_m622.fasta O=ref_cholera_m622.dict
        # https: // software.broadinstitute.org / gatk / documentation / article?id = 1601
        
        #Need to use specific java 7 for gatk version 1.4
        _prepend = "java - jar"
        if "JAVA7" in os.environ:
            _prepend = os.path.join(os.environ["JAVA7"], "java -Xmx8192m -Xms1024m -jar")
            
        p1 = COMPASSCFG['tools']['gatk'].popen(source="jar", prepend=_prepend, append=self.substitutePars(
            "-T VariantAnnotator -I {opts_inbam} -V {opts_out_varsDP} -l INFO -R {opts_inref} -o {opts_out_varsDPAnnot} -A BaseCounts -A GCContent -S LENIENT"),
            stderr=self.stderror, stdout=self.stderror)
        retcode += p1.wait()
        self.executedCommands.append(p1.cmd)
        p1 = COMPASSCFG['tools']['gatk'].popen(source="jar", prepend=_prepend, append=self.substitutePars(
            "-T VariantAnnotator -I {opts_inbam} -V {opts_out_vars} -l INFO -R {opts_inref} -o {opts_out_varsAnnot} -A BaseCounts -A GCContent -S LENIENT"),
            stderr=self.stderror, stdout=self.stderror)
        retcode += p1.wait()
        self.executedCommands.append(p1.cmd)
        return retcode

    def dumpStdError(self):
        self.stderror.seek(0)
        LE.error(self.stderror)

    def clean(self):
        log = self.stderror.name
        self.stderror.close()

        try:
            os.unlink(log)
        except:
            pass
        try:
            os.unlink(self.opts_out_vars)
        except:
            pass
        try:
            os.unlink(self.opts_out_varsDP)
        except:
            pass
        try:
            os.unlink(self.opts_out_pileupDP)
        except:
            pass
        try:
            os.unlink(self.opts_out_indels)
        except:
            pass
        try:
            os.unlink(self.opts_out_varsAnnot)
        except:
            pass
        try:
            os.unlink(self.opts_out_varsDPAnnot)
        except:
            pass
        try:
            os.unlink(self.opts_out_vars + ".idx")
        except:
            pass
        try:
            os.unlink(self.opts_out_varsDP + ".idx")
        except:
            pass
        try:
            os.unlink(self.opts_out_varsAnnot + ".idx")
        except:
            pass
        try:
            os.unlink(self.opts_out_varsDPAnnot + ".idx")
        except:
            pass

    def merge(self, output):
        a = GormVcfMerger(self.opts_out_varsAnnot, self.opts_out_varsDPAnnot, self.opts_out_indels, self.opts_inref,
                          self.opts_out_pileupDP, self.executedCommands)
        out = open(output, "w")
        for i in a.mergeHeaders():
            out.write(i + "\n")
        for i in a.mergeRecords():
            out.write(i + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Compass pileup wrapper')

    parser.add_argument("-o", dest="mpileup_gapopeningerrprob", default="DISABLED",
                        help="Phred-scaled gap open sequencing error probability")
    parser.add_argument("-e", dest="mpileup_gapexterrprob", default="DISABLED",
                        help="Phred-scaled gap extension seq error probability")
    parser.add_argument("-H", dest="mpileup_coeffhomopolerrs", default="DISABLED",
                        help="coefficient for homopolymer errors")
    parser.add_argument("-m", dest="mpileup_mingappedreadsforindels", default="DISABLED",
                        help="minimum gapped reads for indel candidates")
    parser.add_argument("-F", dest="mpileup_minfracgappedreads", default="DISABLED",
                        help="minimum fraction of gapped reads for candidates")
    parser.add_argument("-D", dest="mpileup_outputDPinBCF", default=False, action="store_true",
                        help="output.bam per-sample DP in BCF")
    parser.add_argument("-S", dest="mpileup_outputstrandbiasPval", default=False, action="store_true",
                        help="output.bam per-sample strand bias P-value in BCF")
    parser.add_argument("-M", dest="mpileup_maxmappingquality",
                        default="DISABLED", help="cap mapping quality at INT")
    parser.add_argument("-q", dest="mpileup_skipalnswithmapqltint", default="DISABLED",
                        help="skip alignments with mapQ smaller than INT")
    parser.add_argument("-A", dest="mpileup_useanomalousreadpairs", default=False, action="store_true",
                        help="count anomalous read pairs")  # only for S.aur
    parser.add_argument("-Q", dest="mpileup_minbaseq", default="DISABLED",
                        help="skip bases with baseQ/BAQ smaller than INT")
    parser.add_argument("-E", dest="mpileup_baqcomp", default=False, action="store_true",
                        help="extended BAQ for higher sensitivity but lower specificity")

    parser.add_argument("-c", dest="bcftools_callsnps",
                        default=False, action="store_true", help="SNP calling")
    parser.add_argument("-g", dest="bcftools_callgenotypes", default=False, action="store_true",
                        help="call genotypes at variant sites")
    parser.add_argument("-K", dest="bcftools_keepallaltalleles", default=False, action="store_true",
                        help="keep all possible alternate alleles at variant sites")
    parser.add_argument("-L", dest="bcftools_calcLDforadjacentsites", default=False, action="store_true",
                        help="calculate LD for adjacent sites")
    parser.add_argument("-t", dest="bcftools_scaledsubstmutrate", default="DISABLED",
                        help="scaled substitution mutation rate")
    parser.add_argument("-i", dest="bcftools_indeltosubstratio",
                        default="DISABLED", help="indel-to-substitution ratio")
    parser.add_argument("-p", dest="bcftools_variantifprobltint",
                        default="DISABLED", help="variant if P(ref|D)<FLOAT")
    parser.add_argument("-P", dest="bcftools_typeofprior", default="full",
                        type=lambda x: parser.error("Error in -P [full,cond2,flat]") if x not in ["full", "cond2",
                                                                                                  "flat"] else x,
                        help="type of prior: full, cond2, flat")

    parser.add_argument("-B", dest="inbam", required=True,
                        type=lambda x: parser.error("File does not exist [{0}]".format(x)) if not os.path.isfile(
                            x) else x, help="Input BAM")
    parser.add_argument("-R", dest="ref_id", required=True,
                        help="Path to FASTA reference")  # type=lambda x: parser.error("File does not exist [{0}]".format(x)) if not os.path.isfile(x) else x

    parser.add_argument("-out", dest="output",
                        required=True, help="Output Vcf file")
    parser.add_argument("-outpileup", dest="outpileup",
                        required=True, help="Output mpileup file")
    parser.add_argument("-out_reference", dest="out_reference",
                        help="Reference output.bam path")

    options = parser.parse_args()

    # python g4_pileup.py -o 40 -e 20 -H 100 -m 2 -F 0.002 -D -S -M0 -q 30 -A
    # -Q25 -E -c -g -K -L -t0.01 -i -1 -p0.5 -P full -B ../test/in.bam -R
    # ../test/R00000039.fa -samtools ../mmm/src/samtools/v0.1.18/ -gatk
    # ../mmm/src/gatk/v1.4.21/ -out out.vcf -outpileup pileup.txt

    # try:
    c = Compass_Pileup(gapopeningerrprob=options.mpileup_gapopeningerrprob,
                           gapexterrprob=options.mpileup_gapexterrprob,
                           coeffhomopolerrs=options.mpileup_coeffhomopolerrs,
                           mingappedreadsforindels=options.mpileup_mingappedreadsforindels,
                           minfracgappedreads=options.mpileup_minfracgappedreads,
                           outputDPinBCF=options.mpileup_outputDPinBCF,
                           outputstrandbiasPval=options.mpileup_outputstrandbiasPval,
                           maxmappingquality=options.mpileup_maxmappingquality,
                           skipalnswithmapqltint=options.mpileup_skipalnswithmapqltint,
                           useanomalousreadpairs=options.mpileup_useanomalousreadpairs,
                           minbaseq=options.mpileup_minbaseq, baqcomp=options.mpileup_baqcomp,
                           callsnps=options.bcftools_callsnps,
                           callgenotypes=options.bcftools_callgenotypes,
                           keepallaltalleles=options.bcftools_keepallaltalleles,
                           calcLDforadjacentsites=options.bcftools_calcLDforadjacentsites,
                           scaledsubstmutrate=options.bcftools_scaledsubstmutrate,
                           indeltosubstratio=options.bcftools_indeltosubstratio,
                           variantifprobltint=options.bcftools_variantifprobltint,
                           typeofprior=options.bcftools_typeofprior,
                           inbam=options.inbam, inref=options.ref_id,
                           pileup_out=options.outpileup)
    try:
        c.runPileup()
        
        if c.annotate():
            print "Error with annotation!"
            c.dumpStdError()
            sys.exit(-1)
        c.merge(options.output)
        c.clean()
        LE.info("Finished!")
    
    except:
        c.dumpStdError()
        c.clean()
        dump_exc()
        try:
            pass
            c.dumpStdError()
            c.clean()
        except:
            pass
        dump_exc()
