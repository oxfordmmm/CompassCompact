#!/usr/bin/env python
__author__ = 'laura'
import pysam
import argparse
import os
import ntpath
import subprocess
import shlex
import itertools
import uuid
import re
import shutil
from lib.logerror import LE, dump_exc
from collections import OrderedDict
import math
import array
import sys
from lib.compassconfig import COMPASSCFG
from StringIO import StringIO


class StampyMap:
	def __init__( self, input, output, ref, subrate, keepgoodreads, seqstats, flagstats, alignquals=False, baq=False,
	              se=False, guuid="-" ):
		self.input = input
		self.output = output
		self.origref = ref
		self.ref = ref
		self.subrate = subrate
		self.keepgoodreads = keepgoodreads
		self.alignquals = alignquals
		self.baq = baq
		self.singleend = se
		self.insertsizes = {}
		self.seqstatsfile = seqstats
		self.flagstatsfile = flagstats
		self.sampleguuid = guuid
		bamf = pysam.Samfile(input)
		guuid = bamf.header["RG"][0]["ID"]
		bamf.close()
		
		self.newsams = []
		self.commandsHistory = []
		
		self.path = ntpath.split(os.path.realpath(input))[0]

	
	@staticmethod
	def Median( x ):
		return x[len(x) / 2 + 1] if len(x) % 2 else (x[len(x) / 2 + 1] + x[len(x) / 2]) / 2
	
	@staticmethod
	def stdev( a ):
		mean = sum(a) / len(a)
		res = float(0)
		for i in a:
			res += (i - mean) ** 2
		
		return math.sqrt(res / len(a))
	
	@staticmethod
	def Compute_MAD( a, c=1, axis=None ):
		"""
		Median Absolute Deviation along given axis of an array:

		median(abs(a - median(a)))

		Copied from http://code.google.com/p/agpy/source/browse/trunk/agpy/mad.py
		Downloaded 7-Dec-2012.
		"""
		
		d = StampyMap.Median(a)
		nary = array.array('d')
		for i in a:
			nary.append(float(abs(i - d)) / c)
		nary = list(nary)
		nary.sort()
		mad = StampyMap.Median(nary)
		del nary
		
		return mad
	
	def run( self ):
		self.map()
		self.merge()
		self.markDuplicates()
		self.generateStats()
	
	def merge( self ):
		LE.debug("Doing merge, writing in " + self.output)
		filestomerge = [pysam.Samfile(i) for i in self.newsams]
		
		newHeader = {}
		
		rgs = {}
		for j in itertools.chain(*[i.header["RG"] for i in filestomerge]):
			rgs.setdefault(j["ID"], j)
		
		sqs = OrderedDict()
		for j in itertools.chain(*[i.header["SQ"] for i in filestomerge]):
			sqs.setdefault(j["SN"], j)
		
		newHeader["HD"] = filestomerge[0].header["HD"]
		newHeader["RG"] = rgs.values()
		newHeader["SQ"] = sqs.values()
		newHeader["CO"] = list(itertools.chain(*[i.header["CO"] for i in filestomerge]))
		
		pgs = list(itertools.chain(*[i.header["PG"] for i in filestomerge]))
		for i in pgs:
			newHeader["CO"].append("\t".join([":".join(k) for k in i.items()]))
		newHeader["CO"] = list(set(newHeader["CO"]))
		newHeader["CO"].append("CMD:{0}".format(" ".join(sys.argv)))
		for i in self.commandsHistory:
			newHeader["CO"].append("CMD:{0}".format(i))
		
		outBam = pysam.Samfile(self.output, "wb", header=newHeader)
		for j in filestomerge:
			for i in j:
				outBam.write(i)
		outBam.close()
		pysam.sort(self.output, self.output)
		shutil.move(self.output + ".bam", self.output)
	
	def markDuplicates( self ):
		cmd, stdout, stderr, errcode = COMPASSCFG['tools']['picard'].execute(source='path', file='MarkDuplicates.jar',
		                                                                     prepend='java -jar',
		                                                                     append="I={0} O={0}.dedup METRICS_FILE={1}_metrics.txt ASSUME_SORTED=true VERBOSITY=DEBUG VALIDATION_STRINGENCY=SILENT".format(
			                                                                     self.output, self.input))
		LE.debug(StringIO(stdout), "stdout")
		LE.debug(StringIO(stderr), "stderr")
		
		if errcode:
			LE.critical("MarkDuplicates execution failed {0}".format(errcode))
			raise Exception("MarkDuplicates execution failed {0}".format(errcode))
		
		shutil.move(self.output + ".dedup", self.output)
	
	def clean( self ):
		for i in self.newsams:
			try:
				os.unlink(i)
			except:
				pass
		
		try:
			os.unlink(self.output + ".dedup")
		except:
			pass
	
	def prepareCommand( self, cad ):
		self.commandsHistory.append(cad)
		return shlex.split(cad)
	
	def map( self ):
		sf = pysam.Samfile(self.input)
		rgs = [i['ID'] for i in sf.header["RG"]]
		sf.close()
		self.newsams = [os.path.join(os.getcwd(), str(uuid.uuid4())) + ".sam" for i in range(len(rgs))]
		
		for readgroup, samout in zip(rgs, self.newsams):
			append = "--substitutionrate={0} -g {1} -h {1} -M {2} -o {3} --logfile={3}.log --readgroup=ID:{4} --outputformat=sam -v 3 ".format(
				self.subrate, self.ref, self.input, samout, readgroup)
			
			if self.keepgoodreads: stampyCmd += " --bamkeepgoodreads "
			if self.alignquals: stampyCmd += " --alignquals "
			if self.baq: stampyCmd += " --baq "
			
			cmd, stdout, stderr, errcode = COMPASSCFG["tools"]["stampy"].execute(append=append)
			cmd, stdout, stderr, errcode = "", "", "", 0
			LE.debug(StringIO(stdout), "stdout")
			LE.debug(StringIO(stderr), "stderr")
			
			if errcode:
				LE.critical("Stampy execution failed {0}".format(errcode))
				raise Exception("Stampy execution failed {0}".format(errcode))
			
			self.insertsizes[readgroup] = self.generateSeqStats(samout)
	
	def generateSeqStats( self, inpbam ):
		a = pysam.Samfile(inpbam)
		
		isize = {}
		sizes = []
		
		for i in a:
			if i.mate_is_unmapped: continue
			
			data = isize.setdefault(i.qname, [None, None])
			if i.flag & 64:
				data[0] = i
			if i.flag & 128:
				data[1] = i
			if data[0] and data[1]:
				if data[0].is_unmapped or data[1].is_unmapped: continue
				if data[0].pos > data[1].pos:
					sizes.append(data[0].pos + data[0].alen - data[1].pos)
				else:
					sizes.append(data[1].pos + data[1].alen - data[0].pos)
				del isize[i.qname]
		try:
			sizes.sort()
			median = StampyMap.Median(sizes)
			mad = StampyMap.Compute_MAD(sizes)
			newsizes = [i for i in sizes if i < median + 10 * mad]
			std = StampyMap.stdev(newsizes)
		except:
			return ["0", "0"]
		return [str(int(median)), str(int(std))]
	
	def getFlagstats( self, bam ):
		a = pysam.Samfile(bam, "rb")
		
		tot = qcfail = dup = map = pair = r1 = r2 = proppair = bothmapped = singleton = mapdiffchrs = mapdiffchrsmqg5 = 0
		mates = {}
		
		for i in a:
			tot += 1
			if i.flag & 512: qcfail += 1
			if i.flag & 1024: dup += 1
			if i.flag & 1: pair += 1
			if i.flag & 2: proppair += 1
			if not i.flag & 4: map += 1
			if i.flag & 64: r1 += 1
			if i.flag & 128: r2 += 1
			if not i.flag & 12: bothmapped += 1
			if i.flag & 12 == 8: singleton += 1
			if i.flag & 2048: continue
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
	
	def generateStats( self ):
		flagstatValues = [self.sampleguuid, self.origref, "compass"] + list(self.getFlagstats(self.output))
		
		outfile = open(self.flagstatsfile, "w")
		
		outfile.write("\t".join(
			["sampleguuid", "refid", "comment", "total_reads_num", "failure_num", "duplicates_num", "mapped_num",
			 "mapped_pct", "paried_in_seq_num", "read_1_num", "read_2_num", "properly_paried_num",
			 "properly_paried_pct", "it_self_and_mate_num", "singletons_num", "singletons_pct", "to_diff_chr_num",
			 "to_diff_chr_map_qgeq_5_num"]) + "\n")
		outfile.write("\t".join([str(i) for i in flagstatValues]))
		outfile.close()
		
		outfile = open(self.seqstatsfile, "w")
		outfile.write("\t".join(["seqid", "medianinsertsize", "sdinsertsize"]) + "\n")
		for rg, stats in self.insertsizes.items():
			outfile.write("\t".join([rg] + stats) + "\n")
		
		outfile.close()


if __name__ == "__main__":
	
	parser = argparse.ArgumentParser(description='Stampy wrapper')
	
	parser.add_argument("-b", dest="bam", required=True, help="Bam File")
	parser.add_argument("-r", dest="ref_id", help="Reference Path", default="unknown")
	parser.add_argument("-o", dest="output", help="Output file", required=True)
	parser.add_argument("-g", dest="guuid", help="Sample guuid", default="-")
	parser.add_argument("-ss", dest="seqstats", help="Seqstats Output file", required=True)
	parser.add_argument("-fs", dest="flagstats", help="Flagstats Output file", required=True)
	parser.add_argument("-subrate", dest="subrate", help="Substitution rate [0.01]", default="0.01")
	parser.add_argument("-keepgoodreads", dest="keepgoodreads", help="Use Bwa for mapping", default=False,
	                    action='store_true')
	parser.add_argument("-alignquals", dest="alignquals", help="compute post-alnmt probs", default=False,
	                    action='store_true')
	parser.add_argument("-baq", dest="baq", help="compute base-alignment quality", default=False, action='store_true')
	parser.add_argument("-singleend", dest="singleend", help="Data is single end", default=False, action='store_true')
	
	args = parser.parse_args()
	
	try:
		mapping = StampyMap(args.bam, args.output, args.ref_id, args.subrate, args.keepgoodreads, args.seqstats,
		                    args.flagstats, args.alignquals, args.baq, args.singleend, args.guuid)
		mapping.run()
		mapping.clean()
	except:
		if "mapping" in vars(): mapping.clean()
		dump_exc()
