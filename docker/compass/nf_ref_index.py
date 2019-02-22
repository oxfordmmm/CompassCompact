#!/usr/bin/env python
import sys
import re
import os
from lib.logerror import LE, dump_exc
from Bio import SeqIO
import shlex
from StringIO import StringIO
from Bio.Blast.Applications import NcbiblastnCommandline as blastn
from Bio.Application import ApplicationError as AppError
import numpy as np
from itertools import izip
import uuid
import shutil
import argparse
from lib.compassconfig import *


class RefGen:
	"""
    """
	
	def __init__( self, infasta):
		"""

        Args:
            infasta:
        """
		self.newrefid = os.path.basename(infasta).split('.')[0]
		self.infasta = infasta
		self.contigs = 0
		self.numbp = 0
		self.outdir = os.path.dirname(infasta)
			
		self.outfastapath = os.path.join(
			self.outdir, '{0}_new.fasta'.format(self.newrefid))
		self.outfile = None
		
		self.preprocess()
	
	def preprocess( self ):
		"""

        """
		fasta = SeqIO.parse(self.infasta, "fasta")
		for i in fasta:
			self.contigs += 1
			self.numbp += len(i.seq)
		fasta.close()
	
	def cleanTempFiles( self ):
		"""

        """
		try:
			shutil.rmtree(self.outdir)
		except:
			pass
	
	def cleanAll( self ):
		"""

        """
		self.cleanTempFiles()
		try:
			os.unlink(self.outfile)
		except:
			pass
	
	def generateAll( self ):
     		self.Fix_Fasta_Headers()
        	self.Create_Indexes()
        	self.Make_Repeat_Mask_Txt()
        	LE.info("Everything went OKAY!")
	
	def createFile( self, path ):
		"""

        Args:
            path:

        Returns:

        """
		if not os.path.exists(path):
			os.mkdir(path)
		if not os.path.isdir(path):
			raise Exception("{0} is not a directory")
		
		self.outfile = os.path.join(path, self.newrefid + ".tar.bz2")
		
		try:
			self.execCommand(
				"tar -cj -C {0} -f {1} .".format(self.outdir, self.outfile))
		except Exception, e:
			self.cleanAll()
			raise e
		
		return self.outfile
	
	#Copy output.bam from tmp folder to location defined by output.bam parameter
	def CopyOutput(self, path):
		try:
			shutil.copytree(self.outdir, path)
		except Exception as e:
			raise e
	
	def execCommand( self, cmd ):
		"""

        Args:
            cmd:
        """
		origcmd = cmd
		LE.debug("Running command: [{0}]".format(cmd))
		cmd = shlex.split(cmd)
		p = subprocess.Popen(cmd, stdout=subprocess.PIPE,
		                     stderr=subprocess.PIPE)
		stdout, stderr = p.communicate()
		errcode = p.wait()
		LE.debug(StringIO(stdout))
		if errcode:
			LE.error(StringIO(stderr))
			raise Exception(
				"CMD [{0}] exit with status [{1}]".format(origcmd, errcode))
	
	def cleanUpExecution( self, cmd, stdout, stderr, errcode ):
		"""

        Args:
            cmd:
            stdout:
            stderr:
            errcode:
        """
		LE.debug(StringIO(stdout))
		if errcode:
			LE.error(StringIO(stderr))
			raise Exception(
				"CMD [{0}] exit with status [{1}]".format(cmd, errcode))
	
	def Fix_Fasta_Headers( self ):
		"""
        Create correctly formatted fasta file. Contigs must be in the form REFID[, REFID-2, REFID-3,...].
        """
		
		LE.info('Creating master fasta file {0}.'.format(self.outfastapath))
		
		# helper function to reformat each fasta record on the fly
		def _fixed_records():
			"""

            """
			for i, contig in enumerate(SeqIO.parse(self.infasta, 'fasta'), 1):
				correct_name = contig.id = self.newrefid + \
				                           ('-{0}'.format(i), '')[i == 1]
				if contig.name != correct_name:
					contig.name = correct_name
					contig.description = '{0} {1} {2}'.format(
						correct_name, self.newrefid, contig.description)
				yield contig
		
		SeqIO.write(_fixed_records(), self.outfastapath, 'fasta')
		return
	
	def Index_Fasta( self ):
		"""
        Create/update .fasta.fai file using samtools faidx.
        """
		LE.info('Creating index {0}.fai.'.format(self.outfastapath))
		self.cleanUpExecution(
			*COMPASSCFG['tools']['samtools'].execute(append="faidx {0}".format(self.outfastapath)))
	
	def Create_Symlinks( self ):
		"""
        Create symbolic links to REFID.fasta as REFID.fa (plus indexes); required for some third-party tools.
        """
		prefix = os.path.abspath(os.path.join(self.outdir, self.newrefid))
		outfastapath = os.path.basename(prefix + '.fasta')
		outsymlink = prefix + '.fa'
		for suffix in ('', '.fai'):
			src = outfastapath + suffix
			dst = outsymlink + suffix
			os.symlink(src, dst)
		return
	
	def Create_Indexes( self ):
		"""
        Create stampy and BWA index files.
        """
		prefix = os.path.join(self.outdir, self.newrefid)
		outfastapath = prefix + '.fasta'
		stidxpath = prefix + '.stidx'
		sthashpath = prefix + '.sthash'
		bwaidxpath = outfastapath + '.pac'
		nhrpath = prefix + '.nhr'
		
		self.cleanUpExecution(*COMPASSCFG['tools']['stampy'].execute(
		    append="-G {prefix} {outfastapath}".format(prefix=prefix, outfastapath=outfastapath)))
		self.cleanUpExecution(
		    *COMPASSCFG['tools']['stampy'].execute(append="-g {prefix} -H {prefix}".format(prefix=prefix)))
		# BWA
		self.cleanUpExecution(
			*COMPASSCFG['tools']['bwa'].execute(append="index -a is {outfastapath}".format(outfastapath=outfastapath)))
		
		# blasdb
		self.cleanUpExecution(*COMPASSCFG['tools']['blast'].execute(source='path', file="makeblastdb",
		                                                            append="-dbtype nucl -in {outfastapath} -title {refid} -out {prefix}".format(
			                                                            outfastapath=outfastapath, refid=self.newrefid,
			                                                            prefix=prefix)))
		
	def Make_Repeat_Mask_Txt( self, word_size=17, gapopen=5, e_thresh=0.0001, perc_identity=90, gapextend=2,
	                          min_length=75 ):
		"""
        Run blastn on contigs in input fasta file against database dbname. Parameters set to NCBI recommended defaults for blastn.
        """
		outfastapath = os.path.join(
			self.outdir, '{0}.fasta'.format(self.newrefid))
		prefix = os.path.join(self.outdir, self.newrefid)
		maskpath = prefix + '_repmask.array'
		regionspath = prefix + '_repregions.array'
		statspath = prefix + '.stats'
		
		blastn_cline = blastn(cmd=COMPASSCFG['tools']['blast']['path'] + "blastn", db=prefix, query=outfastapath,
		                      dust='no', word_size=word_size, gapopen=gapopen, gapextend=gapextend, evalue=e_thresh,
		                      perc_identity=perc_identity,
		                      outfmt='"6 qseqid sseqid pident length qstart qend sstart send"')
		try:
			blast_out, blast_err = blastn_cline()
			assert not blast_err
		except (AppError, AssertionError) as err:
			raise Exception(
				'Erro: Blast failed during construction of repeat mask : {0}'.format(err))
		
		repmask_fp = open(maskpath, 'w')
		repregions_fp = open(regionspath, 'w')
		total_bp = 0
		repetitive_bp = 0
		num_regions = 0
		
		# each blast_rec is result from one query sequence (contig)
		blast_stream = StringIO(blast_out)
		prev_header = None
		for contig_count, contig in enumerate(SeqIO.parse(outfastapath, 'fasta'), 1):
			if prev_header != contig.name:
				repregions_fp.write('>{0}\n'.format(contig.name))
				prev_header = contig.name
			total_bp += len(contig)
			repmask = np.zeros(len(contig), dtype=np.bool)
			try:
				fields = blast_stream.next().split()
			except StopIteration:
				fields = None
			while fields and fields[0] == contig.name:
				contig_name, match_name = fields[:2]
				hit_perc_ident = float(fields[2])
				hit_length, q_start, q_end, s_start, s_end = (
					int(x) for x in fields[3:])
				(x1, y1), (x2, y2) = sorted(
					((q_start, q_end), sorted((s_start, s_end))))
				if hit_length >= min_length and (contig_name != match_name or not (x2 <= x1 <= y2 and x2 <= y1 <= y2)):
					repmask[q_start - 1:q_end] = True
				try:
					fields = blast_stream.next().split()
				except StopIteration:  # end of blast hits
					fields = None
			# output.bam repmask as 1 and 0, 100 per line
			repmask_fp.write('>{0}\n'.format(contig.name))
			for i in xrange(0, len(repmask), 100):
				j = min(i + 100, len(repmask))
				repmask_fp.write('{0}\n'.format(''.join(str(i)
				                                        for i in repmask[i:j].astype(int))))
			# identify postitions of repetitive regions (runs of 1s in the
			# repmask array)
			# 0-based numbering
			region_starts = list(np.where(repmask[1:] > repmask[:-1])[0] + 1)
			region_ends = list(np.where(repmask[1:] < repmask[:-1])[0] + 1)
			# special case: full blast hit for this contig against another
			# contig
			if repmask.all():
				region_starts = [0]
				region_ends = [len(repmask)]
			# fix ends, in case regions start from the first position in the
			# sequence or end at the last
			if region_starts and ((not region_ends) or (region_starts[-1] > region_ends[-1])):
				region_ends.append(len(repmask))
			if region_ends and ((not region_starts) or (region_starts[0] > region_ends[0])):
				region_starts = [0] + region_starts
			repregions_fp.writelines('{0}\t{1}\n'.format(
				rs, re) for rs, re in izip(region_starts, region_ends))
			repetitive_bp += repmask.sum()
			num_regions += len(region_starts)
		
		repmask_fp.close()
		repregions_fp.close()
		pct_repetitive = '{0:.2f}'.format(
			(float(repetitive_bp) / total_bp) * 100)
		LE.debug(
			'Info: Repetitive regions for all of {0}: {1}/{2} bp ({3}%)'.format(self.newrefid, repetitive_bp, total_bp,
			                                                                    pct_repetitive))
		
		# save result summary
		statsvalues = '\t'.join((self.newrefid, self.newrefid, str(contig_count), str(total_bp), str(repetitive_bp),
		                         str(num_regions), pct_repetitive))
		with open(statspath, 'w') as o:
			o.write('refid\trefcd\tcontigs\tnumbp\trepetitivebp\trepregions\trepetitivepct\n{values}\n'.format(
				values=statsvalues))
		return
		

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Compass pileup wrapper')
	
	parser.add_argument("-r", dest="ref", required=True,
	                    help="reference Fasta file")
	
	options = parser.parse_args()
	
	a = RefGen(options.ref)
	a.generateAll()

