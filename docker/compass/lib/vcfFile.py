# ============================================================================ #
# Carlos del Ojo Elias
# carlos.delojoelias@ndm.ox.ac.uk
# Juny 2013
#
# Time indexing a file 31s [~4.292.129 records]
# Time for getParsedRecord(0,5.000.000): 1m26s
# ============================================================================ #

from bigtxtfile import BigTxtFile
import re
from copy import deepcopy
import array
from extendibleArray import ExtendibleArray
from genomemap import GenomeMap
from StringIO import StringIO

### TODO 
# WARN PARSING ERRORS
# Test using dictionaries for fields so you will not have ordered attributes (maybe it works)
# Exception managing
# Documentation

class VCFFile(BigTxtFile):
	'''VCFFile(filename) -> VCFfile'''

	### We consider two kind of headers, normal and complex
	### Complex ones are those which accept multiple fields: ##FORMAT=<ID=DP,Number=1,Type=Integer,Des...

	### Headers containing multiple fields
	RecordHeaders=["##INFO","##FILTER","##FORMAT","##ALT","##contig","##SAMPLE","##PEDIGREE"]

	### Regular expressions to parse the headers
	re_normalHeader=re.compile("##([a-z0-9_-]+)=(.*)",re.I)
	re_complexHeader=re.compile("##([a-z0-9_-]+)=<(.*)>",re.I)
	re_varnameQuot=re.compile('''^[a-z0-9]+=['"]''',re.I)

	fld_CHROM=0
	fld_POS=1
	fld_ID=2
	fld_REF=3
	fld_ALT=4
	fld_QUAL=5
	fld_FILTER=6
	fld_INFO=7
	fld_FORMAT=8
	fld_DATA=9
	fld_PARSED_INFO=10
	fld_PARSED_DATA=11

	class NoIDinComplexHeader(Exception):
		pass

	class RecordNotExistent(Exception): 
		pass

	class VcfComplexHeader:
		'''Class that abstracts complex headers like: ##FORMAT=<ID=DP,Number=1,Type=Integer,Des...
		It requires a string containig the header 
		'''

		def __init__(self,line):
			# Get the fileds by unsing a regexp
			vals=VCFFile.re_complexHeader.findall(line)

			if vals:
				tag,info=vals[0]
				info=info.split(",")
				info2=[]
				for i in info:
					if not info2:
						info2.append(i)
					elif (VCFFile.re_varnameQuot.findall(info2[-1]) and info2[-1].strip()[-1] not in "'\""):
						info2[-1]+=","+i
					else:
						info2.append(i)
				info2=[i.split("=",1) for i in info2]

			self.tag=tag

			# sortedkeys contain the original order for the keys found in the input string
			self.sortedkeys=[i[0] for i in info2]
			# info contains the information asociated to avery key
			self.info=dict(info2)

			# Every complex header must have an ID field
			if not "ID" in self.info:
				raise VCFFile.NoIDinComplexHeader()

		def __setitem__(self,k,v):
			if k not in self.sortedkeys: self.sortedkeys.append(k)
			self.info[k]=v

		def __getitem__(self,k):
			return self.info[k]

		def __str__(self):
			return "##"+self.tag+"=<"+",".join(["=".join([i,self.info[i]]) for i in self.sortedkeys])+">"

	class VcfComplexHeaderCollection:
		'''VcfComplexHeaderCollection abstracts a collection of different complex headers'''
		def __init__(self):
			self.pointers={}
			self.sortedEntries=[]

		def addComplexHeader(self,line):
			try:
				vch=VCFFile.VcfComplexHeader(line)
			except VCFFile.NoIDinComplexHeader:
				return											########### REVISE

			if not vch.tag in self.pointers or not vch["ID"] in self.pointers[vch.tag]:
				self.sortedEntries.append(vch)

			self.pointers.setdefault(vch.tag,{})[vch["ID"]]=vch

		def exists(self,tag,ID):
			return tag in self.pointers and ID in self.pointers[tag]

		def copy(self,tag,ID,tag2,ID2):
			assert self.exists(tag,ID) 
			new=deepcopy(self.pointers[tag][ID])
			new["ID"]=ID2
			if not self.exists(tag2,ID2):
				self.sortedEntries.append(new)
			self.pointers.setdefault(tag2,{})[ID2]=new

		def __iter__(self):
			for i in self.sortedEntries:
				yield i

		def __getitem__(self,k):
			return self.pointers[k]

	def __init__(self,fp):
		self.snpIndex= None
		self.indelIndex=None
		
		self.normalHeaders={}
		self.complexHeaders=VCFFile.VcfComplexHeaderCollection()
		self.indels=GenomeMap()
		BigTxtFile.__init__(self,fp,header_preffix="#",split="\t")
		
	def _initialise(self):
		self.parseHeaders()

	def parseHeaders(self):
		for i in self.header:
			if i.startswith("##"):
				self.addHeader(i)
	
	def addHeader(self,i):
			if any([i.startswith(j) for j in VCFFile.RecordHeaders]):
				self.complexHeaders.addComplexHeader(i)	
			else: 
				self.addNormalHeader(i)

	def addNormalHeader(self,line):
		vals=VCFFile.re_normalHeader.findall(line)
		if vals: 
			self.normalHeaders.setdefault(vals[0][0],set()).add(vals[0][1])

	def getNormalHeaders(self):
		return self.normalHeaders

	def _header(self):
		for i,j in self.normalHeaders.items():
			yield "##"+str(i)+"="+str(j)
		for i in self.complexHeaders:
			yield str(i)


	def indexRecords(self):
		self.snpIndex={}
		self.indelIndex={}
		
		self.fp.seek(self.bodypos)
		pos=self.bodypos
		for i in self.fp:
			if i[0]!="#":
				len_record=len(i)
				i=i.split() #Split by tab \t
				refpos=int(i[1])

				# We are indexing SNPS and INDELS as they can be in the same position of the same contig
				if "INDEL" in i[7]:
					if i[0] not in self.indelIndex:
						self.indelIndex[i[0]]=ExtendibleArray('l',defvalue=-1)
					self.indelIndex[i[0]].insert(refpos,pos)
				else:
					if i[0] not in self.snpIndex:
						self.snpIndex[i[0]]=ExtendibleArray('l',defvalue=-1)
					self.snpIndex[i[0]].insert(refpos,pos)
				yield i
			
			pos+=len_record


	def cacheRegion(self,positions,chromosome=None):
		self.cachedVcfContent=StringIO()
		positions=set([str(i) for i in positions])
		if not chromosome:
			for i in self._body():
				if i[1] in positions:
					self.cachedVcfContent.write("\t".join(i)+"\n")
		else:
			for i in self._body():
				if i[0]==chromosome and i[1] in positions:
					self.cachedVcfContent.write("\t".join(i)+"\n")

		self.fp.close()
		self.fp=self.cachedVcfContent
		self.fp.seek(0)
		self.bodypos=0
		for i in self.indexRecords():
			pass

	def iterContigsPos(self):
		contigs=self.snpIndex.keys()
		contigs.sort()
		for c in contigs:
			for pos in range(len(self.snpIndex[c])):
				if self.snpIndex[c][pos]>=0: yield c,pos

	def isindel(self,contig,pos):
		try:
			return self.indelIndex[contig][pos]!=-1
		except: return False

	def getValue(self,regexp,typ,mod=None,indel=False):

	# mod can be a method/lambdfa function to modify or preprocess the data

	#	'c'	char	character	1
	#	'b'	signed char	int	1
	#	'B'	unsigned char	int	1
	#	'u'	Py_UNICODE	Unicode character	2 (see note)
	#	'h'	signed short	int	2
	#	'H'	unsigned short	int	2
	#	'i'	signed int	int	2
	#	'I'	unsigned int	long	2
	#	'l'	signed long	int	4
	#	'L'	unsigned long	long	4
	#	'f'	float	float	4
	#	'd'	double	float	8

		destobj=lambda: ExtendibleArray(typ)
		dest={}

		regexp=re.compile(regexp)

		if not mod:
			if typ in 'hHiIlL':
				for i in self._body():
					if i[0] not in dest: dest[i[0]]=destobj()
					pos=int(i[1])
					res=regexp.findall(i[7])
					if not indel and "INDEL" in i[7]: continue
					if res: dest[i[0]].insert(pos,int(res[0]))
			elif typ in 'fd':
				for i in self._body():
					if i[0] not in dest: dest[i[0]]=destobj()
					pos=int(i[1])
					res=regexp.findall(i[7])
					if not indel and "INDEL" in i[7]: continue
					if res: dest[i[0]].insert(pos,float(res[0]))
			else:
				for i in self._body():
					if i[0] not in dest: dest[i[0]]=destobj()
					pos=int(i[1])
					res=regexp.findall(i[7])
					if not indel and "INDEL" in i[7]: continue
					if res: dest[i[0]].insert(pos,res[0])
		else:
			for i in self._body():
				if i[0] not in dest: dest[i[0]]=destobj()
				pos=int(i[1])
				res=regexp.findall(i[7])
				if not indel and "INDEL" in i[7]: continue
				if res: dest[i[0]].insert(pos,mod(res[0]))

		return dest

	def getRecord(self,contig,pos,indel=False):
		try:
			if indel==True:
				pos=self.indelIndex[contig][pos]
			else:
				pos=self.snpIndex[contig][pos]
			self.fp.seek(pos)
			record=self.fp.readline()[:-self.lenCR].split(self.split)
			
			return record

		except:
			raise VCFFile.RecordNotExistent
		
	# def getParsedRecord(self,contig,posic,indel=False):
	# 	#Fix splitting when QUALS4 has ; character. If so, the will be a case like: "QUALS4=;GC=HGGHT"
	# 	record=self.getRecord(contig,posic,indel)
	# 	if (not indel) and ("QUALS4" in record[self.fld_INFO]):
	# 		q = re.compile(r";QUALS4=[\w.*;<=>:\"\'-@?]+")
	# 		_tmp = q.search(record[self.fld_INFO]).group()[1:]
	# 		_rec = [i.split("=",1) for i in q.sub("",record[self.fld_INFO]).split(";") if "=" in i]
	# 		_rec.append(_tmp.split("=",1))
	# 		record.append(dict(_rec))
	# 	else:
	# 		record.append(dict([i.split("=",1) for i in record[self.fld_INFO].split(";") if "=" in i]))
	# 	record.append(dict(zip(record[self.fld_FORMAT].split(":"),record[self.fld_DATA].split(":"))))
	# 	return record
	
	def getParsedRecord( self, contig, posic, indel=False ):
			record = self.getRecord(contig, posic, indel)
			record.append(dict([i.split("=", 1) for i in record[self.fld_INFO].split(";") if "=" in i]))
			record.append(dict(zip(record[self.fld_FORMAT].split(":"), record[self.fld_DATA].split(":"))))
			return record
	
	def parsedRecords(self):
		for record in self._body():
			record.append(dict([i.split("=") for i in record[self.fld_INFO].split(";") if "=" in i]))
			record.append(dict(zip(record[self.fld_FORMAT].split(":"),record[self.fld_DATA].split(":"))))
			yield record


	def getChromosomes(self):
		chrs={}
		for i,j in self.snpIndex.items():
			chrs.setdefault(i,0)
			chrs[i]=max(chrs[i],len(j))
		for i,j in self.indelIndex.items():
			chrs.setdefault(i,0)
			chrs[i]=max(chrs[i],len(j))

		return chrs
		


if __name__=="__main__":
#	a=VCFFile("tmp/C00002880_R00000003.mapcall.7.merged.vcf")

#
#	for i in range(5000000):
#		try:
#			a.getParsedRecord("R00000003",i)
#		except VCFFile.RecordNotExistent: pass

	# a=VCFFile("/Users/thanhlv/Pipeline/compass/0aa59cfd-12f2-4cd4-9d13-52c1e931a320.annotvcf.vcf")
	a=VCFFile("/Users/thanhlv/Pipeline/compass/experiments/28_S25_merged.annotvcf.vcf.gz")
	for i in a.indexRecords():
		pass
	
	x = a.getRecord('NC_012578.1', 98977)
	
	# print a.getParsedRecord('NC_012578.1', 98977)
 
#	import code
#	code.interact(local=locals())
