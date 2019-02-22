# ============================================================================ #
# Carlos del Ojo Elias
# carlos.delojoelias@ndm.ox.ac.uk
# Juny 2013
# ============================================================================ #
from vcfFile import VCFFile

class GormVcf(VCFFile):

	INT4=lambda x: [int(i) for i in x.split(",")]
	FLOAT=lambda x: int(x) if not "." in x else float(x)
	INT=lambda x: int(x)

	# def FLOAT(x):
	# 	if x.strip() != '' and x.strip() != "None":
	# 		if not "." in x:
	# 			k = int(x)
	# 		else:
	# 			k = float(x)
	# 	else:
	# 		k = -99999
		# return k
	
	FORMATTER=(
		('ABQ4',FLOAT),
		('BaseCounts',INT4),
		('BaseCounts4',INT4),
		('DM4',FLOAT),
		('DM4L',FLOAT),
		('DP',INT),
		('DP4',INT4),
		('DPT4L',FLOAT),
		('DZ4',FLOAT),
		('DZ4L',FLOAT),
		('GC',FLOAT),
		('MQ',FLOAT),
		('MQ4',FLOAT),
		('PCALL4',FLOAT),
		('PCONS4',FLOAT),
		('SBR',INT)
	)


	UNCALLABLECHAR="N"
	APPROXCOVZDEPTHCHAR="-"

	# Print IUP codes with a heterozygous genotype in the input VCF file
	APPROXCOVZDEPTHCOUNT=0
	# Consider site to have approx. zero coverage if number of spanning reads DP<=threshold, and print char in FASTA file
	USEIUPFORHETS=False

	#Phred-scaled quality score for the assertion made in ALT for a variant site >= threshold [not is_indel and QUAL>=FLOAT]
	MINSNPQ=25

	#Phred-scaled quality score for the assertion made in ALT for start of an indel >= threshold [is_indel and QUAL>=FLOAT]
	MININDELQ=25

	NUC2IUP = {
		'A'    : 'A', # Adenine
		'C'    : 'C', # Cytosine
		'G'    : 'G', # Guanine
		'T'    : 'T', # Thymine
		'AG'   : 'R', # puRine
		'CT'   : 'Y', # pYrimidine
		'GT'   : 'K', # Keto
		'AC'   : 'M', # aMino
		'CG'   : 'S', # Strong
		'AT'   : 'W', # Weak
		'CGT'  : 'B', # Not A
		'AGT'  : 'D', # Not C
		'ACT'  : 'H', # Not G
		'ACG'  : 'V', # Not T
		'ACGT' : 'N'  # aNy
	}

	
	def __init__(self,fp):
		VCFFile.__init__(self,fp)
		
		self.fixHeaders()
		
	def fixHeaders(self):
		if self.complexHeaders.exists("FORMAT","GT"):
			self.complexHeaders["FORMAT"]["GT"]["Description"]='''"Genotype of the sample encoded as A/B (unphased) or A|B (phased) where 0 is the ref allele, 1 is the first alt allele, 2 is the second alt allele, X is all remaining alleles and . indicates a missing allele"'''
		if self.complexHeaders.exists("FORMAT","DP"):
			self.complexHeaders["FORMAT"]["DP"]["Description"]='''"Total number of high-quality bases, equivalent to sum(INFO:DP4)"'''
		if self.complexHeaders.exists("INFO","BaseCounts"):
			self.complexHeaders["INFO"]["BaseCounts"]["Description"]='''"Counts of each base (A,C,G,T) contributing to absolute depth INFO:DP"'''
		if self.complexHeaders.exists("INFO","DP"):
			self.complexHeaders["INFO"]["DP"]["Description"]='''"Raw read depth, the absolute depth that includes all covering reads, including those with a gap in the read at this site"'''
		if self.complexHeaders.exists("INFO","GC"):
			self.complexHeaders["INFO"]["GC"]["Type"]='Float'
		if self.complexHeaders.exists("INFO","MQ"):
			self.complexHeaders.copy("INFO","MQ","INFO","MQ4")
			self.complexHeaders["INFO"]["MQ4"]["Description"]='''"Root-mean-square mapping quality of covering reads that pass filtering thresholds"'''

	def _header(self):
		allowedTags={'INFO': ['ABQ4', 'BaseCounts', 'BaseCounts4', 'BASES4', 'BCALL', 'DM', 'DM4', 'DM4L', 'DP', 'DP4', 
							'DPT4L', 'DZ', 'DZ4', 'DZ4L', 'GC', 'INDEL', 'Iq', 'MQ', 'MQ4', 'PCALL', 'PCALL4', 'PCONS', 
							'PCONS4', 'QUALS4', 'SBR', 'Vq'],
					 'FORMAT': ['GT', 'DP']}

		res=[]

		for i,j in self.normalHeaders.items():
			for k in j:
				res.append( "##"+i+"="+k)
		res.sort()

		res2=[]
		for i in self.complexHeaders:
			if i.tag in allowedTags and i["ID"] not in allowedTags[i.tag]: continue
			else: res2.append(str(i))

		res2.sort()

		pos=None
		for i in range(len(res)):
			if res[i].startswith("##fileformat"):
				pos=i
				break

		if pos==None: res.insert(0,"##fileformat=VCFv4.1")
		else: 
			res.insert(0,res.pop(pos))

		for i in res+res2:
			yield i

	def getParsedRecord(self,contig,pos,indel=False):
		r=VCFFile.getParsedRecord(self,contig,pos,indel)
		for i,j in self.FORMATTER:
			if i in r[self.fld_PARSED_INFO]:
				# print r[self.fld_PARSED_INFO] #Debug
				r[self.fld_PARSED_INFO][i]=j(r[self.fld_PARSED_INFO][i])
		return r

	def iterBody(self,parsed=False):
		if not parsed:
			for i in self._body():
				yield i
		else:
			for record in self._body():
				record.append(dict([i.split("=") for i in record[self.fld_INFO].split(";") if "=" in i]))
				record.append(dict(zip(record[self.fld_FORMAT].split(":"),record[self.fld_DATA].split(":"))))
				for i,j in self.FORMATTER:
					if i in record[self.fld_PARSED_INFO]:
						record[self.fld_PARSED_INFO][i]=j(record[self.fld_PARSED_INFO][i])
				record[self.fld_QUAL]=float(record[self.fld_QUAL])
				yield record


	def getAdditionalInfo (self,contig,pos,pileupResource):
		'''Generates extra info mixing information from a mpileup output'''
		indelCall=None
		indel=None
		baseCall=self.APPROXCOVZDEPTHCHAR

		snp=self.getParsedRecord(contig,pos)

		if self.isindel(contig,pos):
			indel=self.getParsedRecord(contig,pos,indel=True)

		# GTBASES contain the bases of the genotype (different alleles)
		gt=set([int(i) for i in snp[self.fld_PARSED_DATA].setdefault("GT","./.") if i in "0123456789"])
		bases=[snp[self.fld_REF].upper()]+snp[self.fld_ALT].upper().split(",")
		snp[self.fld_PARSED_DATA]["GTBASES"]="".join([bases[i] for i in gt if bases[i] in ('A','C','T','G',GormVcf.UNCALLABLECHAR,GormVcf.APPROXCOVZDEPTHCHAR)])

		# GTDEPTH contains the number of reads contributing to the alleles in GTBASES
		bcounts=snp[self.fld_PARSED_INFO]["BaseCounts"]
		snp[self.fld_PARSED_DATA]["GTDEPTH"]=sum([bcounts["ACGT".index(i)] for i in snp[self.fld_PARSED_DATA]["GTBASES"] if i in "ACGT"])

		baseCall=self.basecall(snp)
		if indel:
			indelCall=self.basecall(indel,True)

		# GTHQDEPTH contains the number of high quality reads contributing to the alleles in GTBASES
		snp[self.fld_PARSED_DATA]["BCALLDEPTH"]=sum([bcounts["ACGT".index(i)] for i in baseCall if i in "ACGT"])

		if not (contig,pos) in pileupResource:
			snp[self.fld_PARSED_DATA]["GTHQDEPTH"]=0
			return [snp,indel,GormVcf.APPROXCOVZDEPTHCHAR,GormVcf.APPROXCOVZDEPTHCHAR,[],[],[0,0,0,0],False]


		hqbases=[ i for i in pileupResource.basesAndQuals(contig,pos) if i[0]!="*" and i[1]>=self.MINSNPQ ]
		if not hqbases:
			snp[self.fld_PARSED_DATA]["GTHQDEPTH"]=0
			return [snp,indel,baseCall,indelCall,[],[],[0,0,0,0],True]

		bases4,quals4=zip(*hqbases)
		capitalbases4="".join(bases4).upper()
		basecounts4=[capitalbases4.count("A"),capitalbases4.count("C"),capitalbases4.count("G"),capitalbases4.count("T")]

		# GTHQDEPTH contains the number of high quality reads contributing to the alleles in GTBASES
		bhqcounts=basecounts4
		snp[self.fld_PARSED_DATA]["GTHQDEPTH"]=sum([bhqcounts["ACGT".index(i)] for i in snp[self.fld_PARSED_DATA]["GTBASES"] if i in "ACGT"])

		return [snp,indel,baseCall,indelCall,bases4,quals4,basecounts4,True]


	def basecall(self, record,isindel=False):

		'''
		Record indexing is a requirement

		Return the base call inferred from the information for a normal or indel site.
		
		If the absolute number of spanning reads (DP) is almost zero (less than or equal
		to the approx covz depth count from the gorm.ini file) then return a base call
		of the uncallablechar from the gorm.ini file.
		
		If the alternate allele is the same as the reference and there are no spanning reads,
		then return the approxcovzdepthchar.
		
		If we are trying to call an INDEL and there is no genotype (GT contains A or B
		is '.' for one of the A/B genotypes), then assume the call is the same as the
		reference.
		
		If we have a homozygous call, then output the base or the IUP code.
		
		Otherwise, we must have a heterozygous call. For variants, output the IUP code
		if requested, or the uncallablechar.
		
		Note: Calls to this function for INDEL sites should always set useiupforhets=False
		and isindel=True.
		'''

		ref=record[VCFFile.fld_REF]
		alt=record[VCFFile.fld_ALT]
		gt=record[VCFFile.fld_PARSED_DATA].setdefault("GT",'./.')
		dp=int(record[VCFFile.fld_PARSED_INFO].setdefault("DP",0))

		if not ref or not alt or not gt: return self.UNCALLABLECHAR
		if dp<self.APPROXCOVZDEPTHCOUNT: return self.APPROXCOVZDEPTHCHAR
		if alt == '.' and not dp: return self.APPROXCOVZDEPTHCHAR

		gt_alleles = [ref] + alt.split(',')
		X, Y = gt.split("/")

		# This is to catch INDEL lines which have alt='.' and no FORMAT:GT tag
		if isindel and (X == '.' or Y == '.'): return ref
		
		# This is to catch homozygous calls
		if X==Y: return gt_alleles[int(X)]

		# This is for heterozygous calls
		if self.USEIUPFORHETS and not isindel:
			nuc_str = [gt_alleles[int(X)],gt_alleles[int(Y)]]
			nuc_str.sort()
			return self.NUC2IUP.setdefault("".join(nuc_str),self.APPROXCOVZDEPTHCHAR)

		return self.UNCALLABLECHAR



if __name__=="__main__":
#	a=GormVcf("tmp/C00002880_R00000003.mapcall.5.varannotDP.vcf")
#	b=GormVcf("tmp/C00002880_R00000003.mapcall.5.varannot.vcf")
#	c=GormVcf("tmp/C00002880_R00000003.mapcall.6.indels.vcf")
#
#	for i in b._header():
#		a.addHeader(i)
#
#	for i in c._header():
#		a.addHeader(i)
#			
#	for i in a:
#		print i

	from mpileup import PileupReader

	# p=PileupReader("tmp/C00002880_R00000003.mapcall.3.mpileup.txt",index=True)

	# a=GormVcf("tmp/C00002880_R00000003.mapcall.7.merged.vcf")
	# a=VCFFile("/Users/thanhlv/Pipeline/compass/0aa59cfd-12f2-4cd4-9d13-52c1e931a320.annotvcf.vcf")
	a=VCFFile("/Users/thanhlv/Pipeline/compass/0aa59cfd-12f2-4cd4-9d13-52c1e931a320.annotvcf.vcf")
	
	# print a.getRecord('NC_012578.1',19)
	# for i in range(10,20,1):
	# 	print a.getParsedRecord('NC_012578.1',i)
	
	#
	# print a.parseHeaders()
	# print a.snpIndex
	# i = 0
	# for i in a.indexRecords():
	# 	# pass
	# 	print i
	# 	i += 1
	# 	if i >
	# print a.getChromosomes()
	# for i in range(5000000):
	# 	try:
	# 		a.getAdditionalInfo("R00000003",i,p)
	# 	except VCFFile.RecordNotExistent: pass
