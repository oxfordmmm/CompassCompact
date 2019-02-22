## Gorm4 - FastaFile
#
# This class implements a FastA wrapper
# When the object is created, it parses wht whole fasta file indexing positions
# for every contig
#
# Features:
#	* Information per contig -> length
#	* Fast access to any position
#	* Subsequence extraction
#
# ============================================================================ #
# Carlos del Ojo Elias
# carlos.delojoelias@ndm.ox.ac.uk
# Juny 2013
# ============================================================================ #

import re
import gzip
import bz2
from logerror import LE, dump_exc


class FastaFile():
    '''FastaFile (filename) -> FastaFile()

        Detects automatically wether it is .gz or .bz2
        It is strongly recommended to use uncompressed streams
    '''

    def __init__(self, fil):
        self.fafile = fil
        extension = self.fafile.split(".")[-1].lower()
        tempf = None
        self.chromosomes = {}
        self.chlengths = {}

        # Detects wether i
        if extension in ['gz', 'bz2']:
            tempf = tempfile.NamedTemporaryFile("w+")
            if extension == "gz":
                self.chf = gzip.GzipFile(self.fafile)
            elif extension == "bz2":
                self.chf = bz2.BZ2File(self.fafile)
        else:
            self.chf = open(fil)

        currentChromosome = None
        pos = 0
        firstpart = re.compile("([^ \t,]+).*")  # TODO: explain search what?

        # Indexin every chromosome and the line length (column width)
        i = self.chf.readline()
        leni = 0
        while i:
            if tempf:
                tempf.write(i)
            if i.startswith(">"):
                currentChromosome = firstpart.findall(i[1:].strip())[0]
                self.chlengths[currentChromosome] = 0
                pos = 0
                leni = 0
            elif i.strip():
                self.chlengths[currentChromosome] += len(i.strip())
                if len(i) != leni:
                    leni = len(i)
                    self.chromosomes.setdefault(currentChromosome, []).append(
                        [pos, self.chf.tell() - len(i), leni, len(i.strip())])
                pos += len(i.strip())
            i = self.chf.readline()

        if tempf:
            self.chf.close()
            self.chf = tempf

            #		### Keep the first word of the chromosome, avoiding further description found in the ID
            #		for i in self.chromosomes.keys():
            #			short=firstpart.findall(i)[0]
            #			self.chromosomes.setdefault(short,self.chromosomes[i])
            #			self.chlengths.setdefault(short,self.chlengths[i])
            #
            #		for i in oldkeys:
            #			del self.chromosomes[i]
            #			del self.chlengths[i]

    def getChromosomes(self):
        '''Returns a dictionary of the chromosome names and its lengths'''
        return self.chlengths

    def __iter__(self):
        '''Iterates over the whole file'''
        for i in self.chromosomes:
            yield (i, self.getFasta(i))

    def getFasta(self, ch, reversecomp=False):
        '''Returns a generator for the fasta seqence (including the >identifier)
            Requires a chromosome identifier and allows to iterate over the revese
            complemented sequence
        '''

        realpos, filepos = self.getLowPos(ch, 0)
        self.chf.seek(filepos)
        yield (">" + ch + "\n")

        if not reversecomp:
            while True:
                res = self.chf.readline()
                if res.startswith(">") or not res:
                    break
                yield res
        else:
            for i in range(self.chlengths[ch] - 100, -1, -100):
                yield str(Seq(self.getChunk(ch, i, 100)).reverse_complement())
            yield str(Seq(self.getChunk(ch, 0, self.chlengths[ch] % 100)).reverse_complement()) + "\n"

    def getLowPos(self, chrom, pos):
        for i in range(len(self.chromosomes[chrom])):
            if self.chromosomes[chrom][i][0] > pos:
                i -= 1
                break

        start, fpos, linelen, dnalen = self.chromosomes[chrom][i]

        jump = ((pos - start) / dnalen)
        return jump * dnalen + start, jump * linelen + fpos

    def getChunk(self, chrom, pos, size):
        '''Returns a chunk of any chromosome given the chrom ID, position and spanning length '''
        if pos < 0:
            pos = self.chlengths[chrom] + pos
        realpos, filepos = self.getLowPos(chrom, pos)
        vsize = size + pos - realpos
        self.chf.seek(filepos)
        data = ""
        while True:
            l = self.chf.readline()
            if not l:
                break
            if l[0] == ">":
                break
            if not l:
                break
            data += l.strip()
            if len(data) > vsize:
                break
        return data[pos - realpos:pos - realpos + size]

    def close(self):
        self.chf.close()

    def getRandomChunk(self, ch, chunk_sample_len):
        lgth = self.chlengths[ch]
        if lgth <= chunk_sample_len:
            chunk_sample = self.getChunk(ch, 0, lgth)
        else:
            pos = random.randint(0, lgth - chunk_sample_len)
            chunk_sample = self.getChunk(ch, pos, chunk_sample_len)
        return chunk_sample
