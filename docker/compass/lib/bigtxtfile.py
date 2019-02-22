# Gorm4 - BigTxtFile reader
#
# This class abstracts bigfile reading, detecting carriage return length and
# avoding strip(). It can split by fields if the file has fielded information
# It detects header if it is provided in the constructor
#
# Provides an abstract method initialise() to do further operations anfter
# reading the header in case the class is inherited
#
# ============================================================================ #
# Carlos del Ojo Elias
# carlos.delojoelias@ndm.ox.ac.uk
# Juny 2013
# ============================================================================ #

import itertools
from gzip import GzipFile


class BigTxtFile:
    def __init__(self, fp, inibuff=100, header_preffix=None, split=None):
        '''
BigTextFile(file_obj[, inibuff=100, header_preffix=None, split=None]) -> BigTextFile

fp - must be a file object (seek+tell methods)

inibuff - Initial buffer in lines to be read and detect CRlen

header_preffix - If provided inibuff will be discarded and all the lines
starting with header_preffix will be treated as a different part of the file
(header)

split - If provided, every line after header will be split by the provided
delimiter'''

        assert type(fp) in [str, file]
        if type(fp) == str:
            try:
                self.fp = GzipFile(fp)
                self.fp.readline()
                self.fp.seek(0)
            except:
                self.fp = open(fp)
        else:
            self.fp = fp

        self.header = []
        self.inibuff = []
        self.lenCR = 1
        self.split = split

        self.bodypos = 0

        if header_preffix:
            for i in self.fp:
                self.bodypos += len(i)
                if len(i.strip()) != len(i):
                    self.lenCR = len(i) - len(i.strip())
                if i.startswith(header_preffix):
                    self.header.append(i.strip())
                else:
                    self.bodypos -= len(i)
                    break

        else:
            for i in self.fp:
                self.inibuff.append(i)
                inibuff -= 1
                if inibuff < 1:
                    break

        if self.inibuff:
            self.lenCR = len(self.inibuff[0]) - len(self.inibuff[0].strip())

        self._initialise()

    def _initialise(self):
        '''Abstract method in case you want postprocess headers or indexing'''
        pass

    def _body(self):
        '''Generator for each line of the body (whole file if header not provided)'''
        self.fp.seek(self.bodypos)
        if self.split:
            for i in self.fp:
                yield i[:-self.lenCR].split(self.split)
        else:
            for i in self.fp:
                yield i[:-self.lenCR]

    def _header(self):
        '''Generator for the header records'''
        for i in self.header:
            yield i

    def __iter__(self):
        '''Generator for the whole file'''
        for i in itertools.chain(self._header(), self._body()):
            yield i
