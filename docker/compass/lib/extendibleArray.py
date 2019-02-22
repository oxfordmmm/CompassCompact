## Gorm4 - ExtendibleArray
#
# This class extends array class,but it lets to add elements at any position
# in the array, expanding it automatically.
#
# ============================================================================ #
# Carlos del Ojo Elias
# carlos.delojoelias@ndm.ox.ac.uk
# Juny 2013
# ============================================================================ #

import array
from genomemap import ContigMap


class ExtendibleArray(array.array):
    class DefValue:
        '''Provides a fast zero iterator the the array.extend method'''

        def __init__(self, n, defvalue=0):
            self.n = n
            self.defvalue = defvalue

        def __iter__(self):
            k = self.n
            while k:
                yield self.defvalue
                k -= 1

    def __init__(self, typ, defvalue=0):
        self.reallen = 0
        self.lenarray = self.__len__()
        self.defval = defvalue
        self.setPositions = ContigMap()
        array.array.__init__(self, typ)

    def insert(self, pos, value):
        if pos >= self.lenarray:
            self.extend(ExtendibleArray.DefValue(
                pos - self.lenarray + 1000, self.defval))
            self.lenarray += pos - self.lenarray + 1000
        if pos >= self.reallen:
            self.reallen = pos + 1
        self[pos] = value
        self.setPositions.flagSite(pos)

    def __len__(self):
        return self.reallen

    def is_set(self, pos):
        return self.setPositions[pos]

# def __getitem__(self,pos):
#		if pos>=self.lenarray: return 0
#		return array.array.__getitem__(self,pos)
