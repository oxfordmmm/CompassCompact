import array


class ContigMap:
    def __init__(self):
        self.gmap = array.array("L")
        self.sizearray = 0

    def flagSite(self, pos):
        posinElement = pos & 31
        element = pos >> 5
        if pos >= self.sizearray:
            extraNbits = pos - (self.sizearray)
            extraNbits = extraNbits - posinElement + 32
            self.gmap.extend([0] * (extraNbits >> 5))
            self.sizearray += 32 * (extraNbits >> 5)

        self.gmap[element] |= 1 << posinElement

    def __getitem__(self, pos):
        if pos >= self.sizearray:
            return False
        posinElement = pos & 31
        element = pos >> 5
        return self.gmap[element] & (1 << posinElement) != 0

    def __iter__(self):
        k = 0
        for i in self.gmap:
            j = 0
            while j < 32:
                if i & (1 << j):
                    yield k
                k += 1
                j += 1


class GenomeMap:
    '''		The categories are defined as:
            nilspanning                    allsites - inDP
            yesspanning_nilmapped            DPsites_with_DP=0
            yesspanning_yesmapped_yeshiqual   DP4sites_with_DP4
            yesyesspanning_yesmapped_nilhiqual   yesspanning_yesmapped - yesspanning_yesmapped_yeshiqual
    '''

    def __init__(self):
        self.contigs = {}

    def addSite(self, contig, pos):
        if contig not in self.contigs:
            self.contigs.setdefault(contig, ContigMap())
        self.contigs[contig].flagSite(pos)

    def getSite(self, contig, pos):
        if contig not in self.contigs:
            return False
        return self.contigs[contig][pos]

    def __getitem__(self, k):
        return self.contigs[k]
