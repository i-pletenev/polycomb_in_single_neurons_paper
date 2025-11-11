
chroms = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14',
       'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chrX'] + ['trans']

class Chromosome:
    chromosomes = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 
                   'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14',
                   'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21',
                   'chrX', 'trans']
    
    @classmethod
    def get_chr_order(cls):
        return {cls.chromosomes[i]: i for i in range(len(cls.chromosomes))}
        
    def __init__(self, chr):
        if isinstance(chr, int):
            chr = 'chr' + str(chr)
        if isinstance(chr, str) and chr.startswith('chrom'):
            chr = 'chr' + chr.replace('chrom', '')
        self.chr = chr

    def __hash__(self):
        return hash(self.chr)

    def __eq__(self, other):
        if isinstance(other, str):
            other = Chromosome(other)
        order = self.get_chr_order()
        return order[self.chr] == order[other.chr]

    def __gt__(self, other):
        if isinstance(other, str):
            other = Chromosome(other)
        order = self.get_chr_order()
        return order[self.chr] > order[other.chr]

    def __ge__(self, other):
        if isinstance(other, str):
            other = Chromosome(other)
        order = self.get_chr_order()
        return order[self.chr] >= order[other.chr]

    def __str__(self):
        return str(self.chr)

    def __repr__(self):
        return str(self.chr)
        
class Loop:
    fields = ('chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2')
    
    def __init__(self, *args):
        if len(args) == 6:
            chrom1, start1, end1, chrom2, start2, end2 = args
        elif len(args) == 1:
            args = args[0]
            if len(args) == 2:
                chrom1, start1, end1 = args[0]
                chrom2, start2, end2 = args[1]
            else:
                raise ValueError(args)
        else:
            raise ValueError(args)
    
        chrom1 = Chromosome(chrom1)
        start1 = int(start1)
        end1 = int(end1)
        chrom2 = Chromosome(chrom2)
        start2 = int(start2)
        end2 = int(end2)

        # make left anchor be always lower
        if (chrom2 < chrom1) or (chrom1 == chrom2 and start2 < start1) or (chrom1 == chrom2 and start1 == start2 and end2 < end1):
            chrom1, start1, end1, chrom2, start2, end2 = chrom2, start2, end2, chrom1, start1, end1

        self.chrom1 = chrom1
        self.start1 = start1
        self.end1 = end1
        self.chrom2 = chrom2
        self.start2 = start2
        self.end2 = end2        

    def __getitem__(self, key):
        return getattr(self, key)

    def __hash__(self):
        return sum(hash(str(getattr(self, field))) for field in self.fields)

    def __eq__(self, other):
        if (not isinstance(self, Loop)) or (not isinstance(other, Loop)):
            return str(self) == str(other)
        #return any(getattr(self, filed) != getattr(other, filed) for filed in self.fields)
        return hash(self) == hash(other)

    def __le__(self, other):
        if (not isinstance(self, Loop)) or (not isinstance(other, Loop)):
            return str(self) <= str(other)
        for filed in self.fields:
            self_f = getattr(self, filed)
            other_f = getattr(other, filed)
            if self_f < other_f: return True
            elif self_f == other_f: continue
            elif self_f > other_f: return False
        else:
            return True

    def __lt__(self, other):
        if (not isinstance(self, Loop)) or (not isinstance(other, Loop)):
            return str(self) < str(other)
        for filed in self.fields:
            self_f = getattr(self, filed)
            other_f = getattr(other, filed)
            if self_f < other_f: return True
            elif self_f == other_f: continue
            elif self_f > other_f: return False
        else:
            return False

    def __gt__(self, other):
        if (not isinstance(self, Loop)) or (not isinstance(other, Loop)):
            return str(self) > str(other)
        for filed in self.fields:
            self_f = getattr(self, filed)
            other_f = getattr(other, filed)
            if self_f > other_f: return True
            elif self_f == other_f: continue
            elif self_f < other_f: return False
        else:
            return False

    def __ge__(self, other):
        if (not isinstance(self, Loop)) or (not isinstance(other, Loop)):
            return str(self) >= str(other)
        for filed in self.fields:
            self_f = getattr(self, filed)
            other_f = getattr(other, filed)
            if self_f > other_f: return True
            elif self_f == other_f: continue
            elif self_f < other_f: return False
        else:
            return True
        
    def __str___(self):
        return fr"(('{self.chrom1}', {self.start1/1000000}, {self.end1/1000000}),('{self.chrom2}', {self.start2/1000000}, {self.end2/1000000}))"

    def __str__(self):
        return fr"(('{self.chrom1}', {self.start1}, {self.end1}), ('{self.chrom2}', {self.start2}, {self.end2}))"
    
    def __repr__(self):
        return fr"(('{self.chrom1}', {self.start1}, {self.end1}), ('{self.chrom2}', {self.start2}, {self.end2}))"

    def __tuple__(self):
        return ((self.chrom1, self.start1, self.end1),(self.chrom2, self.start2, self.end2))

    def tuple(self):
       return self.__tuple__()

    def chroms(self, n=1):
        if n == 1:
            return self.chrom1
        elif n == 2:
            return self.chrom2
        else:
            return self.chrom1, self.chrom2

    def is_trans(self):
        return self.chrom1 != self.chrom2

    def what_type(self):
        return 'trans' if self.is_trans() else 'cis'

    def get_chrom(self):
        if self.chrom1 == self.chrom2:
            return self.chrom1
        else:
            return 'trans'