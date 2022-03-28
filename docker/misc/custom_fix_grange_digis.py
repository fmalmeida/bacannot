import copy
import logging

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord


class Grange:
    def __init__(self, genome_name, chrom, start, end, strand, genome_seq, genome_len, circular=True):
        if start <= 0:
            start = 1
        if end <= 0:
            end = 1
        if start == 0 and end > start:
            start += 1
            if end < genome_len:
                end += 1
        if start <= 0 or end <= 0 or start > genome_len or end > genome_len: 
            raise ValueError("GRange: start or end position is out of the range.")
        self.genome_name = genome_name
        self.chr = chrom
        self.start = start
        self.end = end
        self.strand = strand
        self.genome_seq = genome_seq
        self.genome_len = genome_len
        self.circular = circular
        self.width = self.__len__()

    def set_start(self, start):
        self.start = start
        self.width = self.__len__()

    def set_end(self, end):
        self.end = end
        self.width = self.__len__()

    def get_flank_range(self, flank):
        if flank + flank + self.width > self.genome_len:
            max_flank = int((self.genome_len - self.width) / 2)
        else:
            max_flank = flank
        flank_start = self.start - max_flank
        if flank_start <= 0:
            flank_start = flank_start + self.genome_len if self.circular else 1
        flank_end = self.end + max_flank
        if flank_end > self.genome_len:
            flank_end = flank_end - self.genome_len if self.circular else self.genome_len

        return Grange(self.genome_name, self.chr, flank_start, flank_end, self.strand,
                      self.genome_seq, self.genome_len, self.circular)

    def shift_left(self, size):
        self.start -= size
        if self.start <= 0:
            self.start = self.start + self.genome_len if self.circular else 1
        self.end -= size
        if self.end <= 0:
            self.end = self.end + self.genome_len if self.circular else 1
        self.width = self.__len__()

    def remap_offsets(self, left_offset, right_offset):
        if self.strand == "+":
            self.end = self.start + right_offset - 1
            self.start = self.start + left_offset - 1
            if self.start > self.genome_len:
                self.start = self.start - self.genome_len if self.circular else self.genome_len
            if self.end > self.genome_len:
                self.end = self.end - self.genome_len if self.circular else self.genome_len
        else:
            self.start = self.end - right_offset
            self.end = self.end - left_offset
            if self.start <= 0:
                self.start = self.start + self.genome_len if self.circular else 1
            if self.end <= 0:
                self.end = self.end + self.genome_len if self.circular else 1
        self.width = self.__len__()

    def has_overlap(self, other, ignore_strand=False, flank=0):
        return self.get_overlap_length(other, ignore_strand, flank) > 0

    def get_overlap_length(self, other, ignore_strand=False, flank=0):
        other_range = copy.copy(other)
        new_range = self.get_flank_range(flank)
        size = min(new_range.start, other_range.start)-1
        new_range.shift_left(size)
        other_range.shift_left(size)
        if new_range.start <= other_range.start:
            overlap = min(new_range.end-other_range.start+1, other.width)
        else:
            overlap = min(other_range.end-new_range.start+1, new_range.width)

        if self.strand != other.strand and not ignore_strand:
            overlap = 0

        # print(self.start, self.end, self.strand, other.start, other.end, other.strand, flank, overlap)
        return overlap

    def is_inside(self, other, ignore_strand=False):
        return self.get_overlap_length(other, ignore_strand, 0) == self.width

    def get_flank_lengths(self, flank):
        new_range = self.get_flank_range(flank)
        left_flank = self.start - new_range.start
        if left_flank <= 0:
            left_flank += self.genome_len
        right_flank = new_range.end - self.end
        if right_flank <= 0:
            right_flank += self.genome_len
        return left_flank, right_flank

    def get_sequence(self, flank=0, protein=False):

        new_range = self.get_flank_range(flank)

        if new_range.start <= new_range.end:
            seq = self.genome_seq[new_range.start-1:new_range.end]
        else:  # element crossing the genome boundary
            seq = self.genome_seq[new_range.start-1:self.genome_len] + self.genome_seq[0:new_range.end]

        if self.strand == '-':
            seq = seq.reverse_complement()
        if protein:
            seq = seq.translate(table=11)

        return SeqRecord(seq, id=self.genome_name, description='')

    def __str__(self):
        return "{} {} {} {} {}".format(self.genome_name, self.chr, self.start, self.end, self.strand)

    def __len__(self):
        out_len = self.end - self.start + 1
        if out_len <= 0:  # element crossing the genome boundary
            out_len += self.genome_len
        return out_len