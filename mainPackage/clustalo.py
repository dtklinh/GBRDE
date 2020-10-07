"""
Simple interface to Clustal Omega.
"""
from csb.bio.sequence import SequenceTypes, SequenceCollection
import numpy as np


try:
    import Bio
except:
    raise Exception('requires Biophython')

import tempfile, os


from csb.core import OrderedDict
from csb.bio.sequence import Sequence, AbstractSequence
from csb.bio.io import SequenceParser

KYTE_DOOLITTLE = \
{"I":4.5,
 "V":4.2,
 "L":3.8,
 "F":2.8,
 "C":2.5,
 "M":1.9,
 "A":1.8,
 "G":-0.4,
 "T":-0.7,
 "W":-0.9,
 "S":-0.8,
 "Y":-1.3,
 "P":-1.6,
 "H":-3.2,
 "E":-3.5,
 "Q":-3.5,
 "D":-3.5,
 "N":-3.5,
 "K":-3.9,
 "R":-4.5}

class AminoAcidStats(object):

    alphabet = 'ACDEFGHIKLMNPQRSTVWY'

    @classmethod
    def composition(cls, author='beals'):
        """
        Amino acid frequencies observed in the sequence data-bases.

        The index corresponds to different publications:

        M. Beals, L. Gross, S. Harrell
        www.tiem.utk.edu/~harrell/webmodules/aminoacid.htm

        Dayhoff (1978)

        Jones, D.T. Taylor, W.R. & Thornton, J.M.(1991)
        CABIOS 8:275-282
        """
        
        if author.lower() == 'beals':
            freqs = {"A": 0.074, "R": 0.042, "C": 0.033, "G": 0.074,
                     "H": 0.029, "N": 0.044, "D": 0.059, "E": 0.058,
                     "Q": 0.037, "I": 0.038, "L": 0.076, "K": 0.072,
                     "M": 0.018, "F": 0.040, "P": 0.050, "S": 0.081,
                     "T": 0.062, "W": 0.013, "Y": 0.033, "V": 0.068}

        elif author.lower() == 'dayhoff':
            freqs = {"L": 0.085, "A": 0.087, "G": 0.089, "S": 0.070,
                     "V": 0.065, "E": 0.050, "T": 0.058, "K": 0.081,
                     "I": 0.037, "D": 0.047, "R": 0.041, "P": 0.051,
                     "N": 0.040, "Q": 0.038, "F": 0.040, "Y": 0.030,
                     "M": 0.015, "H": 0.034, "C": 0.033, "W": 0.010}

        elif author.lower() == 'jones':
            freqs = {"L": 0.091, "A": 0.077, "G": 0.074, "S": 0.069,
                     "V": 0.066, "E": 0.062, "T": 0.059, "K": 0.059,
                     "I": 0.053, "D": 0.052, "R": 0.051, "P": 0.051,
                     "N": 0.043, "Q": 0.041, "F": 0.040, "Y": 0.032,
                     "M": 0.024, "H": 0.023, "C": 0.020, "W": 0.014}

        else:
            msg = 'Author {} unknown'
            raise ValueError(msg.format(author))
        
        return freqs

    def composition_vector(cls, author='beals'):
        """
        Returns a numeric array of data-base frequencies where the
        indices of the vector correspond to the position in the AA-
        alphabet.
        In order to control the db-type etc. use the index (see doc-
        umentation of dbFrequencies).
        """
        
        freqs = self.composition(author)
        freqs = np.array(map(freqs.__getitem__, cls.alphabet))

        return freqs / freqs.sum()
    
class Alignment(object):
    """
    A simple alignment object that behaves like an ordered dictionary
    and has a few more alignment specific functions.
    """
    gap = '-'

    def __init__(self, seqs=None, enforce_length=False):

        self._seqs = OrderedDict()
        self._length = 0
        self._enforce_length = bool(enforce_length)
        self._filename = None
        
        if seqs is not None:
            for seq in seqs: self.add(seq)

    def from_fasta(self, filename):
        
        parser = SequenceParser()
        for seq in parser.read(filename):
            self.add(seq)

        self._filename = filename
        
    def __iter__(self):
        return iter(self._seqs.values())

    def __str__(self):
        return self.to_fasta()
            
    @property
    def size(self):
        return len(self._seqs)
    @property
    def length(self):
        return self._length
    @property
    def ids(self):
        return self._seqs.keys()

    def matching_columns(self, *ids):
        """
        Returns a list of column indices where none of the specified
        sequences has a gap.
        """
        seqs = [self[i] for i in ids]
        length = min([seq.length for seq in seqs])
        seqs = [seq.sequence for seq in seqs]

        matches = []

        for k in range(length):

            col = [seq[k] for seq in seqs]
            if Alignment.gap in col: continue

            matches.append(k)

        return matches

    def column2number(self, a):
        """
        Maps alignment columns to sequence number for non-gap
        residues.
        """
        seq = self[a]
        indexmap = {}

        j = 0

        for i in range(len(seq)):
            if seq.sequence[i] == Alignment.gap: continue
            indexmap[i] = j
            j +=1

        return indexmap
            
    def matching_numbers(self, *ids):
        """
        Returns of list of tuples of the matching columns converted
        to the respective sequence numbers.
        """
        matches = self.matching_columns(*ids)
        maps = [self.column2number(i) for i in ids]
        numbers = [map(m.__getitem__,matches) for m in maps]
        
        return zip(*numbers)
    
    def add(self, seq, name=None):
        """
        Add a sequence
        """
        if type(seq) is str:
            if name is None:
                name = 'seq<{}>'.format(self.size)
            seq = Sequence(name, name, seq)
            
        elif not isinstance(seq, AbstractSequence):
            msg = 'Sequence must be provided either as instance of CSB Sequence or string'
            raise TypeError(msg)

        if seq.id in self._seqs :
            msg = 'Sequence <{}> already in alignment'
            raise KeyError(msg.format(seq.id))
        
        self._seqs[seq.id] = seq
        self._length = max(self._length, seq.length)
        
    def _getbyname(self, name):

        if not name in self._seqs:
            msg = 'Unknown sequence: {}'
            raise KeyError(msg.format(name))

        return self._seqs[name]

    def _getbyindex(self, index):

        if index > self.size:
            msg = 'Only {} sequences contained in alignment'
            raise IndexError(msg.format(index))

        return self._seqs.values()[index]
        
    def __getitem__(self, name_or_index):

        if type(name_or_index) is int:
            return self._getbyindex(name_or_index)
        elif type(name_or_index) is str:
            return self._getbyname(name_or_index)
        else:
            msg = 'Sequences are retrieved either by name or index'
            return KeyError(msg)

    @staticmethod
    def _insert_breaks(seq, length=60):
        s = [seq[i*length:(i+1)*length] for i in range(len(seq)/length+1)]
        return '\n'.join(s)

    def to_fasta(self, filename=None):

        out = ''
        for seq in self:
            out += '>{0}\n{1}\n'.format(seq.id, Alignment._insert_breaks(seq.sequence))

        if filename is None:
            return out

        with open(filename,'w') as f:
            f.write(out)
        
if __name__ == '__main__':

    aln = Alignment()
    aln.add('asdfasdf','a')
    aln.add('ALKSDSLKF', 'b')
    




class ClustalO(object):

    cmd = 'clustalo -i {infile} -o {outfile}'
    prefix = 'seq_{0}'

    @classmethod
    def _run(cls, seqfile):

        from Bio.Align.Applications import ClustalOmegaCommandline

        cmdline = ClustalOmegaCommandline(cls.cmd.split()[0], infile=seqfile, outfile='tmp_alignment')
        cmdline()
        outfile = 'tmp_alignment'
        if not os.path.exists(outfile):
            raise FunctionError('clustalo failed; {} missing'.format(outfile))

        return outfile

    @classmethod
    def run(cls, seqfile, cleanup=True):

        outfile1 = cls._run(seqfile)

        parser = SequenceParser()
        ali = parser.parse_file('tmp_alignment')

        if cleanup: os.unlink('tmp_alignment')

        return ali

    def align(self, sequences, cleanup=True):

        counter = 0
        seqs = []

        for sequence in sequences:

            if type(sequence) == str:
                sequence = Sequence(self.prefix.format(counter),
                                    self.prefix.format(counter),
                                    sequence,
                                    type=SequenceTypes.Protein)
                counter += 1
            elif not isinstance(sequence, Sequence):
                raise TypeError(str(sequence))

            seqs.append(sequence)

        seqfile = tempfile.mktemp() + '.fasta'
        self._tempfile = seqfile

        seqs = SequenceCollection(seqs)
        seqs.to_fasta(seqfile)

        A = Alignment()
        for seq in self.__class__.run(seqfile, cleanup):
            iden = seq.id.replace('P1;','')
            sequence = seq.sequence.replace('*','')
            A.add(Sequence(iden,iden,sequence))

        if cleanup:
            os.unlink(seqfile)

        return A


def align(sequences):
    clustalo = ClustalO()
    return clustalo.align(sequences)


if __name__ == '__main__':

    a = 'ARTTNMRAR'
    b = 'ARSTNVKA'    

    sequences = [a,b]

    clustalo = ClustalO()
    A = clustalo.align(sequences, True)
