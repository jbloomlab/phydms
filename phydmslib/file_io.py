"""
Module for input / output from files.
"""


import sys
import os
import re
import time
import tempfile
import platform
import importlib
import math
import cStringIO
import Bio.SeqIO
import phydmslib



def Versions():
    """Returns a string with version information.

    You would call this function if you want a string giving detailed information
    on the version of ``phydms`` and the associated packages that it uses.
    """
    s = [\
            'Version information:',
            '\tTime and date: %s' % time.asctime(),
            '\tPlatform: %s' % platform.platform(),
            '\tPython version: %s' % sys.version.replace('\n', ' '),
            '\tphydms version: %s' % phydmslib.__version__,
            ]
    for modname in ['Bio', 'cython', 'dms_tools', 'scipy']:
        try:
            v = importlib.import_module(modname).__version__
            s.append('\t%s version: %s' % (modname, v))
        except ImportError:
            s.append('\t%s cannot be imported into Python' % modname)
    return '\n'.join(s)


def ReadCodonAlignment(fastafile, checknewickvalid):
    """Reads codon alignment from file.
    
    *fastafile* is the name of an existing FASTA file.

    *checknewickvalid* : if *True*, we require that names are unique and do
    **not** contain spaces, commas, colons, semicolons, parentheses, square brackets,
    or single or double quotation marks.
    If any of these disallowed characters are present, raises an Exception.
    
    Reads the alignment from the *fastafile* and returns the aligned
    sequences as a list of 2-tuple of strings *(header, sequence)*
    where *sequence* is upper case.

    If the terminal codon is a stop codon for **all** sequences, then
    this terminal codon is trimmed. Raises an exception if the sequences
    are not aligned codon sequences that are free of stop codons (with
    the exception of a shared terminal stop) and free of ambiguous nucleotides.
    """
    dnamatch = re.compile('^[ATCG]{3}$')
    gapmatch = re.compile('^-+^')
    seqs = [(seq.description.strip(), str(seq.seq).upper()) for seq in Bio.SeqIO.parse(fastafile, 'fasta')]
    if not seqs:
        raise ValueError("%s failed to specify any sequences" % fastafile)
    assert isinstance(seqs[0][0], str) and isinstance(seqs[0][1], str)
    if not all([len(seq) == len(seqs[0][1]) for (head, seq) in seqs]):
        raise ValueError("All sequences in %s are not of the same length; they are not be properly aligned" % fastafile)
    if (len(seqs[0][1]) < 3) or (len(seqs[0][1]) % 3 != 0):
        raise ValueError("The length of the sequences in %s is %d which is not divisible by 3; they are not valid codon sequences" % (fastafile, len(seqs[0][1])))
    terminalcodon = []
    codons_by_position = dict([(icodon, []) for icodon in range(len(seq) // 3)])
    for (head, seq) in seqs:
        assert len(seq) % 3 == 0
        for icodon in range(len(seq) // 3):
            codon = seq[3 * icodon : 3 * icodon + 3]
            codons_by_position[icodon].append(codon)
            if dnamatch.search(codon):
                aa = str(Bio.Seq.Seq(codon).translate())
                if aa == '*':
                    if icodon + 1 != len(seq) // 3:
                        raise ValueError("In %s, sequence %s, the non-terminal codon %d is a stop codon: %s" % (fastafile, head, icodon + 1, codon))
            elif codon == '---':
                aa = '-'
            else:
                raise ValueError("In %s, sequence %s, codon %d is invalid: %s" % (fastafile, head, icodon + 1, codon))
        terminalcodon.append(aa)
    for (icodon, codonlist) in codons_by_position.items():
        if all([codon == '---' for codon in codonlist]):
            raise ValueError("In %s, all codons are gaps at position %d" % (fastafile, icodon + 1))
    if all([aa in ['*', '-'] for aa in terminalcodon]):
        if len(seq) == 3:
            raise ValueError("The only codon is a terminal stop codon for the sequences in %s" % fastafile)
        seqs = [(head, seq[ : -3]) for (head, seq) in seqs]
    elif any([aa == '*' for aa in terminalcodon]):
        raise ValueError("Only some sequences in %s have a terminal stop codon. Either all or none must have a terminal stop" % fastafile)
    if any([gapmatch.search(seq) for (head, seq) in seqs]):
        raise ValueError("In %s, at least one sequence is entirely composed of gaps" % fastafile)
    if checknewickvalid:
        if len(set([head for (head, seq) in seqs])) != len(seqs):
            raise ValueError("The headers in %s are not all unique" % fastafile)
        disallowedheader = re.compile('[\s\:\;\(\)\[\]\,\'\"]')
        for (head, seq) in seqs:
            if disallowedheader.search(head):
                raise ValueError("There is an invalid character in the following header in %s:\n%s" % (fastafile, head))
    return seqs

    


if __name__ == '__main__':
    import doctest
    doctest.testmod()
