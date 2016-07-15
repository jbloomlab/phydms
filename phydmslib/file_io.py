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
import Bio.Phylo
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
    for modname in ['Bio', 'cython', 'dms_tools', 'scipy', 'matplotlib']:
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
    for (i, (head1, seq1)) in enumerate(seqs):
        for (head2, seq2) in seqs[i + 1 : ]:
            assert len(seq1) == len(seq2)
            ndiffs = sum(map(lambda (x, y): 1 if (x != y and x != '-' and y != '-') else 0, zip(seq1, seq2)))
#            if ndiffs < 1:
#                raise ValueError("The alignment in {0} has duplicate sequences:\n{1}\n{2}\nPlease remove one of these so that all sequences are unique at non-gap positions.".format(fastafile, head1, head2))
    return seqs


def ReadStringencyBySite(infile):
    """Reads ``*_stringencybsite.txt`` files created by ``phydms``.

    *infile* can be either a readable file-like boject (assumed to be
    at its start) or a string giving the name of an existing file.

    The returned value is a dictionary *stringencybysite* keyed by string
    site number with element *stringencybysite[site]* being a dictionary keyed
    by the three string keys *stringency*, *P*, and *dLnL* each
    specifying the value for that property.

    >>> f = cStringIO.StringIO()
    >>> f.write("# site stringency_ratio P dLnL\\n1 0.017 2.45e-05 8.9\\n2 0.079 0.0005 6.1\\n")
    >>> f.seek(0)
    >>> stringencybysite = ReadStringencyBySite(f)
    >>> set(stringencybysite.keys()) == set(['1', '2'])
    True
    >>> '0.017' == '%.3f' % stringencybysite['1']['stringency']
    True
    >>> '0.0000245' == '%.7f' % stringencybysite['1']['P']
    True
    >>> '8.9' == '%.1f' % stringencybysite['1']['dLnL']
    True
    >>> '0.079' == '%.3f' % stringencybysite['2']['stringency']
    True
    >>> '0.0005' == '%.4f' % stringencybysite['2']['P']
    True
    >>> '6.1' == '%.1f' % stringencybysite['2']['dLnL']
    True
    """
    if isinstance(infile, str):
        with open(infile) as f:
            lines = f.readlines()
    else:
        lines = infile.readlines()
    stringencybysite = {}
    for line in lines:
        if line[0] == '#' or line.isspace() or not line:
            continue
        entries = line.split()
        assert len(entries) == 4, "Unexpected number of entries in line:\n%s" % line
        (site, stringency, P, dLnL) = entries
        assert site not in stringencybysite, "Duplicate site %d" % site
        stringency = float(stringency)
        assert 0 <= stringency, "Invalid stringency: %g" % stringency
        P = float(P)
        assert 0 <= P <= 1, "Invalid P: %g" % P
        dLnL = float(dLnL)
        stringencybysite[site] = {'stringency':stringency, 'P':P, 'dLnL':dLnL}
    return stringencybysite


def ReadOmegaBySite(infile):
    """Reads ``*_omegabysite.txt`` files created by ``phydms``.
   
    *infile* can be either a readable file-like object (assumed to be
    at its start) or a string giving the name of an existing file.

    The returned value is a dictionary *omegabysite* keyed by string
    site number with element *omegabysite[site]* being a dictionary keyed
    by the three string keys *omega*, *P*, and *dLnL* each specifying
    the value for that property.
    
    >>> f = cStringIO.StringIO()
    >>> f.write("# site omega P dLnL\\n1 2.73 0.0012 4.6\\n2 0.03 0.09 1.2\\n")
    >>> f.seek(0)
    >>> omegabysite = ReadOmegaBySite(f)
    >>> set(omegabysite.keys()) == set(['1', '2'])
    True
    >>> '2.73' == '%.2f' % omegabysite['1']['omega']
    True
    >>> '0.0012' == '%.4f' % omegabysite['1']['P']
    True
    >>> '4.6' == '%.1f' % omegabysite['1']['dLnL']
    True
    >>> '0.03' == '%.2f' % omegabysite['2']['omega']
    True
    >>> '0.09' == '%.2f' % omegabysite['2']['P']
    True
    >>> '1.2' == '%.1f' % omegabysite['2']['dLnL']
    True
    """
    if isinstance(infile, str):
        with open(infile) as f:
            lines = f.readlines()
    else:
        lines = infile.readlines()
    omegabysite = {}
    for line in lines:
        if line[0] == '#' or line.isspace() or not line:
            continue
        entries = line.split()
        assert len(entries) == 4, "Unexpected number of entries in line:\n%s" % line
        (site, omega, P, dLnL) = entries
        assert site not in omegabysite, "Duplicate site %d" % site
        omega = float(omega)
        assert 0 <= omega, "Invalid omega: %g" % omega
        P = float(P)
        assert 0 <= P <= 1, "Invalid P: %g" % P
        dLnL = float(dLnL)
        omegabysite[site] = {'omega':omega, 'P':P, 'dLnL':dLnL}
    return omegabysite


def SetMinBrLen(intree, outtree, minbrlen):
    """Sets all branch lengths to be >= a value.

    Reads tree from *intree* newick file, adjusts all
    branch lengths to be >= *minbrlen*, then writes to
    *outtree*.
    """
    assert minbrlen >= 0, "minbrlen must be >= 0"
    tree = Bio.Phylo.read(intree, 'newick')
    for node in tree.get_terminals() + tree.get_nonterminals():
        if node != tree.root:
            node.branch_length = max(minbrlen, node.branch_length)
    with open(outtree, 'w') as f:
        f.write(tree.format('newick').strip())


if __name__ == '__main__':
    import doctest
    doctest.testmod()
