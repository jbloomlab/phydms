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
import io
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
    for modname in ['Bio', 'scipy', 'matplotlib']:
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
    **not** contain spaces, commas, colons, semicolons, parentheses, square 
    brackets, or single or double quotation marks.
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
            ndiffs = sum([1 if (x_y[0] != x_y[1] and x_y[0] != '-' and x_y[1] != '-') else 0 for x_y in zip(seq1, seq2)])
#            if ndiffs < 1:
#                raise ValueError("The alignment in {0} has duplicate sequences:\n{1}\n{2}\nPlease remove one of these so that all sequences are unique at non-gap positions.".format(fastafile, head1, head2))
    return seqs


def ReadOmegaBySite(infile):
    """Reads ``*_omegabysite.txt`` files created by ``phydms``.
   
    *infile* can be either a readable file-like object (assumed to be
    at its start) or a string giving the name of an existing file.

    The returned value is a dictionary *omegabysite* keyed by string
    site number with element *omegabysite[site]* being a dictionary keyed
    by the three string keys *omega*, *P*, and *dLnL* each specifying
    the value for that property.
    
    >>> f = io.StringIO()
    >>> n = f.write(u"# site omega P dLnL\\n1 2.73 0.0012 4.6\\n2 0.03 0.09 1.2\\n")
    >>> n = f.seek(0)
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


def ReadPreferences(f):
    """Reads the amino-acid preferences written by `dms_tools`.
    
    This is an exact copy of the same code from 
    `dms_tools.file_io.ReadPreferences`. It is copied because
    `dms_tools` is currently only compatible with `python2`, and
    we needed something that also works with `python3`.

    *f* is the name of an existing file or a readable file-like object.

    The return value is the tuple: *(sites, wts, pi_means, pi_95credint, h)*
    where *sites*, *wts*, *pi_means*, and *pi_95credint* will all
    have the same values used to write the file with *WritePreferences*,
    and *h* is a dictionary with *h[r]* giving the site entropy (log base
    2) for each *r* in *sites*.

    See docstring of *WritePreferences* for example usage.
    """
    charmatch = re.compile('^PI_([A-z\*\-]+)$')
    if isinstance(f, str):
        f = open(f)
        lines = f.readlines()
        f.close()
    else:
        lines = f.readlines()
    characters = []
    sites = []
    wts = {}
    pi_means = {}
    pi_95credint = {}
    h = {}
    for line in lines:
        if line.isspace():
            continue
        elif line[0] == '#' and not characters:
            entries = line[1 : ].strip().split()
            if len(entries) < 4:
                raise ValueError("Insufficient entries in header:\n%s" % line)
            if not (entries[0] in ['POSITION', 'SITE'] and entries[1][ : 2] == 'WT' and entries[2] == 'SITE_ENTROPY'):
                raise ValueError("Not the correct first three header columns:\n%s" % line)
            i = 3
            while i < len(entries) and charmatch.search(entries[i]):
                characters.append(charmatch.search(entries[i]).group(1))
                i += 1
            if i  == len(entries):
                pi_95credint = None
                linelength = len(characters) + 3
            else:
                if not len(entries) - i == len(characters):
                    raise ValueError("Header line does not have valid credible interval format:\n%s" % line)
                if not all([entries[i + j] == 'PI_%s_95' % characters[j] for j in range(len(characters))]):
                    raise ValueError("mean and credible interval character mismatch in header:\n%s" % line)
                linelength = 2 * len(characters) + 3
        elif line[0] == '#':
            continue
        elif not characters:
            raise ValueError("Found data lines before encountering a valid header")
        else:
            entries = line.strip().split()
            if len(entries) != linelength:
                raise ValueError("Line does not have expected %d entries:\n%s" % (linelength, line))
            r = entries[0]
            assert r not in sites, "Duplicate site of %s" % r
            sites.append(r)
            wts[r] = entries[1]
            assert entries[1] in characters or entries[1] == '?', "Character %s is not one of the valid ones in header. Valid possibilities: %s" % (entries[1], ', '.join(characters))
            h[r] = float(entries[2])
            pi_means[r] = dict([(x, float(entries[3 + i])) for (i, x) in enumerate(characters)])
            if pi_95credint != None:
                pi_95credint[r] = dict([(x, (float(entries[3 + len(characters) + i].split(',')[0]), float(entries[3 + len(characters) + i].split(',')[1]))) for (i, x) in enumerate(characters)])
    return (sites, wts, pi_means, pi_95credint, h)


if __name__ == '__main__':
    import doctest
    doctest.testmod()
