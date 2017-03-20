"""
Module for input / output from files.
"""


import sys
import os
import re
import time
import platform
import importlib
import math
import random
import io
import Bio.Seq
import Bio.SeqIO
import pandas
import phydmslib
import phydmslib.constants


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
    for modname in ['Bio', 'cython', 'numpy', 'scipy', 'matplotlib',
            'natsort', 'sympy', 'six', 'pandas', 'pyvolve', 'statsmodels',
            'weblogolib', 'PyPDF2']:
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

    Read aligned sequences in this example:

    >>> seqs = [('seq1', 'ATGGAA'), ('seq2', 'ATGAAA')]
    >>> f = io.StringIO()
    >>> n = f.write(u'\\n'.join(['>{0}\\n{1}'.format(*tup) for tup in seqs]))
    >>> n = f.seek(0)
    >>> a = ReadCodonAlignment(f, True)
    >>> seqs == a
    True

    Trim stop codons from all sequences in this example:

    >>> seqs = [('seq1', 'ATGTAA'), ('seq2', 'ATGTGA')]
    >>> f = io.StringIO()
    >>> n = f.write(u'\\n'.join(['>{0}\\n{1}'.format(*tup) for tup in seqs]))
    >>> n = f.seek(0)
    >>> a = ReadCodonAlignment(f, True)
    >>> [(head, seq[ : -3]) for (head, seq) in seqs] == a
    True

    Read sequences with gap:

    >>> seqs = [('seq1', 'ATG---'), ('seq2', 'ATGAGA')]
    >>> f = io.StringIO()
    >>> n = f.write(u'\\n'.join(['>{0}\\n{1}'.format(*tup) for tup in seqs]))
    >>> n = f.seek(0)
    >>> a = ReadCodonAlignment(f, True)
    >>> [(head, seq) for (head, seq) in seqs] == a
    True

    Premature stop codon gives error:

    >>> seqs = [('seq1', 'TGAATG'), ('seq2', 'ATGAGA')]
    >>> f = io.StringIO()
    >>> n = f.write(u'\\n'.join(['>{0}\\n{1}'.format(*tup) for tup in seqs]))
    >>> n = f.seek(0)
    >>> a = ReadCodonAlignment(f, True) # doctest: +IGNORE_EXCEPTION_DETAIL
    Traceback (most recent call last):
    ValueError:
    """
    codonmatch = re.compile('^[ATCG]{3}$')
    gapmatch = re.compile('^-+^')
    seqs = [(seq.description.strip(), str(seq.seq).upper()) for seq
            in Bio.SeqIO.parse(fastafile, 'fasta')]
    assert seqs, "{0} failed to specify any sequences".format(fastafile)

    seqlen = len(seqs[0][1])
    if not all([len(seq) == seqlen for (head, seq) in seqs]):
        raise ValueError(("All sequences in {0} are not of the same length; "
                "they must not be properly aligned").format(fastafile))
    if (seqlen < 3) or (seqlen % 3 != 0):
        raise ValueError(("The length of the sequences in {0} is {1} which "
                "is not divisible by 3; they are not valid codon sequences"
                ).format(fastafile, seqlen))

    terminalcodon = []
    codons_by_position = dict([(icodon, []) for icodon in range(seqlen // 3)])
    for (head, seq) in seqs:
        assert len(seq) % 3 == 0
        for icodon in range(seqlen // 3):
            codon = seq[3 * icodon : 3 * icodon + 3]
            codons_by_position[icodon].append(codon)
            if codonmatch.search(codon):
                aa = str(Bio.Seq.Seq(codon).translate())
                if aa == '*':
                    if icodon + 1 != len(seq) // 3:
                        raise ValueError(("In {0}, sequence {1}, non-terminal "
                                "codon {2} is stop codon: {3}").format(
                                fastafile, head, icodon + 1, codon))
            elif codon == '---':
                aa = '-'
            else:
                raise ValueError(("In {0}, sequence {1}, codon {2} is invalid: "
                        "{3}").format(fastafile, head, icodon + 1, codon))
        terminalcodon.append(aa)

    for (icodon, codonlist) in codons_by_position.items():
        if all([codon == '---' for codon in codonlist]):
            raise ValueError(("In {0}, all codons are gaps at position {1}"
                    ).format(fastafile, icodon + 1))

    if all([aa in ['*', '-'] for aa in terminalcodon]):
        if len(seq) == 3:
            raise ValueError(("The only codon is a terminal stop codon for "
                    "the sequences in {0}").format(fastafile))
        seqs = [(head, seq[ : -3]) for (head, seq) in seqs]
    elif any([aa == '*' for aa in terminalcodon]):
        raise ValueError(("Only some sequences in {0} have a terminal stop "
                "codon. All or none must have terminal stop.").format(fastafile))

    if any([gapmatch.search(seq) for (head, seq) in seqs]):
        raise ValueError(("In {0}, at least one sequence is entirely composed "
                "of gaps.").format(fastafile))

    if checknewickvalid:
        if len(set([head for (head, seq) in seqs])) != len(seqs):
            raise ValueError("Headers in {0} not all unique".format(fastafile))
        disallowedheader = re.compile('[\s\:\;\(\)\[\]\,\'\"]')
        for (head, seq) in seqs:
            if disallowedheader.search(head):
                raise ValueError(("Invalid character in header in {0}:"
                        "\n{2}").format(fastafile, head))

    return seqs


def readPrefs(prefsfile, minpref=0, avgprefs=False, randprefs=False,
        seed=1, sites_as_strings=False):
    """Read preferences from file with some error checking.

    Args:
        `prefsfile` (string or readable file-like object)
            File holding amino-acid preferences. Can be
            comma-, space-, or tab-separated file with column
            headers of `site` and then all one-letter amino-acid
            codes, or can be in the more complex format written
            `dms_tools`. Must be prefs for consecutively numbered
            sites starting at 1. Stop codon prefs can be present
            (stop codons are indicated by ``*``); if so they are
            removed and prefs re-normalized to sum to 1.
        `minpref` (float >= 0)
            Adjust all preferences to be >= this number.
        `avgprefs`, `randprefs` (bool)
            Mutually exclusive options specifying to average or
            randomize prefs across sites.
        `seed` (int)
            Seed used to sort random number generator for `randprefs`.
        `sites_as_strings` (bool)
            By default, the site numers are coerced to integers.
            If this option is `True`, then they are kept as strings.

    Returns:
        `prefs` (dict)
            `prefs[r][a]` is the preference of site `r` for amino-acid `a`.
            `r` is an `int` unless `sites_as_strings=True`.
    """
    assert minpref >= 0, 'minpref must be >= 0'

    aas = set(phydmslib.constants.AA_TO_INDEX.keys())

    try:
        df = pandas.read_csv(prefsfile, sep=None, engine='python')
        pandasformat = True
    except ValueError:
        pandasformat = False
    if pandasformat and (set(df.columns) == aas.union(set(['site'])) or
            set(df.columns) == aas.union(set(['site', '*']))):
        # read valid preferences as data frame
        sites = df['site'].tolist()
        prefs = {}
        for r in sites:
            rdf = df[df['site'] == r]
            prefs[r] = {}
            for aa in df.columns:
                if aa != 'site':
                    prefs[r][aa] = float(rdf[aa])
    else:
        # try reading as dms_tools format
        prefs = phydmslib.file_io.readPrefs_dms_tools_format(prefsfile)[2]
        sites = list(prefs.keys())

    # error check prefs
    if not sites_as_strings:
        try:
            sites = [int(r) for r in sites]
        except ValueError:
            raise ValueError("sites not int in prefsfile {0}".format(prefsfile))
        assert (min(sites) == 1 and max(sites) - min(sites) == len(sites) - 1),\
                "Sites not consecutive starting at 1"
        prefs = dict([(int(r), rprefs) for (r, rprefs) in prefs.items()])
    else:
        sites = [str(r) for r in sites]
        prefs = dict([(str(r), rprefs) for (r, rprefs) in prefs.items()])

    assert len(set(sites)) == len(sites), "Non-unique sites in prefsfiles"
    assert all([all([pi >= 0 for pi in rprefs.values()]) for rprefs in
            prefs.values()]), "prefs < 0 in prefsfile {0}".format(prefsfile)
    for r in list(prefs.keys()):
        rprefs = prefs[r]
        assert sum(rprefs.values()) - 1 <= 0.01, (
            "Prefs in prefsfile {0} don't sum to one".format(prefsfile))
        if '*' in rprefs:
            del rprefs['*']
        assert aas == set(rprefs.keys()), ("prefsfile {0} does not include "
                "all amino acids at site {1}").format(prefsfile, r)
        rsum = float(sum(rprefs.values()))
        prefs[r] = dict([(aa, pi / rsum) for (aa, pi) in rprefs.items()])
    assert set(sites) == set(prefs.keys())

    # Iteratively adjust until all prefs exceed minpref after re-scaling.
    for r in list(prefs.keys()):
        rprefs = prefs[r]
        iterations = 0
        while any([pi < minpref for pi in rprefs.values()]):
            rprefs = dict([(aa, max(1.1 * minpref,
                    pi)) for (aa, pi) in rprefs.items()])
            newsum = float(sum(rprefs.values()))
            rprefs = dict([(aa, pi / newsum) for (aa, pi) in rprefs.items()])
            iterations += 1
            assert iterations <= 3, "minpref adjustment not converging."
        prefs[r] = rprefs

    if randprefs:
        assert not avgprefs, "randprefs and avgprefs are incompatible"
        random.seed(seed)
        sites = sorted([r for r in prefs.keys()])
        prefs = [prefs[r] for r in sites]
        random.shuffle(sites)
        prefs = dict(zip(sites, prefs))
    elif avgprefs:
        avg_prefs = dict([(aa, 0.0) for aa in aas])
        for rprefs in prefs.values():
            for aa in aas:
                avg_prefs[aa] += rprefs[aa]
        for aa in aas:
            avg_prefs[aa] /= float(len(prefs))
        for r in list(prefs.keys()):
            prefs[r] = avg_prefs

    return prefs


def readPrefs_dms_tools_format(f):
    """Reads the amino-acid preferences written by `dms_tools`.

    This is an exact copy of the same code from
    `dms_tools.file_io.ReadPreferences`. It is copied because
    `dms_tools` is currently only compatible with `python2`, and
    we needed something that also works with `python3`.

    *f* is the name of an existing file or a readable file-like object.
    It should be in the format written by `dms_tools`.

    The return value is the tuple: *(sites, wts, pi_means, pi_95credint, h)*
    where *sites*, *wts*, *pi_means*, and *pi_95credint* will all
    have the same values used to write the file with *WritePreferences*,
    and *h* is a dictionary with *h[r]* giving the site entropy (log base
    2) for each *r* in *sites*.
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

def readDivPressure(fileName):
    """Reads in diversifying pressures from some file.

    Scale diversifying pressure values so absolute value of the max value is 1,
    unless all values are zero.

    Args:
        `fileName` (string or readable file-like object)
            File holding diversifying pressure values. Can be
            comma-, space-, or tab-separated file. The first column
            is the site (consecutively numbered, sites starting
            with one) and the second column is the diversifying pressure values.

    Returns:
        `divPressure` (dict keyed by ints)
            `divPressure[r][v]` is the diversifying pressure value of site `r`.
    """
    try:
        df = pandas.read_csv(fileName, sep=None, engine='python')
        pandasformat = True
    except ValueError:
        pandasformat = False
    df.columns = ['site', 'divPressureValue']
    scaleFactor = max(df["divPressureValue"].abs())
    if scaleFactor > 0:
        df["divPressureValue"] = [x / scaleFactor for x in df["divPressureValue"]]
    assert len(df['site'].tolist()) == len(set(df['site'].tolist())),"There is at least one non-unique site in {0}".format(fileName)
    assert max(df["divPressureValue"].abs()) <= 1, "The scaling produced a diversifying pressure value with an absolute value greater than one."
    sites = df['site'].tolist()
    divPressure = {}
    for r in sites:
        divPressure[r] = df[df['site'] == r]["divPressureValue"].tolist()[0]
    return divPressure


if __name__ == '__main__':
    import doctest
    doctest.testmod()
