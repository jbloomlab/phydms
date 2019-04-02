"""Functions for performing model adequacy tests.

Written by Sarah Hilton.
"""

from phydmslib.constants import *
import phydmslib.models
import pandas as pd
import scipy
import math
import multiprocessing
from statsmodels.sandbox.stats.multicomp import multipletests


# global variable
amino_acids = list(INDEX_TO_AA.values())
amino_acids.sort()


def calc_aa_frequencies(alignment):
    """Calculate amino-acid frequencies from a codon alignment.

    Args:
        `alignment` (list)
            Alignment of codon sequences as a list of tuples, (seq_id, seq)
    Returns:
        `pandas` dataframe of amino-acid frequencies by site

    >>> answer = pd.DataFrame({"site": [1, 2], "A": [0.0, 0.0],\
                               "C": [0.0, 0.0], "D": [0.0, 0.0],\
                               "E": [0.0, 0.0], "F": [0.0, 0.0],\
                               "G": [0.0, 0.0], "H": [0.0, 0.0],\
                               "I": [0.0, 0.0], "K": [0.0, 0.0],\
                               "L": [0.0, 0.0], "M": [0.0, 0.0],\
                               "N": [0.0, 0.0], "P": [0.0, 0.0],\
                               "Q": [0.0, 0.0], "R": [0.0, 0.0],\
                               "S": [0.0, 0.0], "T": [0.0, 0.0],\
                               "V": [0.0, 0.0], "W": [0.0, 0.0],\
                               "Y": [0.0, 0.0]}, columns=["site","A","C",\
                                                          "D","E","F","G",\
                                                          "H","I","K","L",\
                                                          "M","N","P","Q",\
                                                          "R","S","T","V",\
                                                          "W","Y"])
    >>> align1 = [("seq_1", "ATGATG"), ("seq_2", "CTTATG")]
    >>> align2 = [("seq_1", "ATGATG"), ("seq_2", "CTT---")]
    >>> answer1 = answer.copy()
    >>> answer1[["L", "M"]] = pd.DataFrame({"L": [0.5, 0.0],\
                                            "M": [0.5, 1.0]})
    >>> answer1.equals(calc_aa_frequencies(align1))
    True
    >>> answer1.equals(calc_aa_frequencies(align2))
    True

    """
    # Read in the alignnment
    assert scipy.all(scipy.array([len(s[1]) % 3 for s in alignment]) == 0),\
        "At least one sequence in the alignment is not a multiple of 3."
    seqlength = len(alignment[0][1]) // 3
    df = {k: [0 for x in range(seqlength)] for k in list(INDEX_TO_AA.values())}

    # count amino acid frequencies
    for seq in alignment:
        for i in range(seqlength):
            codon = seq[1][3 * i: 3 * i + 3]
            if codon != '---':
                df[CODONSTR_TO_AASTR[codon]][i] += 1
    df = pd.DataFrame(df)

    # Normalize the dataframe
    assert not scipy.any(df.sum(axis=1) == 0), ("Attempting to normalize a "
                                                "site by an amino acid count"
                                                " of zero. Does the alignment"
                                                " have an all gap column?")
    df = df.div(df.sum(axis=1), axis=0)
    assert scipy.allclose(df.sum(axis=1), 1, atol=0.005)

    # Final formatting
    aa = [x for x in INDEX_TO_AA.values()]
    aa.sort()  # ABC order
    final_cols = ["site"]
    final_cols.extend(aa)
    df["site"] = [x+1 for x in range(len(df))]
    df = df[final_cols]
    return df


def prefDistance(pi1, pi2, distmetric, check_input=True):
    """Compute the distance between two arrays of preferences.

    Args:
        `pi1` and `pi2` (array-like)
            Two arrays of preferences.
        `distmetric` (string)
            Distance metric to use. Can be:
                - `half_sum_abs_diff`: half sum absolute value of difference
                - `JensenShannon`: square root of Jensen-Shannon divergence
                - `RMSD`: root mean square distances

    Returns:
        The distance between `pi1` and `pi2`.

    >>> pi1 = scipy.array([0.5, 0.2, 0.3])
    >>> pi2 = scipy.array([0.2, 0.4, 0.4])
    >>> scipy.allclose(prefDistance(pi1, pi1, 'half_sum_abs_diff'), 0)
    True
    >>> scipy.allclose(prefDistance(pi1, pi1, 'JensenShannon'), 0)
    True
    >>> scipy.allclose(prefDistance(pi1, pi1, 'RMSD'), 0)
    True
    >>> scipy.allclose(prefDistance(pi1, pi2, 'half_sum_abs_diff'), 0.3)
    True
    >>> scipy.allclose(prefDistance(pi1, pi2, 'JensenShannon'), 0.2785483)
    True
    >>> scipy.allclose(prefDistance(pi1, pi2, 'RMSD'), 0.2160245)
    True

    """
    if check_input:
        assert len(pi1) == len(pi2)
        assert scipy.allclose(pi1.sum(), 1, atol=0.005)
        assert scipy.allclose(pi2.sum(), 1, atol=0.005)
        assert scipy.all(pi1 >= 0)
        assert scipy.all(pi2 >= 0)

    if distmetric == 'half_sum_abs_diff':
        dist = (scipy.fabs(pi1 - pi2)).sum() / 2.0

    elif distmetric == 'JensenShannon':
        dist = math.sqrt(divJensenShannon(pi1, pi2))
    elif distmetric == "RMSD":
        dist = math.sqrt(scipy.square(pi1 - pi2).mean())
    else:
        raise ValueError('Invalid `distmetric` {0}'.format(distmetric))

    return dist


def divJensenShannon(p1, p2):
    """Jensen-Shannon divergence between two distributions.

    The logarithms are taken to base 2, so the result will be
    between 0 and 1.
    Args:
        `p1` and `p2` (array-like)
            The two distributions for which we compute divergence.
    Returns:
        The Jensen-Shannon divergence as a float.
    >>> p1 = scipy.array([0.5, 0.2, 0.2, 0.1])
    >>> p2 = scipy.array([0.4, 0.1, 0.3, 0.2])
    >>> p3 = scipy.array([0.0, 0.2, 0.2, 0.6])
    >>> scipy.allclose(divJensenShannon(p1, p1), 0, atol=1e-5)
    True
    >>> scipy.allclose(divJensenShannon(p1, p2), 0.035789, atol=1e-5)
    True
    >>> scipy.allclose(divJensenShannon(p1, p3), 0.392914, atol=1e-5)
    True

    """
    def _kldiv(a, b):
        with scipy.errstate(all='ignore'):
            kl = a * scipy.log2(a / b)
            kl = scipy.nansum(kl)
        return kl

    m = 0.5 * (p1 + p2)

    return 0.5 * (_kldiv(p1, m) + _kldiv(p2, m))


def calc_stationary_state_freqs(model):
    """Calculate the stationary state amino-acids frequencies of a model.

    Args:
        `model` (phydmslib.models.ExpCM)

    Returns:
        frequencies (`numpy.ndarray` of floats)
            The stationary state amino-acid frequencies.
            frequencies[r][a] is the statioanry state frequence of amino acid
            `a` at site `r`.

    """
    def _calc_aa_freq(aa, ss):
        """Calc the frequency of a single site/amino acid pair."""
        codon_indices = scipy.where(CODON_TO_AA == aa)[0]
        return ss[codon_indices].sum()

    frequencies = []
    for r in range(model.nsites):
        aminoacid_ss = scipy.array([_calc_aa_freq(aa, model.stationarystate[r])
                                    for aa in range(N_AA)])
        frequencies.append(aminoacid_ss)
    return scipy.array(frequencies)


def make_expcm(model_fname, prefs):
    """Make an ExpCM from a model params file."""
    params = pd.read_csv(model_fname, engine="python", sep=" = ", header=None)
    params = dict(zip(params[0], params[1]))
    params["phiT"] = 1 - sum([params[x] for x in params.keys()
                             if x.startswith("phi")])
    phi = scipy.array([params["phi{0}".format(INDEX_TO_NT[x])]
                       for x in range(N_NT)])
    return phydmslib.models.ExpCM(prefs, kappa=params["kappa"],
                                  omega=params["omega"], beta=params["beta"],
                                  mu=0.3, phi=phi, freeparams=['mu'])


def make_YNGKP_M0(model_fname, nsites):
    """Make an YNGKP_M0 from a model params file."""
    params = pd.read_csv(model_fname, engine="python", sep=" = ", header=None)
    params = dict(zip(params[0], params[1]))
    e_pw = scipy.zeros((3, N_NT))
    for key in params.keys():
        if key.startswith("phi"):
            p = int(key[-2])
            w = int(NT_TO_INDEX[key[-1]])
            e_pw[p][w] = params[key]
        elif key == "kappa":
            kappa = params[key]
        elif key == "omega":
            omega = params[key]
        else:
            raise ValueError("Unexpected parameter {0}".format(key))
    for p in range(3):
        e_pw[p][3] = 1 - e_pw[p].sum()
    return phydmslib.models.YNGKP_M0(e_pw, nsites, kappa=kappa, omega=omega,
                                     mu=1.0, freeparams=['mu'])


def calculate_pvalue(simulation_values, true_value, seed=False):
    """Calculate pvalue based on simuation distribution.

    The p value is defined as (# simulations greater than true + 1) /
    (# simulations +1).

    In the case where there is at least one simulation with the exact
    same value as the true value, the number of "tied" simulations
    which will be recorded as "greater" will be randomly determined.

    Args:
        `simulation_values` (list)
            List of simulation values.
        `true_value` (`float`)
            True value to calculate p value for.
        `seed` (`int` or `False`)
            Seed used to randomly break ties.
    Returns:
        The p value as a float.
    >>> true = 10
    >>> print("{:1.1f}".format(calculate_pvalue([1, 2, 3, 4], true)))
    0.2
    >>> print("{:1.1f}".format(calculate_pvalue([11, 12, 13, 14], true)))
    1.0
    >>> print("{:1.1f}".format(calculate_pvalue([3, 4, 12, 9], true)))
    0.4
    >>> print("{:1.1f}".format(calculate_pvalue([1, 10, 10, 11], true, 1)))
    0.6

    """
    if seed is not False:
        scipy.random.seed(seed)
    assert len(simulation_values) >= 2, "Must have at least two simulations."
    greater = scipy.sum(scipy.greater(simulation_values, true_value))
    tie_breaker = scipy.sum(scipy.equal(simulation_values, true_value))
    if tie_breaker >= 1:
        tie_breaker = scipy.random.randint(tie_breaker, size=1)[0]
    pvalue = (greater + tie_breaker + 1) / (len(simulation_values) + 1)
    assert 0 <= pvalue <= 1.0, "pvalue is > 1.0 or < 0.0"
    return pvalue


def run_simulations(tree, model, nsim=1000, seed=0, ncpus=1):
    """Simualate a batch of alignments.

    Args:
        `tree` (`Bio.Phylo` tree object)
            Tree to simulate alignments along. Branch lengths are assumed
            to be in units of average codon substitutions / site.
        `model` (`phydmslib.model` object)
            Model to simulate alignments under.
        `nsim` (`int`)
            Number of replicate simualations.
        `seed` (`int`)
            Random seed for simulations. Each simulation will be `seed + i`
            where `i` is the replicate index.
        `ncpus` (`int`)
            Number of cpus for the simulations.
    Returns:
        List of alignments (length `nsim`). Each element is a `list` of
        `tuples` in the form (seq_id, seq).

    """
    # build the simulator object
    simulator = phydmslib.simulate.Simulator(tree, model)
    seeds = [seed+i for i in range(nsim)]

    # run the simulations
    if ncpus == 1:
        simulations = map(simulator.simulate, seeds)
        simulations = list(simulations)
    elif ncpus > 1:  # use multiprocessing
        manager = multiprocessing.Manager()
        sims = manager.dict()
        seed_batches = scipy.array_split(seeds, ncpus)
        assert len(seed_batches) == ncpus
        processes = []
        for i, seed_batch in enumerate(seed_batches):
            p = multiprocessing.Process(target=runSimulator, args=(simulator,
                                                                   seed_batch,
                                                                   i,
                                                                   sims))
            processes.append(p)
            p.start()
        for process in processes:
            process.join()
        # unpack the simulations from `Manager`
        simulations = []
        for key in sims.keys():
            simulations.extend(sims[key])
    else:
        raise ValueError("Unexpected number of cpus ({0})".format(ncpus))
    len(simulations) == nsim
    return simulations


def runSimulator(simulator, seed_list, i, return_dict):
    """Run `simulator.simulate` for a list of seeds."""
    seed_list = [int(x) for x in seed_list]
    simulations = list(map(simulator.simulate, seed_list))
    return_dict[i] = simulations


def process_simulations(simulations, nsites):
    """Reformat simulations as site-specific amino-acid frequencies.

    Args:
        `simulations` (`list`)
            List of alignments (length `nsim`). Each element is a `list` of
            `tuples` in the form (seq_id, seq).
        `nsites` (`int`)
            Number of sites in an alignment.

    Returns:
        Simulation amino-acid frequencies in the shape `(nsites, nsim, N_AA)`
        with `simulations_by_site[r][i][x]` as the frequency of amino acid `x`
        at site `r` in simulation replicate `i`.

    """
    # process
    simulations_by_site = [[] for x in range(nsites)]
    for rep in simulations:
        sim_freqs = calc_aa_frequencies(rep)
        sim_freqs = sim_freqs[amino_acids].values
        for site in range(nsites):
            simulations_by_site[site].append(sim_freqs[site])
    return simulations_by_site


def calc_distances_pvalues(simulations_by_site, model, alignment, seed,
                           metrics=["JensenShannon", "half_sum_abs_diff", "RMSD"]):
    """Calculate distances and pvalues.

    Args:
        `simulations_by_site` (`list`)
            Simulation amino-acid frequencies in the shape
            `(nsites, nsim, N_AA)` with `simulations_by_site[r][i][x]` as the
            frequency of amino acid `x` at site `r` in simulation replicate `i`
        `model` (`phydmslib.model` object )
        `alignment` (list)
            Alignment of codon sequences as a list of tuples, (seq_id, seq)
        `seed` (`int`)
            Random seed to break ties. The random seed for a given site is
            `seed` + `site`.
        `metrics` (`list` of `str`)
            Distance metrics. Must be allowed in `prefDistance`.

    Returns:
        `pandas` dataframe of length `nsites` with columns `site`, `pvalue`,
        and `metric`

    """
    # model statioanry state frequencies
    ss_freqs = calc_stationary_state_freqs(model)
    nsites = model.nsites
    assert len(simulations_by_site) == nsites

    # natural alignment frequencies
    alignment_freqs = calc_aa_frequencies(alignment)
    alignment_freqs["alignment"] = (alignment_freqs[amino_acids]
                                    .apply(lambda r: tuple(r), axis=1)
                                    .apply(scipy.array))
    alignment_freqs = scipy.array(alignment_freqs["alignment"].tolist())
    assert len(alignment_freqs) == nsites

    # pvalues
    pvalues = []
    for site in range(nsites):
        sims = simulations_by_site[site]
        ss = ss_freqs[site]
        natural = alignment_freqs[site]
        for metric in metrics:
            # distance from model to simulations
            sim = scipy.array([prefDistance(sim, ss, metric, check_input=False)
                               for sim in sims])
            # distance from model to natural alignment
            nat_dist = prefDistance(natural, ss, metric, check_input=False)
            pvalue = calculate_pvalue(sim, nat_dist, seed=seed + site)
            pvalues.append((site+1, pvalue, metric))
    pvalues = pd.DataFrame(pvalues, columns=['site', 'pvalue', 'metric'])
    return pvalues


def calc_qvalues(pvalues):
    """Calculate qvalues.

    Args:
        `pvalues` (`pandas` dataframe)
            `pandas` dataframe of length `nsites` with columns `site`,
            `pvalue`, and `metric`

    Returns:
        `pandas` dataframe of length `nsites` with columns `site`, `pvalue`,
         `metric`, and `qvalue`

    """
    final = []
    for name, group in pvalues.groupby("metric"):
        group = group.sort_values(by="pvalue", ascending=True)
        group["qvalue"] = multipletests(group["pvalue"].tolist(),
                                        method='fdr_bh')[1]
        final.append(group)
    return pd.concat(final)


def calc_ncpus(cpus, njobs, min_jobs_per_cpu=1000):
    """Calcluate the number of CPUs to use.

    Args:
        `cpus` (`int`)
            Number of CPUs available. If -1, then a number was not specified
            by the user.
        `njobs` (`int`)
            Number of jobs to run.
        `min_jobs_per_cpu` (`int`)
            Minimum number of jobs per CPU.

    Returns:
        `ncpus` (`int`)
            number of CPUs
        `msg` (`str`)
            message for logger
    >>> print(calc_ncpus(10, 10000, 1000)[0])
    10
    >>> print(calc_ncpus(10, 10000, 8000)[0])
    1

    """
    if cpus == -1:
        ncpus = multiprocessing.cpu_count()
    else:
        ncpus = cpus
    assert ncpus >= 1, "{0} CPUs specified".format(ncpus)
    if (njobs/ncpus) < min_jobs_per_cpu:
        ncpus = max(1, njobs // min_jobs_per_cpu)
        msg = ('Using {0} CPU(s) to ensure >= {1} simulations / CPU'
               .format(ncpus, min_jobs_per_cpu))
    else:
        msg = 'Using {0} CPUs.'.format(ncpus)
    return ncpus, msg


if __name__ == '__main__':
    import doctest
    doctest.testmod()
