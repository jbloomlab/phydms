"""Profiles ``phydms`` on a small-ish data set.

Runs ``phydms`` with a series of common options,
using the ``--profile`` option to profile the
computational time for the likelihood maximization.
"""


import os
import glob
import time
import subprocess
import multiprocessing


def main():
    """Main body of script."""

    profiledir = './NP_profile/'
    if not os.path.isdir(profiledir):
        os.mkdir(profiledir)
    print("Profiling results will be written to {0}".format(profiledir))

    alignment = './NP_data/NP_alignment.fasta'
    tree = './NP_data/NP_tree.newick'
    prefs = './NP_data/NP_prefs.tsv'
    print("Using alignment, tree, and prefs in:\n\t{0}".format('\n\t'.join(
            [alignment, tree, prefs])))

    ncpus = max(1, multiprocessing.cpu_count() - 1)
    print("Will run programs using {0} CPUs.".format(ncpus))

    # run in a multiprocessing pool to run several at same time
    pool = {}
    for (name1, model, baseargs) in [
            ('YNGKP_M0', 'YNGKP_M0', []),
            ('YNGKP_M5', 'YNGKP_M5', []),
            ('ExpCM', 'ExpCM_{0}'.format(prefs), []),
            ('ExpCM_fitphi', 'ExpCM_{0}'.format(prefs), ['--fitphi']),
            ('ExpCM_gammaomega', 'ExpCM_{0}'.format(prefs), ['--gammaomega']),
            ]:
        for (brlen, brlenargs) in [
                ('brlen-scale', ['--brlen', 'scale']),
                ('brlen-optimize', ['--brlen', 'optimize']),
                ]:
            args = baseargs + brlenargs
#            args.append('--omegabysite')
            name = '{0}_{1}'.format(name1, brlen)
            outprefix = '{0}/{1}'.format(profiledir, name)
            for f in glob.glob('{0}*'.format(outprefix)):
                os.remove(f)
            cmd = ['phydms', alignment, tree, model, outprefix, '--profile'] + args
            pool[name] = multiprocessing.Process(target=subprocess.check_call, 
                    args=(cmd,), kwargs={'stdout':subprocess.PIPE, 
                    'stderr':subprocess.PIPE})
    started = dict([(name, False) for name in pool])
    completed = dict([(name, False) for name in pool])

    while not all(completed.values()):
        nrunning = (list(started.values()).count(True) - 
                list(completed.values()).count(True))
        if nrunning < ncpus:
            for (name, p) in pool.items():
                if (not started[name]) and (not completed[name]):
                    p.start()
                    started[name] = True
                    print("Started profiling run for {0}".format(name))
                    break
        for (name, p) in pool.items():
            if started[name] and (not completed[name]) and (
                    not p.is_alive()):
                completed[name] = True
                print("Completed profiling run for {0}".format(name))
        time.sleep(5)

    print("Program complete.")


if __name__ == '__main__':
    main() # run the script
