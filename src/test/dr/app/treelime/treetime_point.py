# TreeTime as a point estimate, for TopologySearchStudy's pipeline Fb: the treetime command line as it is, but with
# zero-length branches kept (prune_short=False; the command line always prunes them into polytomies) and, with
# --keep-polytomies, no polytomy resolution, so the binary input tree comes back binary: TreeTime then only roots
# and dates it. Run with the Python of TreeTime's environment, with treetime's own arguments.
import sys
from treetime import wrappers
from treetime.__main__ import main

run_timetree = wrappers.run_timetree


def run_timetree_unpruned(myTree, params, outdir, tree_suffix='', prune_short=True, method_anc='probabilistic'):
    return run_timetree(myTree, params, outdir, tree_suffix=tree_suffix, prune_short=False, method_anc=method_anc)


wrappers.run_timetree = run_timetree_unpruned

if __name__ == '__main__':
    sys.argv[0] = 'treetime'
    sys.exit(main())
