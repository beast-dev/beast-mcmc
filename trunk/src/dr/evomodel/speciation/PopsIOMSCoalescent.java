package dr.evomodel.speciation;

import dr.inference.model.Likelihood;

/**
 * User: Graham Jones
 * Date: 10/05/12
 */

public class PopsIOMSCoalescent extends Likelihood.Abstract {
    private final PopsIOSpeciesTreeModel piostm;
    private final PopsIOSpeciesBindings piosb;


    public PopsIOMSCoalescent(PopsIOSpeciesBindings piosb, PopsIOSpeciesTreeModel piostm) {
        super(piostm);
        this.piostm = piostm;
        this.piosb = piosb;

        piostm.addModelRestoreListener(this);

        final PopsIOSpeciesBindings.GeneTreeInfo[] trees = piosb.getGeneTrees();
        for(PopsIOSpeciesBindings.GeneTreeInfo geneTree : trees) {
            geneTree.tree.addModelListener(this);
        }

    }


    @Override
    protected double calculateLogLikelihood() {
        for (int i = 0; i < piosb.numberOfGeneTrees(); i++) {
            if (!piosb.geneTreeFitsInNetwork(i, piostm)) {
                return Double.NEGATIVE_INFINITY;
            }
        }
        // grjtodo-oneday JH has compatible flags for efficiency. I'm checking
        // every time.

        double logl = 0;
        for(int i = 0; i < piosb.numberOfGeneTrees(); i++) {
            final double v = piosb.geneTreeLogLikelihood(i, piostm);
            assert ! Double.isNaN(v);
            logl += v;
        }
        return logl;
    }


    @Override
    protected boolean getLikelihoodKnown() {
        return false;
    }
}
