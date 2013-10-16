package dr.evomodel.speciation;

import dr.inference.model.Likelihood;

/**
 * @author  Graham  Jones
 * Date: 10/05/2012
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

        piosb.fillSpeciesTreeWithCoalescentInfo(piostm);
        return piostm.logLhoodAllGeneTreesInSpeciesTree();
    }


    @Override
    protected boolean getLikelihoodKnown() {
        return false;
    }
}
