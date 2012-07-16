package dr.evomodel.speciation;

import dr.inference.model.Likelihood;

/**
 * User: Graham Jones
 * Date: 10/05/12
 */
public class PopsIOSpeciesTreePrior  extends Likelihood.Abstract {
    private SpeciationModel sppm;
    private PopsIOSpeciesTreeModel piostm;

    public PopsIOSpeciesTreePrior(SpeciationModel sppm, PopsIOSpeciesTreeModel piostm) {
        super(sppm);
        this.piostm = piostm;
        this.sppm = sppm;
    }

    @Override
    protected double calculateLogLikelihood() {
        double lhood = sppm.calculateTreeLogLikelihood(piostm);
        return lhood;
    }

    @Override
    protected boolean getLikelihoodKnown() {
        return false;
    }
}
