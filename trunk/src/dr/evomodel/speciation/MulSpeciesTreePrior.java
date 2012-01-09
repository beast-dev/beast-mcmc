package dr.evomodel.speciation;

import dr.inference.model.Likelihood;



public class MulSpeciesTreePrior  extends Likelihood.Abstract {
	MulSpeciesTreeModel mulsptree;
	SpeciationModel prior;

	public MulSpeciesTreePrior(SpeciationModel prior, MulSpeciesTreeModel mulsptree) {
		super(prior);
		this.mulsptree = mulsptree;
		this.prior = prior;
	}

	
	@Override
	protected boolean getLikelihoodKnown() {
		return false;
	}

	
	
	@Override
	protected double calculateLogLikelihood() {
		double lhood = prior.calculateTreeLogLikelihood(mulsptree);
		return lhood;
	}
}
