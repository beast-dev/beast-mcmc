package dr.evomodel.speciation;

import dr.inference.model.Likelihood;


/**
 * 
 * Calculates prior likelihood for an allopolyploid network.
 * 
 * @author Graham Jones
 *         Date: 01/07/2011
 */





public class AlloppNetworkPrior extends Likelihood.Abstract {
	AlloppSpeciesNetworkModel asnm;
	AlloppNetworkPriorModel prior;

	public AlloppNetworkPrior(AlloppNetworkPriorModel prior, AlloppSpeciesNetworkModel asnm) {
		super(prior);
		this.asnm = asnm;
		this.prior = prior;
	}

	@Override
	protected boolean getLikelihoodKnown() {
		return false;
	}

	@Override
	protected double calculateLogLikelihood() {
		double rate = prior.getRate().getParameterValue(0);
		double heights[] = asnm.getEventHeights();
		int s = asnm.getTipCount();
		double loglhood = 0.0;
		if (s > 1) {
			loglhood += Math.log(s*(s-1)); 
		}
		loglhood += s*Math.log(rate);
		for (int i = 0; i < heights.length; i++) {
			loglhood -= rate*(heights[i]);
		}			
		// grjtodo morethanonetree. Currently for one tet tree case, simple version.
		return loglhood;
	}	
}
