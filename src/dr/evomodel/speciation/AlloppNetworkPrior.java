package dr.evomodel.speciation;

import dr.evolution.tree.SimpleNode;
import dr.evolution.tree.SimpleTree;
import dr.evolution.tree.Tree;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;


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
	BirthDeathGernhard08Model bdgm;
	double beta = 0.25; // two diploids make allotetraploids at this rate relative
	                            // to rate that one diploid speciates.
	double rho = 1.0; // for hybridization, not (yet) for birth-death-sample. Must be <= 1
	double gamma = 0.8; // mu/lambda
	// One day might estimate some of these values...
	

	public AlloppNetworkPrior(AlloppNetworkPriorModel prior, AlloppSpeciesNetworkModel asnm) {
		super(prior);
		this.asnm = asnm;
		this.prior = prior;
		
		// 2011-09-17 I only use BirthDeathGernhard08Model for testing against.
		Parameter birthDiffRateParameter = prior.getRate();
		Parameter relativeDeathRateParameter = new Parameter.Default(gamma);
		BirthDeathGernhard08Model.TreeType treetype = BirthDeathGernhard08Model.TreeType.ORIENTED;
		dr.evolution.util.Units.Type unittype = dr.evolution.util.Units.Type.SUBSTITUTIONS;
		bdgm = new BirthDeathGernhard08Model(
				birthDiffRateParameter, relativeDeathRateParameter, null,
                treetype, unittype);
		
		SimpleNode [] test = new SimpleNode[9];
		for (int i=0; i < 9; i++) {
			test[i] = new SimpleNode();
		}
		test[0].setId("a");
		test[1].setId("b");
		test[2].addChild(test[0]);
		test[2].addChild(test[1]);		
		test[3].setId("c");
		test[4].setId("d");
		test[5].addChild(test[3]);
		test[5].addChild(test[4]);		
		test[6].setId("e");
		test[7].addChild(test[5]);
		test[7].addChild(test[6]);		
		test[8].addChild(test[2]);
		test[8].addChild(test[7]);		
		test[2].setHeight(0.01);
		test[5].setHeight(0.02);
		test[7].setHeight(0.03);
		test[8].setHeight(0.04);
		SimpleTree testtree = new SimpleTree(test[4]);
		double z = bdgm.calculateTreeLogLikelihood(testtree);
		double w = loglhoodDiploidTree(testtree);
		if (Math.abs(z - w) >= 1e-14) {
			System.err.println("AlloppNetworkPrior.calculateLogLikelihood() numerical error?");
		}
		
	}
	
	

	@Override
	protected boolean getLikelihoodKnown() {
		return false;
	}

	@Override
	protected double calculateLogLikelihood() {
		
		double loglhood = 0.0;
		
		// normal likelihood for diploid tree, with uniform prior on origin. 
		AlloppLeggedTree ditree = asnm.getHomoploidTree(AlloppSpeciesNetworkModel.DITREES, 0);
		loglhood += bdgm.calculateTreeLogLikelihood(ditree);
		double z = loglhoodDiploidTree(ditree);
		if (Math.abs(z - loglhood) >= 1e-12) {
			System.err.println("AlloppNetworkPrior.calculateLogLikelihood() numerical error?");
		}

		int noftetratrees = asnm.getNumberOfTetraTrees();
		for (int i = 0; i < noftetratrees; i++) {
			AlloppLeggedTree tetratree = asnm.getHomoploidTree(AlloppSpeciesNetworkModel.TETRATREES, i);
			// likelihood conditioned on origin = hybridization time for tetraploid trees 
			loglhood += loglhoodLeggedTree(tetratree);
			// not-principled likelihood for hybridization time
			loglhood += loglikelihoodHybridizationTime(ditree, tetratree);
		}
		
		return loglhood;
	}
	
	
	private double loglikelihoodHybridizationTime(AlloppLeggedTree ditree, AlloppLeggedTree tetratree) {
		double rooth = ditree.getRootHeight();
		int ndips = ditree.getExternalNodeCount();
		double delta = prior.getRate().getParameterValue(0);
		double hybh = tetratree.getHybridHeight();
		double lhood = beta * delta / (1.0 - gamma);
		double extantdips = ndips / rho;
		double extantdippairs = 0.5 * extantdips * (extantdips-1);
		lhood *= (hybh/rooth) + ((rooth-hybh)/rooth) * extantdippairs;
		double loglhood = Math.log(lhood);
		return loglhood;
	}
	
	
	private double loglhoodDiploidTree(Tree ditree) {
		int ntips = ditree.getExternalNodeCount();
		double z = Math.log(ntips);
		double delta = prior.getRate().getParameterValue(0);
		z += (ntips-1) * Math.log(delta);
		z -= (ntips-2) * Math.log(1.0-gamma);
 		for (int i = 0; i < ditree.getInternalNodeCount(); i++) {
			double t = ditree.getNodeHeight(ditree.getInternalNode(i));
			if (ditree.isRoot(ditree.getInternalNode(i))) {
			z -= delta * t;
			z -= Math.log( 1.0 - gamma*(expdx(delta, t)) );
			}
			z += logp1(delta, t);
		}
		return z;
	}
	
	
	
	
	private double loglhoodLeggedTree(AlloppLeggedTree tetratree) {
		int ntips = tetratree.getExternalNodeCount();
		if (ntips == 1) {
			return 0.0;
		}
		
		double z = 0.0;
		double h = tetratree.getHybridHeight();
		double delta = prior.getRate().getParameterValue(0);
		double y = logq1(delta, h);
		for (int i = 0; i < tetratree.getInternalNodeCount(); i++) {
			double t = tetratree.getNodeHeight(tetratree.getInternalNode(i));
			z += logp1(delta, t);
			z -= y;
		}
		return z;
	}
	
	
	
	private double expdx(double delta, double x) {
		return Math.exp(-delta * x);
	}
	
	
	private double oneminusexpdx(double delta, double x) {
		assert x >= 0.0;
		double y = delta * x;
		if (y > 1e-6) {
			return 1.0 - expdx(delta, x);
		} else {
			return y - 0.5*y*y;
		}
	}
	
	
	private double logq1(double delta, double x) {
		double z = 0.0;
		z += Math.log(1.0 - gamma);
		z += Math.log(oneminusexpdx(delta, x));
		z -= Math.log(delta);
		z -= Math.log(1.0 - gamma * expdx(delta, x));
		return z;
	}
	
	
	
	private double logp1(double delta, double x) {
		double z = 0.0;
		z +=  2.0 * Math.log(1.0 - gamma);
		z -= delta * x;
		z -= 2 * Math.log(1.0 - gamma * expdx(delta, x));
		return z;
	}
		
	
}
