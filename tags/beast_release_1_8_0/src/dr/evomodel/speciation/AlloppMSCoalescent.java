package dr.evomodel.speciation;


import dr.evolution.util.Units;
import dr.inference.model.Likelihood;


/**
 * Computes coalescent log-likelihood of a set of gene trees embedded inside a 
 * allopolyploid species network.
 *
 * @author Graham Jones
 *         Date: 03/06/2011
 */


/*
 * AlloppMSCoalescent, AlloppSpeciesNetworkModel, AlloppLeggedTree,
 * AlloppDiploidHistory, and AlloppSpeciesBindings collaborate closely.
 * 
 * AlloppSpeciesNetworkModel contains an AlloppDiploidHistory and an array of
 * AlloppLeggedTree's for representing the tetraploid trees. It also contains a
 * multiply labelled tree as an alternative representation.
 * 
 * AlloppSpeciesBindings represents the species-individual-genome structure
 * of the data, and contains the gene trees.
 * 
 * This module, AlloppMSCoalescent, uses AlloppSpeciesNetworkModel and 
 * AlloppSpeciesBindings to calculate P(g_i|S) = probability that gene tree
 * g_i fits into network S.
 * 
 * apsp.geneTreeFitsInNetwork calls
 *   FitsInNetwork (inside a gene tree) calls
 *     subTreeFitsInNetwork (recursive, inside a gene tree) calls
 *        CoalescenceIsCompatible to test a single coalescence against mullab tree
 *        
 * apsp.geneTreeLogLikelihood calls
 *   TreeLogLikelihood  (inside a gene tree) calls
 *     clearCoalescences (in mullab tree) and
 *     recordCoalescence (in mullab tree) and
 *     recordLineageCounts (in mullab tree) and
 *     to set up data in mullab tree nodes, then
 *     mullabTreeLogLikelihood which calls
 *       mullabSubTreeLogLikelihood (recursive, in mullab tree) calls
 *         limbLogLike calls
 *           limbLinPopIntegral
 *        
 */




public class AlloppMSCoalescent extends Likelihood.Abstract implements Units {
    private final AlloppSpeciesNetworkModel asnetwork;
    private final AlloppSpeciesBindings apsp;
	
    
    
    public AlloppMSCoalescent(AlloppSpeciesBindings apspecies, AlloppSpeciesNetworkModel apspnetwork) {
        super(apspnetwork);
        apsp = apspecies;
        asnetwork = apspnetwork;
        
        asnetwork.addModelListener(this);
        apsp.addModelListeners(this);
    }

    
    
    @Override
	protected double calculateLogLikelihood() {
    	for (int i = 0; i < apsp.numberOfGeneTrees(); i++) {
    		if (!apsp.geneTreeFitsInNetwork(i, asnetwork)) {
    			return Double.NEGATIVE_INFINITY;
    		}
    	}
	    // grjtodo-oneday JH has compatible flags for efficiency. I'm checking
	    // every time.
    	
        double logl = 0;
        for(int i = 0; i < apsp.numberOfGeneTrees(); i++) {
            final double v = apsp.geneTreeLogLikelihood(i, asnetwork);
            assert ! Double.isNaN(v);
            logl += v;
        }
        return logl;
    }


	@Override
	protected boolean getLikelihoodKnown() {
		return false;
	}

    
	public Type getUnits() {
		return asnetwork.getUnits();
	}

	public void setUnits(Type units) {
		// TODO Auto-generated method stub
        // one day may allow units other than substitutions

	}


	
}
