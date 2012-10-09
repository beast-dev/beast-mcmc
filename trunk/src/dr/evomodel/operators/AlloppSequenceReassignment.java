package dr.evomodel.operators;

import dr.evolution.util.Taxon;
import dr.evomodel.speciation.AlloppLeggedTree;
import dr.evomodel.speciation.AlloppSpeciesBindings;
import dr.evomodel.speciation.AlloppSpeciesNetworkModel;
import dr.evomodelxml.operators.AlloppSequenceReassignmentParser;
import dr.inference.operators.OperatorFailedException;
import dr.inference.operators.SimpleMCMCOperator;
import dr.math.MathUtils;

import java.util.ArrayList;


/**
 * 
 * @author Graham Jones
 *         Date: 01/07/2011
 */


public class AlloppSequenceReassignment extends SimpleMCMCOperator {

	private final AlloppSpeciesNetworkModel apspnet;
	private final AlloppSpeciesBindings apsp;


	public AlloppSequenceReassignment(AlloppSpeciesNetworkModel apspnet, AlloppSpeciesBindings apsp, double weight) {
		this.apspnet = apspnet;
		this.apsp = apsp;
		setWeight(weight);
	}	


	public String getPerformanceSuggestion() {
		return "None";
	}

	@Override
	public String getOperatorName() {
		return AlloppSequenceReassignmentParser.SEQUENCE_REASSIGNMENT + "(" + apspnet.getId() +
		"," + apsp.getId() + ")";
	}

	@Override
	public double doOperation() throws OperatorFailedException {
		apspnet.beginNetworkEdit();

        if (MathUtils.nextInt(10) == 0) {
            int tt = MathUtils.nextInt(apspnet.getNumberOfTetraTrees());
            AlloppLeggedTree ttree = apspnet.getTetraploidTree(tt);
            ArrayList<Taxon> sptxs = ttree.getSpeciesTaxons();
            for (Taxon tx : sptxs) {
                int spi = apsp.apspeciesId2index(tx.getId());
                apsp.flipAssignmentsForAllGenesOneSpecies(spi);
            }
            apspnet.flipLegsOfTetraTree(tt);
        } else {
            if (MathUtils.nextInt(2) == 0) {
                apsp.permuteOneSpeciesOneIndivForOneGene();
            } else {
                apsp.permuteSetOfIndivsForOneGene();
            }
        }

		
		apspnet.endNetworkEdit();
        assert apspnet.alloppspeciesnetworkOK();
		return 0;
	}

}


