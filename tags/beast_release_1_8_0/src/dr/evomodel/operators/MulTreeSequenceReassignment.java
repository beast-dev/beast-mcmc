
package dr.evomodel.operators;

import dr.evomodel.speciation.MulSpeciesBindings;
import dr.evomodel.speciation.MulSpeciesTreeModel;
import dr.evomodelxml.operators.MulTreeSequenceReassignmentParser;
import dr.inference.operators.OperatorFailedException;
import dr.inference.operators.SimpleMCMCOperator;
import dr.math.MathUtils;


/**
 * Changes sequence assignments for one gene, one or more individuals.
 * Very similar to AlloppSequenceReassignment.
 * 
 * @author Graham Jones
 *         Date: 20/12/2011
 */


public class MulTreeSequenceReassignment extends SimpleMCMCOperator {

	private final MulSpeciesTreeModel multree;
	private final MulSpeciesBindings mulspb;


	public MulTreeSequenceReassignment(MulSpeciesTreeModel multree, MulSpeciesBindings mulspb, double weight) {
		this.multree = multree;
		this.mulspb = mulspb;
		setWeight(weight);
	}	


	public String getPerformanceSuggestion() {
		return "None";
	}

	@Override
	public String getOperatorName() {
		return MulTreeSequenceReassignmentParser.MULTREE_SEQUENCE_REASSIGNMENT + "(" + multree.getId() +
		"," + mulspb.getId() + ")";
	}

	@Override
	public double doOperation() throws OperatorFailedException {
		multree.beginTreeEdit();
		if (MathUtils.nextInt(2) == 0) {
			mulspb.permuteOneSpeciesOneIndivForOneGene();
		} else {
			mulspb.permuteSetOfIndivsForOneGene();
		}
		
		multree.endTreeEdit();
		return 0;
	}

}


