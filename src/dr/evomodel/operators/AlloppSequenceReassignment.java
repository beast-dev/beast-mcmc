package dr.evomodel.operators;

import dr.evomodel.speciation.AlloppSpeciesBindings;
import dr.evomodel.speciation.AlloppSpeciesNetworkModel;
import dr.evomodelxml.operators.AlloppSequenceReassignmentParser;
import dr.inference.operators.OperatorFailedException;
import dr.inference.operators.SimpleMCMCOperator;
import dr.math.MathUtils;


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
		if (MathUtils.nextInt(2) == 0) {
			apsp.permuteOneSpeciesOneIndivForOneGene();
		} else {
			apsp.permuteSetOfIndivsForOneGene();
		}
		
		apspnet.endNetworkEdit();
		return 0;
	}

}


