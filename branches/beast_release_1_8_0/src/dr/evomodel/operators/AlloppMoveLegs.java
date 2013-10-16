package dr.evomodel.operators;

import dr.evomodel.speciation.AlloppSpeciesBindings;
import dr.evomodel.speciation.AlloppSpeciesNetworkModel;
import dr.evomodelxml.operators.AlloppMoveLegsParser;
import dr.inference.operators.OperatorFailedException;
import dr.inference.operators.SimpleMCMCOperator;


/**
 * 
 * @author Graham Jones
 *         Date: 31/08/2011
 */


// 2012-07 ood
public class AlloppMoveLegs extends SimpleMCMCOperator {

	private final AlloppSpeciesNetworkModel apspnet;
	private final AlloppSpeciesBindings apsp;
	
	public AlloppMoveLegs(AlloppSpeciesNetworkModel apspnet, AlloppSpeciesBindings apsp, double weight) {
		this.apspnet = apspnet;
		this.apsp = apsp;
		setWeight(weight);
	}	

	
	public String getPerformanceSuggestion() {
		return "None";
	}

	@Override
	public String getOperatorName() {
		return AlloppMoveLegsParser.MOVE_LEGS + "(" + apspnet.getId() +
		"," + apsp.getId() + ")";
	}

	@Override
	public double doOperation() throws OperatorFailedException {
		apspnet.beginNetworkEdit();
		apspnet.moveLegs();		
		apspnet.endNetworkEdit();
		return 0;
	}

}
