package dr.evomodel.speciation;

import dr.evolution.util.Units;
import dr.evomodelxml.speciation.AlloppNetworkPriorModelParser;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.inference.model.Variable.ChangeType;


/**
 * 
 * Model for prior of an allopolyploid network.
 * 
 * @author Graham Jones
 *         Date: 01/07/2011
 */




public class AlloppNetworkPriorModel extends AbstractModel  implements Units {
	private Parameter rate;
	private Units.Type units;
	

	public AlloppNetworkPriorModel(Parameter rate, Units.Type units) {
		super(AlloppNetworkPriorModelParser.ALLOPPNETWORKPRIORMODEL);
		this.rate = rate;
		this.units = units;
		addVariable(rate);
	}

	public Parameter getRate() {
		return rate;
	}
	
	@Override
	protected void handleModelChangedEvent(Model model, Object object, int index) {
		if (AlloppSpeciesNetworkModel.DBUGTUNE)
			System.out.println("AlloppNetworkPriorModel.handleModelChangedEvent() " + model.getId());
	}

	@Override
	protected void handleVariableChangedEvent(Variable variable, int index,
			ChangeType type) {
		if (AlloppSpeciesNetworkModel.DBUGTUNE)
			System.out.println("AlloppNetworkPriorModel.handleModelChangedEvent() " + variable.getId());
	}

	@Override
	protected void storeState() {
		// addVariable(rate) deals with this
		if (AlloppSpeciesNetworkModel.DBUGTUNE)
			System.out.println("AlloppNetworkPriorModel.storeState()");

	}

	@Override
	protected void restoreState() {
		// addVariable(rate) deals with this
		if (AlloppSpeciesNetworkModel.DBUGTUNE)
			System.out.println("AlloppNetworkPriorModel.restoreState()");

	}

	@Override
	protected void acceptState() {
		// TODO Auto-generated method stub

	}
	
	public Type getUnits() {
		return units;
	}

	public void setUnits(Type units) {
		// TODO Auto-generated method stub
		
	}	

}
