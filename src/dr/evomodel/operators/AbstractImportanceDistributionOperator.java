/**
 * 
 */
package dr.evomodel.operators;

import dr.evomodel.tree.TreeModel;
import dr.inference.model.Parameter;
import dr.inference.operators.CoercionMode;
import dr.inference.operators.MCMCOperator;
import dr.inference.operators.OperatorFailedException;
import dr.inference.operators.OperatorSchedule;
import dr.inference.operators.ScaleOperator;
import dr.inference.operators.SimpleMCMCOperator;
import dr.inference.operators.SimpleOperatorSchedule;
import dr.inference.operators.UniformOperator;

/**
 * @author Sebastian Hoehna
 *
 */
public abstract class AbstractImportanceDistributionOperator extends SimpleMCMCOperator {

	private int transitions = 0;
	
	private OperatorSchedule schedule;
	
	protected TreeModel tree;
	

	/**
	 * 
	 */
	public AbstractImportanceDistributionOperator(TreeModel tree) {
		super();
		
		this.tree = tree;
		
		init();
	}
	
	private void init(){
		schedule = getOperatorSchedule(tree);
	}
	
	private OperatorSchedule getOperatorSchedule(TreeModel treeModel) {

        ExchangeOperator narrowExchange = new ExchangeOperator(ExchangeOperator.NARROW, treeModel, 10);
        ExchangeOperator wideExchange = new ExchangeOperator(ExchangeOperator.WIDE, treeModel, 3);
        SubtreeSlideOperator subtreeSlide = new SubtreeSlideOperator(treeModel, 10.0, 1.0, true, false, false, false, CoercionMode.COERCION_ON);
        NNI nni = new NNI(treeModel, 10.0);
        WilsonBalding wilsonBalding = new WilsonBalding(treeModel, null, 3.0);
        FNPR fnpr = new FNPR(treeModel, 5.0);

        OperatorSchedule schedule = new SimpleOperatorSchedule();
        schedule.addOperator(narrowExchange);
        schedule.addOperator(wideExchange);
        schedule.addOperator(subtreeSlide);
        schedule.addOperator(nni);
        schedule.addOperator(wilsonBalding);
        schedule.addOperator(fnpr);

        return schedule;
    }
	
	protected double doUnguidedOperation() throws OperatorFailedException{
		int index = schedule.getNextOperatorIndex();
		SimpleMCMCOperator operator = (SimpleMCMCOperator) schedule.getOperator(index);
		
		return operator.doOperation();
	}
	
	/**
     * @return the number of transitions since last call to reset().
     */
    public int getTransitions(){
    	return transitions;
    }

    /**
     * Set the number of transitions since last call to reset(). This is used
     * to restore the state of the operator
     *
     * @param rejected number of rejections
     */
    public void setTransitions(int transitions){
    	this.transitions = transitions;
    }
    
    public double getTransistionProbability() {
        int accepted = getAccepted();
        int rejected = getRejected();
        int transition = getTransitions();
        return (double) transition / (double) (accepted + rejected);
    }
    
    public double getMinimumAcceptanceLevel() {
        return 0.50;
    }

    public double getMaximumAcceptanceLevel() {
        return 1.0;
    }

    public double getMinimumGoodAcceptanceLevel() {
        return 0.75;
    }

    public double getMaximumGoodAcceptanceLevel() {
        return 1.0;
    }

	/* (non-Javadoc)
	 * @see dr.inference.operators.SimpleMCMCOperator#doOperation()
	 */
	@Override
	public abstract double doOperation() throws OperatorFailedException;

	/* (non-Javadoc)
	 * @see dr.inference.operators.SimpleMCMCOperator#getOperatorName()
	 */
	@Override
	public abstract String getOperatorName();

	/* (non-Javadoc)
	 * @see dr.inference.operators.MCMCOperator#getPerformanceSuggestion()
	 */
	public abstract String getPerformanceSuggestion();

}
