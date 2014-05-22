package dr.inference.operators;

import dr.evomodel.continuous.LatentFactorModel;

/**
 * Created with IntelliJ IDEA.
 * User: max
 * Date: 5/22/14
 * Time: 12:49 PM
 * To change this template use File | Settings | File Templates.
 */
public class FactorGibbsOperator extends SimpleMCMCOperator implements GibbsOperator{
    private LatentFactorModel LFM;
    public FactorGibbsOperator(LatentFactorModel LFM, double weight){
        this.LFM=LFM;
        setWeight(weight);
    }


    public int getStepCount() {
        return 0;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public String getPerformanceSuggestion() {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public String getOperatorName() {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public double doOperation() throws OperatorFailedException {
        return 0;  //To change body of implemented methods use File | Settings | File Templates.
    }


}
