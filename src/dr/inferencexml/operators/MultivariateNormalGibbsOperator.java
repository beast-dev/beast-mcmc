package dr.inferencexml.operators;

import com.sun.tools.javac.code.Attribute;
import dr.evomodel.continuous.AbstractMultivariateTraitLikelihood;
import dr.inference.distribution.MultivariateDistributionLikelihood;
import dr.inference.distribution.MultivariateNormalDistributionModel;
import dr.inference.model.CompoundParameter;
import dr.inference.model.MatrixParameter;
import dr.math.matrixAlgebra.Matrix;
import dr.inference.model.Parameter;
import dr.inference.operators.OperatorFailedException;
import dr.inference.operators.SimpleMCMCOperator;
import dr.inference.operators.GibbsOperator;
import dr.math.matrixAlgebra.Vector;

/**
 * Created with IntelliJ IDEA.
 * User: max
 * Date: 2/27/14
 * Time: 4:21 PM
 * To change this template use File | Settings | File Templates.
 */
public class MultivariateNormalGibbsOperator extends SimpleMCMCOperator implements GibbsOperator {
    private MatrixParameter priorPrecision;
    private Parameter priorMean;
    private MatrixParameter likelihoodPrecision;
    private Parameter likelihoodMean;
    private CompoundParameter data;


    MultivariateNormalGibbsOperator(CompoundParameter data, MultivariateNormalDistributionModel likelihood, MultivariateDistributionLikelihood prior, Double weight){
        this.data=data;
//        this.likelihoodPrecision=new MatrixParameter(likelihood.getDistribution().getScaleMatrix());
        double[] priorMean= prior.getDistribution().getMean();
//        setParameterValue(this.priorMean, priorMean);

//        this.priorPrecision=new Matrix(prior.getDistribution().getScaleMatrix());
        setWeight(weight);
    }

    private void setParameterValue(Parameter set, double[] value){
      set.setDimension(value.length);
        for(int i=0; i<value.length; i++)
        {set.setParameterValueQuietly(i,value[i]);}
        set.fireParameterChangedEvent();
    }


    @Override
    public String getOperatorName() {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public String getPerformanceSuggestion() {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public double doOperation() throws OperatorFailedException {
        return 0;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public int getStepCount() {
        return 1;  //To change body of implemented methods use File | Settings | File Templates.
    }
}
