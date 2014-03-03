package dr.inferencexml.operators;

import dr.inference.distribution.MultivariateDistributionLikelihood;
import dr.inference.distribution.MultivariateNormalDistributionModel;
import dr.inference.model.MatrixParameter;
import dr.inference.model.Parameter;
import dr.inference.operators.GibbsOperator;
import dr.inference.operators.OperatorFailedException;
import dr.inference.operators.SimpleMCMCOperator;
import dr.math.distributions.MultivariateNormalDistribution;
import dr.math.matrixAlgebra.IllegalDimension;
import dr.math.matrixAlgebra.Matrix;
import dr.math.matrixAlgebra.Vector;
import dr.util.Attribute;

import java.util.List;

/**
 * Created with IntelliJ IDEA.
 * User: max
 * Date: 2/27/14
 * Time: 4:21 PM
 * To change this template use File | Settings | File Templates.
 */
public class MultivariateNormalGibbsOperator extends SimpleMCMCOperator implements GibbsOperator {
    private Matrix priorPrecision;
    private Vector priorMean;
    private MatrixParameter likelihoodPrecision;
    private Parameter likelihoodMean;
    private MultivariateDistributionLikelihood likelihood;
    private int dim;


    MultivariateNormalGibbsOperator(MultivariateDistributionLikelihood likelihood, MultivariateDistributionLikelihood prior, Double weight) throws IllegalDimension {

        
        MultivariateNormalDistribution tempPrior=(MultivariateNormalDistribution) prior.getDistribution();
        this.priorMean=new Vector(tempPrior.getMean());
        this.priorPrecision=new Matrix(tempPrior.getScaleMatrix());
         MultivariateNormalDistributionModel tempLikelihood=(MultivariateNormalDistributionModel) likelihood.getDistribution();
        this.likelihoodMean=tempLikelihood.getMeanParameter();
        this.likelihoodPrecision=tempLikelihood.getPrecisionMatrixParameter();
        this.likelihood=likelihood;
        this.dim=likelihoodMean.getValues().length;
//        if(dataTemp.contains(MatrixParameter.class))
//        {System.err.print("Well, at least you know it's there...\n");}
//        else{System.err.print("Nope, you screwed up\n");}

        setWeight(weight);
    }


    private void setParameterValue(Parameter set, double[] value){
        set.setDimension(value.length);
        for(int i=0; i<value.length; i++)
        {set.setParameterValueQuietly(i,value[i]);}
        set.fireParameterChangedEvent();
    }

    private double[] getMeanSum(){
        double[] answer=new double[dim];
        List<Attribute<double[]>> dataList = likelihood.getDataList();
        for(Attribute<double[]> d: dataList){
            for(int i=0; i<d.getAttributeValue().length; i++)
            {
                answer[i]+=d.getAttributeValue()[i];
            }
        }
//        System.err.print(answer[0]);
//        System.err.print("\n");
        return answer;}

    private Matrix getPrecision() throws IllegalDimension {
        Matrix currentPrecision=new Matrix(likelihoodPrecision.getParameterAsMatrix());
        return priorPrecision.add(currentPrecision).inverse();
    }

    private Vector getMean() throws IllegalDimension {
        Vector meanSum=new Vector(getMeanSum());
        Matrix workingPrecision=new Matrix(likelihoodPrecision.getParameterAsMatrix());
        Vector meanPart=workingPrecision.product(meanSum);
        meanPart=meanPart.add(priorPrecision.product(priorMean));
        Matrix precisionPart=getPrecision();

        Vector answer=precisionPart.product(meanPart);
//        this.priorPrecision=new Matrix(prior.getDistribution().getScaleMatrix());

//        System.out.print(answer.toComponents()[0]);
        return answer;
    }

//    private Vector getDraws() throws IllegalDimension{
//        double[] rUniform=new double[dim];
//        for(int i=0; i<dim; i++)
//        {rUniform[i]=}
//        Vector draws=new Vector(MultivariateNormalDistribution.);
//        return draws;
//    }

    @Override
    public String getOperatorName() {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public String getPerformanceSuggestion() {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public double doOperation() throws OperatorFailedException{
        double[] draws=null;
        try {
            draws=MultivariateNormalDistribution.nextMultivariateNormalPrecision(getMean().toComponents(), getPrecision().toComponents());
        } catch (IllegalDimension illegalDimension) {
            illegalDimension.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }

        setParameterValue(likelihoodMean, draws);

        return 0;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public int getStepCount() {
        return 1;  //To change body of implemented methods use File | Settings | File Templates.
    }
}
