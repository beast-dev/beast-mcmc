package dr.inference.distribution;

import dr.inference.model.*;
import dr.math.MathUtils;
import dr.math.distributions.*;
import dr.inferencexml.distribution.PointMassMixtureDistributionModelParser;

public class PointMassMixtureDistributionModel extends AbstractModel implements ParametricMultivariateDistributionModel {

    public PointMassMixtureDistributionModel(Parameter weights, CompoundParameter realizedParameters, boolean weightsNormalized) {
        super(PointMassMixtureDistributionModelParser.POINT_MASS_MIXTURE_DISTRIBUTION_MODEL);
        this.weights = weights;
        addVariable(this.weights);
        this.realizedParameters = realizedParameters;
        addVariable(this.realizedParameters);
        this.weightsNormalized = weightsNormalized;
        this.dim = realizedParameters.getParameter(0).getDimension();
    }

    // *****************************************************************
    // Interface MultivariateDistribution
    // *****************************************************************


    public double logPdf(double[] x) {
        return PointMassMixtureDistribution.logPdf(x,getProbabilities(),getRealizedValues());
    }

    public double[][] getScaleMatrix() {
        throw new RuntimeException("getScaleMatrix should not get called");
        //return PointMassMixtureDistribution.getScaleMatrix(getProbabilities(),getRealizedValues());
    }

    public double[] getMean() {
        return PointMassMixtureDistribution.getMean(getProbabilities(),getRealizedValues());
    }

    public String getType() {
        return PointMassMixtureDistribution.TYPE;
    }

    // *****************************************************************
    // Interface Model
    // *****************************************************************

    public void handleModelChangedEvent(Model model, Object object, int index) {
        // no intermediates need to be recalculated...
    }

    public Likelihood getLikelihood() {
        return null;
    }

    protected final void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
    }

    protected void storeState() {
    }

    protected void restoreState() {
    }

    protected void acceptState() {
    }

    //@Override
    //public int getDimension() {
    //return weights.getDimension();
    //}

    //@Override
    //public double[][] getPrecisionMatrix() {
    //    return PointMassMixtureDistribution.getCovarianceMatrix(getProbabilities(),getRealizedValues());
    //}

    // *****************************************************************
    // Interface DensityModel
    // *****************************************************************

    @Override
    public Variable<Double> getLocationVariable() {
        double[] loc = PointMassMixtureDistribution.getMean(getProbabilities(),getRealizedValues());
        Parameter locParam = new Parameter.Default(loc.length);
        for(int i = 0; i < locParam.getSize(); i++){
            locParam.setParameterValue(i,loc[i]);
        }
        return locParam;
    }


    // **************************************************************
    // Private instance variables and functions
    // **************************************************************

    public double[] getProbabilities(){
        if(weightsNormalized){
            return weights.getParameterValues();
        }else{
            double[] w = weights.getParameterValues();
            double normConst = 0;
            for(int i = 0; i < w.length; i++){
                normConst = normConst + w[i];
            }
            for(int j = 0; j < w.length; j++){
                w[j] = w[j]/normConst;
            }
            return w;
        }
    }


    private double[][] getRealizedValues(){
        int numRows = realizedParameters.getParameter(0).getSize();
        int numCols = weights.getSize();
        double[][] mat = new double[numRows][numCols];
        for(int i = 0; i < numRows; i++){
            for(int j = 0; j < numCols; j++){
                mat[i][j] = realizedParameters.getParameterValue(i,j);
            }
        }
        return mat;
    }

    private double[] getWeights(){
        return weights.getParameterValues();
    }


    private Parameter weights;
    private CompoundParameter realizedParameters;
    private boolean weightsNormalized;
    private int dim;

    // RandomGenerator interface
    public double[] nextRandom() {
        int index = MathUtils.randomChoicePDF(getProbabilities());
        return realizedParameters.getParameter(index).getParameterValues();
    }

    public double logPdf(Object x) {
        double[] v = (double[]) x;
        return PointMassMixtureDistribution.logPdf(v,getProbabilities(),getRealizedValues());
    }

    public int getDimension(){
        return dim;
    }

}