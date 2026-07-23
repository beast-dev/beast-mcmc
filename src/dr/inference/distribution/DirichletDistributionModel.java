package dr.inference.distribution;

import dr.inference.model.*;
import dr.inferencexml.distribution.DirichletDistributionModelParser;
import dr.math.GammaFunction;
import dr.math.MathUtils;
import dr.math.distributions.DirichletDistribution;

public class DirichletDistributionModel extends AbstractModel implements ParametricMultivariateDistributionModel {

    private Parameter countsParameter;
    private Parameter dispersionParameter;
    private boolean sumToNumberOfElements;
    public static final double ACCURACY_THRESHOLD = 1E-12;
    public static final boolean DEBUG = false;
    private int dim;

    public DirichletDistributionModel(Parameter countsParam, Parameter dispersionParam, boolean sumToNumberOfElements) {
        super(DirichletDistributionModelParser.DIRICHLET_DISTRIBUTION_MODEL);
        this.sumToNumberOfElements = sumToNumberOfElements;
        this.countsParameter = countsParam;
        addVariable(countsParameter);
        this.dispersionParameter = dispersionParam;
        addVariable(dispersionParameter);
        this.dim = countsParameter.getSize();
    }

    // *****************************************************************
    // Interface MultivariateDistribution
    // *****************************************************************

    public double logPdf(double[] x) {
        if (x.length != countsParameter.getSize()) {
            throw new IllegalArgumentException("data array is of the wrong dimension");
        }

        double returnValue;
        double countsSum = 0;

        for (int i = 0; i < countsParameter.getSize(); i++) {
            //System.err.println("countesParameter.getParameterValue(" + i + "): " + countsParameter.getParameterValue(i));
            //System.err.println("dispersionParameter.getParameterValue(0): " + dispersionParameter.getParameterValue(0));
            countsSum = countsSum + countsParameter.getParameterValue(i)*dispersionParameter.getParameterValue(0);
        }
        returnValue = GammaFunction.lnGamma(countsSum);

        for (int i = 0; i < countsParameter.getSize(); i++) {
            returnValue = returnValue - GammaFunction.lnGamma(countsParameter.getParameterValue(i)*dispersionParameter.getParameterValue(0));
        }

        double countsParameterSum = 0;

        if (sumToNumberOfElements) {
            countsParameterSum = (double)countsParameter.getSize();
        } else {
            countsParameterSum = 1.0;
        }

        returnValue = returnValue - (countsParameter.getSize()) * Math.log(countsParameterSum);

        double parameterSum = 0.0;
        for (int i = 0; i < countsParameter.getSize(); i++) {
            returnValue = returnValue + (countsParameter.getParameterValue(i)*dispersionParameter.getParameterValue(0) - 1) * (Math.log(x[i]) - Math.log(countsParameterSum));
            parameterSum = parameterSum + x[i];
        }

        if (Math.abs(parameterSum - countsParameterSum) > ACCURACY_THRESHOLD) {
            if (DEBUG) {
                System.out.println("Parameters do not sum to " + countsParameterSum);
                for (int i = 0; i < countsParameter.getSize(); i++) {
                    System.out.println("x[" + i + "] = " + x[i]);
                }
                System.out.println("Current sum = " + parameterSum);
            }
            returnValue = Double.NEGATIVE_INFINITY;
        }

        return returnValue;
    }

    public double[][] getScaleMatrix() {
        return null;
    }

    public double[] getMean() {

        double[] returnVal = new double[countsParameter.getSize()];
        double countsSum = 0;
        for (int i = 0; i < countsParameter.getSize(); i++) {
            countsSum = countsSum + countsParameter.getParameterValue(i)*dispersionParameter.getParameterValue(0);
        }
        for(int i = 0; i < countsParameter.getSize(); i++){
            returnVal[i] = (countsParameter.getParameterValue(i)*dispersionParameter.getParameterValue(0))/countsSum;
        }
        return returnVal;
    }

    public int getDimension(){
        return dim;
    }

    public String getType() {
        return DirichletDistribution.TYPE;
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


    // *****************************************************************
    // Interface DensityModel
    // *****************************************************************

    @Override
    public Variable<Double> getLocationVariable() {
        return null;
    }


    // RandomGenerator interface
    public double[] nextRandom() {
        double[] randomGammas = new double[countsParameter.getSize()];
        double sumRandomGammas = 0;
        for(int i = 0; i < countsParameter.getSize(); i++){
            randomGammas[i] = MathUtils.nextGamma(countsParameter.getParameterValue(i)*dispersionParameter.getParameterValue(0),1);
            sumRandomGammas = sumRandomGammas + randomGammas[i];
        }

        double[] returnValue = new double[countsParameter.getSize()];
        for(int i = 0; i < countsParameter.getSize(); i++){
            returnValue[i] = randomGammas[i]/sumRandomGammas;
        }
        return returnValue;
    }

}