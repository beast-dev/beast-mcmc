package dr.inference.distribution;


import dr.inference.model.*;
import dr.inference.model.Parameter;
import dr.inferencexml.distribution.MomentDistributionModelParser;
import dr.math.distributions.RandomGenerator;

//@author Max Tolkoff
public class MomentDistributionModel extends AbstractModelLikelihood implements ParametricMultivariateDistributionModel, RandomGenerator {

    public MomentDistributionModel(Parameter mean, Parameter precision, Parameter cutoff, Parameter data) {
        super(MomentDistributionModelParser.MOMENT_DISTRIBUTION_MODEL);

        this.mean=mean;
        this.precision=precision;
//        this.mean = new DuplicatedParameter(mean);
//        this.mean.addDuplicationParameter(new Parameter.Default(cutoff.getDimension()));
//        DuplicatedParameter precTemp= new DuplicatedParameter(precision);
//        precTemp.addDuplicationParameter(new Parameter.Default(cutoff.getDimension()));
//        this.precision=new DiagonalMatrix(precTemp);
        addVariable(mean);
        mean.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY, 1));
        addVariable(precision);
//        precision.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));
        this.cutoff=cutoff;
        addVariable(cutoff);
        cutoff.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, cutoff.getDimension()));
        addVariable(data);
        this.data=data;
        untruncated=new NormalDistributionModel(mean, precision, true);
        sumKnown=false;
        untruncatedKnown=false;
    }

    private final Parameter mean;
    private final Parameter precision;
//    private final DuplicatedParameter mean;
//    private final DiagonalMatrix precision;
    private final Parameter cutoff;
    private NormalDistributionModel untruncated;
    private double sum;
    private boolean sumKnown;

    private boolean storedSumKnown;
    private double storedSum;
    private boolean untruncatedKnown;
    private boolean storedUntruncatedKnown;
    private NormalDistributionModel storedUntruncated;
    private Parameter data;

    public double logPdf(Parameter data) {
//        untruncatedKnown=false;
//        sumKnown=false;
        checkDistribution();
        if(sumKnown)
            return sum;
        else
        {
            sum=0;
        }
        if(data.getDimension()!=cutoff.getDimension()){
            throw new RuntimeException("Incorrect number of cutoffs");
        }
        for (int i = 0; i <data.getDimension() ; i++) {
            if (Math.sqrt(precision.getParameterValue(0) * cutoff.getParameterValue(i)) > Math.abs(data.getParameterValue(i)) && data.getParameterValue(i)!=0)
                return Double.NEGATIVE_INFINITY;
            else if(data.getParameterValue(i)==0)
                sum+=-1000-Math.log(precision.getParameterValue(0));
            else
                sum+=untruncated.logPdf(data.getParameterValue(i));//(2*untruncated.logPdf(cutoff.getParameterValue(i)));
        }
        sumKnown=true;
        return sum;
        //}

    }

    @Override
    public double logPdf(double[] x) {
        return 0;
    }

    @Override
    public double[][] getScaleMatrix() {
        double[][] temp=new double[1][1];
        temp[0][0]=precision.getParameterValue(0);
        return temp;
//        return precision.getParameterAsMatrix();
    }

    @Override
    public double[] getMean() {
        return mean.getParameterValues();
    }

    @Override
    public String getType() {
        return "Moment Distribution Model";
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {

    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        sumKnown=false;
        if(variable==mean || variable==precision)
        {untruncatedKnown=false;}

    }

    @Override
    protected void storeState() {
        storedSumKnown=sumKnown;
        storedSum=sum;
        storedUntruncated=untruncated;
        storedUntruncatedKnown=untruncatedKnown;
    }

    @Override
    protected void restoreState() {
        sumKnown=storedSumKnown;
        sum=storedSum;
        untruncated=storedUntruncated;
        untruncatedKnown=storedUntruncatedKnown;

    }

    @Override
    protected void acceptState() {

    }

    private NormalDistributionModel createNewDistribution() {
        return new NormalDistributionModel(mean, precision, true);
//        return new NormalDistributionModel(new Parameter.Default(mean.getParameterValue(0)), new Parameter.Default(precision.getParameterValue(0)), true);
    }

    private void checkDistribution() {
        if (!untruncatedKnown) {
            untruncated = createNewDistribution();
            untruncatedKnown = true;
        }
    }

    @Override
    public double[] nextRandom() {
        return new double[0];
    }

    @Override
    public double logPdf(Object x) {
        if(x instanceof Parameter)
         return logPdf((Parameter) x);
        else
            return 0;
    }

    @Override
    public Model getModel() {
        return this;
    }

    @Override
    public double getLogLikelihood() {
        return logPdf(data);
    }

    @Override
    public void makeDirty() {
        sumKnown=false;
        untruncatedKnown=false;

    }
}
