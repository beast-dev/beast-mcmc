package dr.inference.distribution;

import com.sun.tools.doclets.internal.toolkit.util.SourceToHTMLConverter;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.inferencexml.distribution.MomentDistributionModelParser;
import dr.math.UnivariateFunction;
import dr.math.distributions.MultivariateNormalDistribution;
import dr.math.distributions.RandomGenerator;

/**
 * Created by max on 5/13/15.
 */
public class MomentDistributionModel extends AbstractModel implements ParametricMultivariateDistributionModel, RandomGenerator {

    public MomentDistributionModel(Parameter mean, Parameter precision, Parameter cutoff) {
        super(MomentDistributionModelParser.MOMENT_DISTRIBUTION_MODEL);

        this.mean = mean;
        this.precision = precision;
        addVariable(mean);
        mean.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY, 1));
        addVariable(precision);
        precision.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));
        this.cutoff=cutoff;
        addVariable(cutoff);
       cutoff.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));
        untruncated=new NormalDistributionModel(mean, precision, true);
        sumKnown=false;
        untruncatedKnown=false;
    }

    private final Parameter mean;
    private final Parameter precision;
    private final Parameter cutoff;
    private NormalDistributionModel untruncated;
    private double sum;
    private boolean sumKnown;

    private boolean storedSumKnown;
    private double storedSum;
    private boolean untruncatedKnown;
    private boolean storedUntruncatedKnown;
    private NormalDistributionModel storedUntruncated;

    @Override
    public double logPdf(double[] x) {
        System.out.println(x[0]);
        checkDistribution();
//        if(sumKnown)
//            return sum;
//        else
//        {
            sum=0;
//        }
        if(x.length!=cutoff.getDimension()){
            throw new RuntimeException("Incorrect number of cutoffs");
        }
        for (int i = 0; i <x.length ; i++) {
            if (Math.sqrt(precision.getParameterValue(0) * cutoff.getParameterValue(i)) > Math.abs(x[i]) && x[i]!=0)
                return Double.NEGATIVE_INFINITY;
            else if(x[i]==0)
                sum+=-1000-Math.log(precision.getParameterValue(0));
            else
                sum+=untruncated.logPdf(x[i]);//(2*untruncated.logPdf(cutoff.getParameterValue(i)));
        }
//            sumKnown=true;
        return sum;
        //}

    }

    @Override
    public double[][] getScaleMatrix() {
        double[][] temp=new double[1][1];
        temp[0][0]=precision.getParameterValue(0);
        return temp;
    }

    @Override
    public double[] getMean() {
        double[] meanList=new double[1];
        meanList[0]=mean.getParameterValue(0);
        return meanList;
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
        if(variable==mean)
        {untruncatedKnown=false;}

    }

    @Override
    protected void storeState() {
//        storedSumKnown=sumKnown;
//        storedSum=sum;
        storedUntruncated=untruncated;
        storedUntruncatedKnown=untruncatedKnown;
    }

    @Override
    protected void restoreState() {
//        sumKnown=storedSumKnown;
//        sum=storedSum;
        untruncated=storedUntruncated;
        untruncatedKnown=storedUntruncatedKnown;

    }

    @Override
    protected void acceptState() {

    }

    private NormalDistributionModel createNewDistribution() {
        return new NormalDistributionModel(mean, precision, true);
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
         return logPdf(((Parameter) x).getParameterValues());
        else
            return 0;
    }
}
