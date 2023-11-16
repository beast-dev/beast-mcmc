package dr.inference.distribution;

import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.*;
import dr.inference.operators.repeatedMeasures.MultiplicativeGammaGibbsHelper;
import dr.math.distributions.NormalDistribution;
import dr.xml.Reportable;

/**
 * @author Max Tolkoff
 * @author Gabriel Hassler
 * @author Marc Suchard
 */

public class IndependentNormalDistributionModel extends AbstractModelLikelihood implements NormalStatisticsProvider,
        MultiplicativeGammaGibbsHelper, GradientWrtParameterProvider, Reportable {
    Parameter mean;
    Parameter variance;
    Parameter precision;
    Parameter data;
    boolean usePrecision;

    public static String INDEPENDENT_NORMAL_DISTRIBUTION_MODEL = "independentNormalDistributionModel";


    public IndependentNormalDistributionModel(String id, Parameter mean, Parameter variance, Parameter precision, Parameter data) {
        super(id);
        addVariable(mean);
        this.mean = mean;
        if (precision != null) {
            usePrecision = true;
            addVariable(precision);
        } else {
            usePrecision = false;
            addVariable(variance);
        }
        this.precision = precision;
        this.variance = variance;
        this.data = data;
        addVariable(data);
    }


    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {

    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {

    }

    @Override
    protected void storeState() {

    }

    @Override
    protected void restoreState() {

    }

    @Override
    protected void acceptState() {

    }

    @Override
    public Model getModel() {
        return this;
    }

    @Override
    public double getLogLikelihood() {
        double sum = 0;
        for (int i = 0; i < data.getDimension(); i++) {
            double sd = getNormalSD(i);
            sum += NormalDistribution.logPdf(data.getParameterValue(i),
                    mean.getParameterValue(i),
                    sd);
        }
        return sum;
    }

    @Override
    public void makeDirty() {

    }

    public boolean precisionUsed() {
        return usePrecision;
    }

    public Parameter getPrecision() {
        return precision;
    }

    public Parameter getVariance() {
        return variance;
    }

    public Parameter getMean() {
        return mean;
    }

    public Parameter getData() {
        return data;
    }

    @Override
    public double getNormalMean(int dim) {
        return mean.getParameterValue(dim);
    }

    @Override
    public double getNormalSD(int dim) {
        if (usePrecision) {
            return 1 / Math.sqrt(precision.getParameterValue(dim));
        } else {
            return Math.sqrt(variance.getParameterValue(dim));
        }
    }

    @Override
    public double computeSumSquaredErrors(int column) {
        double error = mean.getParameterValue(column) - data.getParameterValue(column);
        return error * error;
    }

    @Override
    public int getRowDimension() {
        return 1;
    }

    @Override
    public int getColumnDimension() {
        return data.getDimension();
    }

    @Override
    public Likelihood getLikelihood() {
        return this;
    }

    @Override
    public Parameter getParameter() {
        return data;
    }

    @Override
    public int getDimension() {
        return data.getDimension();
    }

    @Override
    public double[] getGradientLogDensity() {
        double[] grad = new double[getDimension()];
        for (int i = 0; i < getDimension(); i++) {
            double sd = getNormalSD(i);
            grad[i] = NormalDistribution.gradLogPdf(data.getParameterValue(i), mean.getParameterValue(i), sd);
        }
        return grad;
    }

    @Override
    public String getReport() {
        StringBuilder sb = new StringBuilder(INDEPENDENT_NORMAL_DISTRIBUTION_MODEL + " report:\n");
        sb.append("\tlogLikelihood: " + getLogLikelihood() + "\n");
        sb.append("\tgradient: ");
        double[] grad = getGradientLogDensity();
        for (int i = 0; i < grad.length; i++) {
            sb.append(grad[i] + " ");
        }
        sb.append("\n\n");
        return sb.toString();
    }
}
