package dr.evomodel.continuous;

import dr.inference.distribution.NormalStatisticsProvider;
import dr.inference.distribution.shrinkage.BayesianBridgeLikelihood;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.*;
import dr.inference.operators.repeatedMeasures.MultiplicativeGammaGibbsHelper;
import dr.xml.Reportable;

import java.util.Arrays;

/**
 * @author Gabriel Hassler
 */

public class MatrixShrinkageLikelihood extends AbstractModelLikelihood implements GradientWrtParameterProvider,
        NormalStatisticsProvider, Reportable, MultiplicativeGammaGibbsHelper {

    private final MatrixParameterInterface loadings;
    private final BayesianBridgeLikelihood[] rowPriors;

    private final double[] gradientLogDensity;

    private final CompoundLikelihood likelihood;

    public MatrixShrinkageLikelihood(String name,
                                     MatrixParameterInterface loadings,
                                     BayesianBridgeLikelihood[] rowPriors) {
        super(name);

        this.loadings = loadings;

        this.rowPriors = rowPriors;

        for (BayesianBridgeLikelihood model : rowPriors) {
            addModel(model);
        }

        this.gradientLogDensity = new double[loadings.getDimension()];
        addVariable(loadings);

        this.likelihood = new CompoundLikelihood(Arrays.asList(rowPriors));


    }

    @Override
    public Likelihood getLikelihood() {
        return this;
    }

    public BayesianBridgeLikelihood getLikelihood(int i) {
        return rowPriors[i];
    }

    @Override
    public Parameter getParameter() {
        return loadings;
    }

    @Override
    public int getDimension() {
        return loadings.getDimension();
    }

    @Override
    public double[] getGradientLogDensity() {
        for (int i = 0; i < rowPriors.length; i++) {
            double[] grad = rowPriors[i].getGradientLogDensity();

            int offset = i * loadings.getRowDimension();

            for (int j = 0; j < loadings.getRowDimension(); j++) {
                gradientLogDensity[j + offset] = grad[j];
            }
        }

        return gradientLogDensity;
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
        for (BayesianBridgeLikelihood likelihood : rowPriors) {
            likelihood.handleModelChangedEvent(model, object, index);
        }
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        for (BayesianBridgeLikelihood likelihood : rowPriors) {
            likelihood.handleVariableChangedEvent(variable, index, type);
        }
    }

    @Override
    protected void storeState() {
        for (BayesianBridgeLikelihood likelihood : rowPriors) {
            likelihood.storeState();
        }
    }

    @Override
    protected void restoreState() {
        for (BayesianBridgeLikelihood likelihood : rowPriors) {
            likelihood.restoreState();
        }
    }

    @Override
    protected void acceptState() {
        for (BayesianBridgeLikelihood likelihood : rowPriors) {
            likelihood.acceptState();
        }
    }

    @Override
    public Model getModel() {
        return this;
    }

    @Override
    public double getLogLikelihood() {
        return likelihood.getLogLikelihood();
    }

    @Override
    public void makeDirty() {
        likelihood.makeDirty();
    }

    @Override
    public double getNormalMean(int dim) {
        return 0; //BayesianBridgeDistribution assumes zero mean
    }

    @Override
    public double getNormalSD(int dim) {
        int row = dim / loadings.getRowDimension();
        int col = dim - row * loadings.getRowDimension();
        double globalScale = rowPriors[row].getGlobalScale().getParameterValue(0);
        double localScale = rowPriors[row].getLocalScale().getParameterValue(col);
        return globalScale * localScale;
    }

    @Override
    public String getReport() {
        StringBuilder sb = new StringBuilder("MatrixShrinkageLikelihood\n");
        int counter = 0;
        for (BayesianBridgeLikelihood like : rowPriors) {
            counter += 1;
            sb.append("\tLikelihood " + counter + ": ");
            sb.append(like.getLogLikelihood());
            sb.append("\n");
        }
        sb.append("Likelihood: ");
        sb.append(getLogLikelihood());
        sb.append("\n");
        return sb.toString();
    }

    @Override
    public double computeSumSquaredErrors(int column) {
        int p = loadings.getRowDimension();

        double sumSquares = 0;
        for (int j = 0; j < p; j++) {
            double localSD = getLikelihood(column).getLocalScale().getParameterValue(j);
            double loadEl = loadings.getParameterValue(j, column);
            double x = loadEl / localSD;
            sumSquares += x * x;
        }
        return sumSquares;
    }

    @Override
    public int getRowDimension() {
        return loadings.getRowDimension();
    }

    @Override
    public int getColumnDimension() {
        return loadings.getColumnDimension();
    }
}