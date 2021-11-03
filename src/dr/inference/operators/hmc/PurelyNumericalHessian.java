package dr.inference.operators.hmc;

import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.hmc.HessianWrtParameterProvider;
import dr.inference.model.HessianProvider;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.math.MachineAccuracy;
import dr.xml.Reportable;

/**
 * @author Marc A. Suchard
 * @author Xiang Ji
 * @author Andy Magee
 */
public class PurelyNumericalHessian implements HessianWrtParameterProvider, Reportable {

    private final GradientWrtParameterProvider gradientProvider;

    public PurelyNumericalHessian(GradientWrtParameterProvider gradientWrtParameterProvider) {
        this.gradientProvider = gradientWrtParameterProvider;
    }

    @Override
    public double[] getDiagonalHessianLogDensity() {

        final int dim = gradientProvider.getDimension();
        double[][] numericalHessian = getNumericalHessianCentral(); //todo: no need to get the full hessian if only need the diagonals
        double[] numericalHessianDiag = new double[dim];

        for (int i = 0; i < dim; i++) {
            numericalHessianDiag[i] = numericalHessian[i][i];
        }
        return numericalHessianDiag;
    }

    @Override
    public double[][] getHessianLogDensity() {
        return getNumericalHessianCentral();
    }

    @Override
    public Likelihood getLikelihood() {
        return gradientProvider.getLikelihood();
    }

    @Override
    public Parameter getParameter() {
        return gradientProvider.getParameter();
    }

    @Override
    public int getDimension() {
        return gradientProvider.getDimension();
    }

    @Override
    public double[] getGradientLogDensity() {
        final int dim = gradientProvider.getDimension();
        double[] gradient = new double[dim];

        final double[] oldValues = gradientProvider.getParameter().getParameterValues();

        double[] likelihoodPlus = new double[dim];
        double[] likelihoodMinus = new double[dim];

        double[] h = new double[dim];
        for (int i = 0; i < dim; i++) {
            h[i] = MachineAccuracy.SQRT_SQRT_EPSILON * (Math.abs(oldValues[i]) + 1.0);
            gradientProvider.getParameter().setParameterValue(i, oldValues[i] + h[i]);
            likelihoodPlus[i] = gradientProvider.getLikelihood().getLogLikelihood();

            gradientProvider.getParameter().setParameterValue(i, oldValues[i] - h[i]);
            likelihoodMinus[i] = gradientProvider.getLikelihood().getLogLikelihood();
            gradientProvider.getParameter().setParameterValue(i, oldValues[i]);
        }

        for (int i = 0; i < dim; i++) {
            gradient[i] = (likelihoodPlus[i] - likelihoodMinus[i]) / (2.0 * h[i]);
        }

        return gradient;
    }

    private double[][] getNumericalHessianCentral() {

        final int dim = gradientProvider.getDimension();
        double[][] hessian = new double[dim][dim];

        final double[] oldValues = gradientProvider.getParameter().getParameterValues();

        double[][] gradientPlus = new double[dim][dim];
        double[][] gradientMinus = new double[dim][dim];

        double[] h = new double[dim];
        for (int i = 0; i < dim; i++) {
            h[i] = MachineAccuracy.SQRT_SQRT_EPSILON * (Math.abs(oldValues[i]) + 1.0);
            gradientProvider.getParameter().setParameterValue(i, oldValues[i] + h[i]);
            gradientPlus[i] = gradientProvider.getGradientLogDensity();

            gradientProvider.getParameter().setParameterValue(i, oldValues[i] - h[i]);
            gradientMinus[i] = gradientProvider.getGradientLogDensity();
            gradientProvider.getParameter().setParameterValue(i, oldValues[i]);
        }

        for (int i = 0; i < dim; i++) {
            for (int j = i; j < dim; j++) {
                hessian[j][i] = hessian[i][j] = (gradientPlus[j][i] - gradientMinus[j][i]) / (4.0 * h[j]) + (gradientPlus[i][j] - gradientMinus[i][j]) / (4.0 * h[i]);
            }
        }

        return hessian;
    }

    @Override
    public String getReport() {
        return GradientWrtParameterProvider.getReportAndCheckForError(this,
                Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY,
                GradientWrtParameterProvider.TOLERANCE);
    }
}

