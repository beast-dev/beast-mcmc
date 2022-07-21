package dr.evomodel.treedatalikelihood.discrete;

import com.github.lbfgs4j.liblbfgs.Function;
import dr.inference.hmc.HessianWrtParameterProvider;
import dr.inference.model.Likelihood;
import dr.inference.operators.hmc.NumericalHessianFromGradient;
import dr.math.distributions.MultivariateNormalDistribution;
import dr.math.matrixAlgebra.Matrix;
import dr.util.Timer;
import dr.util.Transform;
import dr.xml.Reportable;

import java.util.Arrays;

/**
 * @author Andy Magee
 */

public class LaplaceApproximation implements Reportable {

    private final MaximizerWrtParameter maximizer;
    private final boolean diagonal;
    private final boolean estimateKL;
    private final boolean computeML;
    private double[] sdVector;
    private double[][] covarianceMatrix;
    private double[][] precisionMatrix;
    private double time = 0.0;
    private double KL = Double.NaN;
    private double ML = Double.NaN;

    public LaplaceApproximation(MaximizerWrtParameter maximizer, boolean diagonal, boolean KL, boolean ML) {
        this.maximizer = maximizer;
        this.diagonal = diagonal;
        this.estimateKL = KL;
        this.computeML = ML;
    }

    public void approximate() {
        Timer timer = new Timer();
        double[] mu = maximizer.getMinimumPoint(false);
        HessianWrtParameterProvider hessian = new NumericalHessianFromGradient(maximizer.getGradient());

        timer.start();
        precisionMatrix = hessian.getHessianLogDensity();
        for (int i = 0; i < mu.length; i++) {
            for (int j = 0; j < mu.length; j++) {
                precisionMatrix[i][j] = -precisionMatrix[i][j];
            }
        }
        timer.stop();
        time += timer.toSeconds();
        covarianceMatrix = (new Matrix(precisionMatrix).inverse()).toComponents();


        if ( diagonal ) {
            sdVector = new double[mu.length];
            for (int i = 0; i < mu.length; i++) {
                sdVector[i] = Math.sqrt(covarianceMatrix[i][i]);
            }
        }

        MultivariateNormalDistribution distribution = new MultivariateNormalDistribution(mu, precisionMatrix);
        Function function = maximizer.getFunction();
        if ( estimateKL ) {
            int nSamples = 1000;
            KL = 0.0;
            for (int i = 0; i < nSamples; i++) {
                double[] x = distribution.nextMultivariateNormal();
                double lnPosterior = -function.valueAt(x);
                KL += distribution.logPdf(x) - lnPosterior;
            }
            KL /= nSamples;
        }

        if ( computeML ) {
            double logDet = distribution.getLogDet();
            // Kass and Raftery, eqn 4
            ML = mu.length/2.0 * Math.log(2.0 * Math.PI) + 0.5 * logDet - function.valueAt(mu);
        }

    }

    @Override
    public String getReport() {
        StringBuilder sb = new StringBuilder();

        if ( !maximizer.wasExecuted() ) {
            sb.append("Not yet executed.");
        } else {
            sb.append("NB: Laplace approximation is taken in transformed space.\n");
            sb.append("Mean vector: ").append(new dr.math.matrixAlgebra.Vector(maximizer.getMinimumPoint(false))).append("\n");
            if ( diagonal ) {
                sb.append("Stdev vector: ").append(new dr.math.matrixAlgebra.Vector(sdVector)).append("\n");
            } else {
                sb.append("Covariance matrix: ").append(new dr.math.matrixAlgebra.Matrix(covarianceMatrix));
            }
            if (estimateKL) {
                sb.append("KL(approximation||posterior) = ").append(KL).append("\n");
            }
            if (computeML) {
                sb.append("log(marginal likelihood) = ").append(ML).append("\n");
            }
            sb.append("Time to compute Hessian: ").append(time).append("s\n");
        }

        return sb.toString();

    }
}
