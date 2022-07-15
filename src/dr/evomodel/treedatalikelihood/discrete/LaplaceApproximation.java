package dr.evomodel.treedatalikelihood.discrete;

import dr.inference.hmc.HessianWrtParameterProvider;
import dr.inference.model.Likelihood;
import dr.inference.operators.hmc.NumericalHessianFromGradient;
import dr.math.distributions.MultivariateNormalDistribution;
import dr.math.matrixAlgebra.Matrix;
import dr.util.Timer;
import dr.util.Transform;
import dr.xml.Reportable;

/**
 * @author Andy Magee
 */

public class LaplaceApproximation implements Reportable {

    private final MaximizerWrtParameter maximizer;
    private final boolean diagonal;
    private double[] sdVector;
    private double[][] covarianceMatrix;
    private double time = 0.0;

    public LaplaceApproximation(MaximizerWrtParameter maximizer, boolean diagonal) {
        this.maximizer = maximizer;
        this.diagonal = diagonal;
    }

    public void approximate() {
        Timer timer = new Timer();
        double[] mu = maximizer.getMinimumPoint(false);
        HessianWrtParameterProvider hessian = new NumericalHessianFromGradient(maximizer.getGradient());
        if ( diagonal ) {
            timer.start();
            sdVector = hessian.getDiagonalHessianLogDensity();
            timer.stop();
            time += timer.toSeconds();
            for (int i = 0; i < mu.length; i++) {
                sdVector[i] = Math.sqrt(1.0 / (-sdVector[i]));
            }
        } else {
            timer.start();
            covarianceMatrix = hessian.getHessianLogDensity();
            timer.stop();
            time += timer.toSeconds();
            for (int i = 0; i < mu.length; i++) {
                for (int j = 0; j < mu.length; j++) {
                    covarianceMatrix[i][j] = 1.0 / (-covarianceMatrix[i][j]);
                }
            }
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
            sb.append("Time to compute Hessian: ").append(time).append("s\n");
        }

        return sb.toString();

    }
}
