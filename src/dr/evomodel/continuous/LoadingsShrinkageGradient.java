package dr.evomodel.continuous;

import dr.inference.distribution.shrinkage.BayesianBridgeLikelihood;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.*;

import java.util.ArrayList;


/**
 * @author Gabriel Hassler
 * @author Xiang Ji
 * @author Marc A Suchard
 */

public class LoadingsShrinkageGradient implements GradientWrtParameterProvider {

    private final MatrixParameterInterface loadings;
    private final BayesianBridgeLikelihood[] rowPriors;
    private final CompoundLikelihood likelihood;

    public LoadingsShrinkageGradient(String name,
                                     MatrixParameterInterface loadings,
                                     BayesianBridgeLikelihood[] rowPriors) {

        this.loadings = loadings;

        this.rowPriors = rowPriors;

        ArrayList<Likelihood> likelihoods = new ArrayList<>(rowPriors.length);
        for (int i = 0; i < rowPriors.length; i++) {
            likelihoods.add(i, rowPriors[i]);
        }

        this.likelihood = new CompoundLikelihood(likelihoods); //TODO: just construct using compound likelihood??

    }

    @Override
    public Likelihood getLikelihood() {
        return likelihood;
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
        double[] gradientLogDensity = new double[loadings.getDimension()];
        for (int i = 0; i < rowPriors.length; i++) {
            double[] grad = rowPriors[i].getGradientLogDensity();

            int offset = i * loadings.getRowDimension();

            for (int j = 0; j < loadings.getRowDimension(); j++) {
                gradientLogDensity[j + offset] = grad[j];
            }
        }

        return gradientLogDensity;
    }

}

