package dr.evomodel.continuous;

import dr.evolution.tree.NodeRef;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.CompoundParameter;
import dr.inference.model.Model;
import dr.math.distributions.MultivariateNormalDistribution;
import dr.math.matrixAlgebra.IllegalDimension;
import dr.math.matrixAlgebra.Matrix;

import java.util.List;

/**
 * Integrated multivariate trait likelihood that assumes a semi-conjugate prior on the root.
 * The semi-conjugate prior is a multivariate normal distribution with an independent precision
 *
 * @author Marc A. Suchard
 */
public class SemiConjugateMultivariateTraitLikelihood extends IntegratedMultivariateTraitLikelihood {

    public SemiConjugateMultivariateTraitLikelihood(String traitName,
                                                    TreeModel treeModel,
                                                    MultivariateDiffusionModel diffusionModel,
                                                    CompoundParameter traitParameter,
                                                    List<Integer> missingIndices,
                                                    boolean cacheBranches,
                                                    boolean scaleByTime,
                                                    boolean useTreeLength,
                                                    BranchRateModel rateModel,
                                                    Model samplingDensity,
                                                    boolean reportAsMultivariate,
                                                    MultivariateNormalDistribution rootPrior,
                                                    boolean reciprocalRates) {

        super(traitName, treeModel, diffusionModel, traitParameter, missingIndices, cacheBranches, scaleByTime,
                useTreeLength, rateModel, samplingDensity, reportAsMultivariate, reciprocalRates);

        setRootPrior(rootPrior); // Semi-conjugate multivariate normal with own mean and precision
    }

    @Override
    public boolean getComputeWishartSufficientStatistics() {
        return false;  // No need for outer products, as Gibbs sampling of diffusion matrix is not possible
    }

    protected  double calculateAscertainmentCorrection(int taxonIndex) {
        throw new RuntimeException("Ascertainment correction not yet implemented for semi-conjugate trait likelihoods");
    }

    protected double getRescaledLengthToRoot(NodeRef node) {
        double length = 0;
        final NodeRef root = treeModel.getRoot();
        while (node != root) {
            length += getRescaledBranchLength(node);
            node = treeModel.getParent(node);
        }
        return length;
    }

    protected double integrateLogLikelihoodAtRoot(double[] y,
                                                  double[] Ay,
                                                  double[][] AplusB,
                                                  double[][] treePrecision, double rootPrecision) {
        double detAplusB = 0;
        double square = 0;

        // square : (Ay + Bz)' (A+B)^{-1} (Ay + Bz)

        if (dimTrait > 1) {
            for (int i = 0; i < dimTrait; i++) {
                Ay[i] += Bz[i];   // Ay is filled with sum, and original value is destroyed
                for (int j = 0; j < dimTrait; j++) {
                    AplusB[i][j] = treePrecision[i][j] * rootPrecision + rootPriorPrecision[i][j];
                }
            }

            Matrix mat = new Matrix(AplusB);

            try {
                detAplusB = mat.determinant();

            } catch (IllegalDimension illegalDimension) {
                illegalDimension.printStackTrace();
            }

            double[][] invAplusB = mat.inverse().toComponents();

            for (int i = 0; i < dimTrait; i++) {
                for (int j = 0; j < dimTrait; j++)
                    square += Ay[i] * invAplusB[i][j] * Ay[j];
            }
        } else {
            // 1D is very simple
            detAplusB = treePrecision[0][0] * rootPrecision + rootPriorPrecision[0][0];
            Ay[0] += Bz[0];
            square = Ay[0] * Ay[0] / detAplusB;
        }

        double retValue = 0.5 * (logRootPriorPrecisionDeterminant - Math.log(detAplusB) - zBz + square);

        if (DEBUG) {
            System.err.println("(Ay+Bz)(A+B)^{-1}(Ay+Bz) = " + square);
            System.err.println("density = " + retValue);
            System.err.println("zBz = " + zBz);
        }

        return retValue;
    }

    private void setRootPriorSumOfSquares() {
        if (integrateRoot) {
            Bz = new double[dimTrait];
            // z'Bz -- sum-of-squares root contribution
            zBz = computeWeightedAverageAndSumOfSquares(rootPriorMean, Bz, rootPriorPrecision, dimTrait, 1.0);
        } else {
            zBz = 0;
        }
    }

    private void setRootPrior(MultivariateNormalDistribution rootPrior) {
        rootPriorMean = rootPrior.getMean();
        rootPriorPrecision = rootPrior.getScaleMatrix();

        try {
            logRootPriorPrecisionDeterminant = Math.log(new Matrix(rootPriorPrecision).determinant());
        } catch (IllegalDimension illegalDimension) {
            illegalDimension.printStackTrace();
        }
        setRootPriorSumOfSquares();
    }

    protected double[][] computeMarginalRootMeanAndVariance(double[] rootMean, double[][] treePrecision,
                                                            double[][] treeVariance, double rootPrecision) {

        computeWeightedAverageAndSumOfSquares(rootMean, Ay, treePrecision, dimTrait, rootPrecision); // Fills in Ay

        double[][] AplusB = tmpM;

        for (int i = 0; i < dimTrait; i++) {
            Ay[i] += Bz[i];   // Ay is filled with sum, and original value is destroyed
            for (int j = 0; j < dimTrait; j++) {
                AplusB[i][j] = treePrecision[i][j] * rootPrecision + rootPriorPrecision[i][j];
            }
        }
        Matrix mat = new Matrix(AplusB);
        double[][] invAplusB = mat.inverse().toComponents();

        // Expected value: (A + B)^{-1}(Ay + Bz)
        for (int i = 0; i < dimTrait; i++) {
            rootMean[i] = 0.0;
            for (int j = 0; j < dimTrait; j++) {
                rootMean[i] += invAplusB[i][j] * Ay[j];
            }
        }
        return invAplusB;
    }

    protected double[] rootPriorMean;
    protected double[][] rootPriorPrecision;
    protected double logRootPriorPrecisionDeterminant;
    protected double[] Bz;
    private double zBz; // Prior sum-of-squares contribution
}
