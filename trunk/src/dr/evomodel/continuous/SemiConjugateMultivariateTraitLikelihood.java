package dr.evomodel.continuous;

import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.CompoundParameter;
import dr.inference.model.Model;
import dr.math.distributions.MultivariateNormalDistribution;
import dr.math.matrixAlgebra.IllegalDimension;
import dr.math.matrixAlgebra.Matrix;
import dr.math.matrixAlgebra.SymmetricMatrix;
import dr.math.matrixAlgebra.Vector;

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

    protected double integrateLogLikelihoodAtRoot(double[] Ay, double[] Bz,
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

    protected double integrateLogLikelihoodAtRootFromFullTreeMatrix(double[][] treeTraitPrecisionMatrix,
                                                                    double[] tipTraits) {

        double logLikelihood = 0;
        //final int tipCount = treeModel.getExternalNodeCount();
        final int tipCount = countNonMissingTips();

        // 1^t\Sigma^{-1} y + Pz
        double[] mean = Ay;
        for (int i = 0; i < dimTrait; i++) {
            mean[i] = 0;
            for (int j = 0; j < dimTrait; j++) {
                mean[i] += rootPriorPrecision[i][j] * rootPriorMean[j];
            }

            for (int j = 0; j < tipCount; j++) {
                final int rowOffset = j * dimTrait + i;
                for (int k = 0; k < tipCount * dimTrait; k++) {
                    mean[i] += treeTraitPrecisionMatrix[rowOffset][k] * tipTraits[k];
                }
            }
        }

        // 1^t \Sigma^{-1} 1 + P
        double[][] precision = tmpM;
        for (int i = 0; i < dimTrait; i++) {
            for (int j = 0; j < dimTrait; j++) {
                precision[i][j] = rootPriorPrecision[i][j];
                for (int k = 0; k < tipCount; k++) {
                    for (int l = 0; l < tipCount; l++) {
                        precision[i][j] += treeTraitPrecisionMatrix[k * dimTrait + i][l * dimTrait + j];
                    }
                }
            }
        }
        double[] normalizedMean = tmp2;
        double[][] variance = new SymmetricMatrix(precision).inverse().toComponents();
        for (int i = 0; i < dimTrait; i++) {
            normalizedMean[i] = 0.0;
            for (int j = 0; j < dimTrait; j++) {
                normalizedMean[i] += variance[i][j] * mean[j];
            }
        }
        mean = normalizedMean;

        // We know:  y ~ MVN(x, A) and x ~ N(m, B)
        // Therefore p(x | y) = N( (A+B)^{-1}(Ay + Bm), A + B)
        // We want: p( y ) = p( y | x ) p( x ) / p( x | y ) for any value x, say x = 0

        logLikelihood += MultivariateNormalDistribution.logPdf(
                rootPriorMean, new double[rootPriorMean.length], rootPriorPrecision,
                logRootPriorPrecisionDeterminant, 1.0
        );

        logLikelihood -= MultivariateNormalDistribution.logPdf(
                mean, new double[mean.length], precision,
                Math.log(MultivariateNormalDistribution.calculatePrecisionMatrixDeterminate(precision)), 1.0
        );

        if (DEBUG) {
            System.err.println("Mean = " + new Vector(mean));
            System.err.println("Prec = " + new Matrix(precision));
            System.err.println("log density = " + logLikelihood);
        }
        return logLikelihood;
    }

    private void setRootPriorSumOfSquares() {
        if (integrateRoot) {
            // z'Bz -- sum-of-squares root contribution
            for (int i = 0; i < dimTrait; i++) {
                Bz[i] = 0;
                for (int j = 0; j < dimTrait; j++)
                    Bz[i] += rootPriorPrecision[i][j] * rootPriorMean[j];
                zBz += rootPriorMean[i] * Bz[i];
            }
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

    private double zBz; // Prior sum-of-squares contribution

    protected double[][] computeMarginalRootMeanAndPrecision(double[] rootMean, double[][] treePrecision, double rootPrecision) {

        determineSumOfSquares(rootMean, Ay, treePrecision, rootPrecision); // Fills in Ay

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
}
