package dr.evomodel.treedatalikelihood.continuous;

import dr.evomodel.continuous.MultivariateDiffusionModel;
import dr.evomodel.continuous.MultivariateElasticModel;
import dr.inference.model.MatrixParameter;
import dr.inference.model.MatrixParameterInterface;
import dr.inference.model.Parameter;
import dr.evomodel.continuous.ou.OUProcessModel;
import dr.evomodel.treedatalikelihood.continuous.gaussian.GaussianBranchTransitionKernel;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;

/**
 * Tree-side Gaussian branch-transition provider backed by the shared time-series
 * OU branch kernel.
 *
 * <p>This adapter snapshots the current tree OU parameters into a reusable
 * {@link OUProcessModel}. The branch actualization and covariance are then obtained
 * from the same branch-length machinery used by the time-series code. The branch
 * displacement is intentionally not modeled here because tree OU delegates often use
 * branch-specific optima rather than a single stationary mean.</p>
 */
public final class TimeSeriesOUGaussianBranchTransitionProvider
        implements GaussianBranchTransitionProvider {

    private final MultivariateElasticModel elasticModel;
    private final MultivariateDiffusionModel diffusionModel;
    private final int dimension;

    private final MatrixParameterInterface driftMatrix;
    private final MatrixParameter diffusionCovariance;
    private final MatrixParameter initialCovariance;
    private final Parameter zeroMean;
    private final OUProcessModel processModel;
    private final GaussianBranchTransitionKernel kernel;

    private final DenseMatrix64F precisionMatrix;
    private final DenseMatrix64F covarianceMatrix;

    public TimeSeriesOUGaussianBranchTransitionProvider(final MultivariateElasticModel elasticModel,
                                                        final MultivariateDiffusionModel diffusionModel) {
        if (elasticModel == null) {
            throw new IllegalArgumentException("elasticModel must not be null");
        }
        if (diffusionModel == null) {
            throw new IllegalArgumentException("diffusionModel must not be null");
        }
        this.elasticModel = elasticModel;
        this.diffusionModel = diffusionModel;
        this.dimension = elasticModel.getStrengthOfSelectionMatrixParameter().getRowDimension();
        if (dimension != diffusionModel.getDimension()) {
            throw new IllegalArgumentException("Elastic and diffusion dimensions must match");
        }

        this.driftMatrix = elasticModel.getStrengthOfSelectionMatrixParameter();
        this.diffusionCovariance = new MatrixParameter("treeOuBridge.diffusion", dimension, dimension);
        this.initialCovariance = new MatrixParameter("treeOuBridge.initial", dimension, dimension);
        this.zeroMean = new Parameter.Default(0.0);
        setIdentity(initialCovariance);

        this.processModel = new OUProcessModel(
                "treeOuBridge.process",
                dimension,
                driftMatrix,
                diffusionCovariance,
                zeroMean,
                initialCovariance,
                OUProcessModel.defaultCovarianceGradientMethod(driftMatrix));
        this.kernel = processModel.getRepresentation(GaussianBranchTransitionKernel.class);

        this.precisionMatrix = new DenseMatrix64F(dimension, dimension);
        this.covarianceMatrix = new DenseMatrix64F(dimension, dimension);
    }

    public GaussianBranchTransitionKernel getKernel() {
        refreshSnapshots();
        return kernel;
    }

    public OUProcessModel getProcessModel() {
        refreshSnapshots();
        return processModel;
    }

    @Override
    public int getStateDimension() {
        return dimension;
    }

    @Override
    public void getInitialMean(final double[] out) {
        refreshSnapshots();
        kernel.getInitialMean(out);
    }

    @Override
    public void getInitialCovariance(final double[][] out) {
        refreshSnapshots();
        kernel.getInitialCovariance(out);
    }

    @Override
    public void fillBranchTransitionMatrix(final double branchLength, final double[][] out) {
        refreshSnapshots();
        kernel.fillTransitionMatrix(branchLength, out);
    }

    @Override
    public void fillBranchTransitionOffset(final double branchLength, final double[] out) {
        refreshSnapshots();
        kernel.fillTransitionOffset(branchLength, out);
    }

    @Override
    public void fillBranchTransitionCovariance(final double branchLength, final double[][] out) {
        refreshSnapshots();
        kernel.fillTransitionCovariance(branchLength, out);
    }

    private void refreshSnapshots() {
        invertPrecisionIntoCovariance(diffusionModel.getPrecisionMatrix(), covarianceMatrix);
        copyIntoMatrixParameter(covarianceMatrix, diffusionCovariance);
        processModel.fireModelChanged();
    }

    private void invertPrecisionIntoCovariance(final double[][] precision, final DenseMatrix64F covarianceOut) {
        final double[] precisionData = precisionMatrix.data;
        for (int i = 0; i < dimension; ++i) {
            for (int j = 0; j < dimension; ++j) {
                precisionData[i * dimension + j] = precision[i][j];
            }
        }
        covarianceOut.set(precisionMatrix);
        CommonOps.invert(covarianceOut);
    }

    private static void copyIntoMatrixParameter(final double[][] source, final MatrixParameter destination) {
        for (int i = 0; i < source.length; ++i) {
            for (int j = 0; j < source[i].length; ++j) {
                destination.setParameterValueQuietly(i, j, source[i][j]);
            }
        }
        destination.fireParameterChangedEvent();
    }

    private static void copyIntoMatrixParameter(final DenseMatrix64F source, final MatrixParameter destination) {
        final int dimension = source.numRows;
        final double[] data = source.data;
        for (int i = 0; i < dimension; ++i) {
            for (int j = 0; j < dimension; ++j) {
                destination.setParameterValueQuietly(i, j, data[i * dimension + j]);
            }
        }
        destination.fireParameterChangedEvent();
    }

    private static void setIdentity(final MatrixParameter matrix) {
        for (int i = 0; i < matrix.getRowDimension(); ++i) {
            for (int j = 0; j < matrix.getColumnDimension(); ++j) {
                matrix.setParameterValueQuietly(i, j, i == j ? 1.0 : 0.0);
            }
        }
        matrix.fireParameterChangedEvent();
    }
}
