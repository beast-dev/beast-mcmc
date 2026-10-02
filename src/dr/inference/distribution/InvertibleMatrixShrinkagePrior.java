package dr.inference.distribution;

import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.AbstractModelLikelihood;
import dr.inference.model.InvertibleMatrixParametrization;
import dr.inference.model.Likelihood;
import dr.inference.model.MatrixParameter;
import dr.inference.model.Model;
import dr.inference.model.NativeMatrixGradientProvider;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;

/**
 * Prior on a square matrix R that both repels R from singularity and bounds its
 * overall scale:
 * <p>
 * p(R) &prop; |det R|^lambda exp(-lambda * c / 2 ||R||_F^2), lambda &gt; 0
 * <p>
 * log p(R) = C + lambda * log|det R| - lambda * c / 2 ||R||_F^2
 * <p>
 * with default c = 1 for backward compatibility and optional c = dim(R) for
 * dimension-scaled shrinkage.
 * <p>
 * The dense gradient is d/dR log p(R) = lambda * (R^-T - c R). When R exposes native
 * HMC coordinates through {@link NativeMatrixGradientProvider}, the dense
 * matrix gradient is pulled back to those native coordinates.
 */
public class InvertibleMatrixShrinkagePrior extends AbstractModelLikelihood
        implements GradientWrtParameterProvider {

    private final MatrixParameter matrix;
    private final InvertibleMatrixParametrization invertibleMatrix;
    private final NativeMatrixGradientProvider nativeGradientProvider;
    private final Parameter lambda;
    private final Parameter targetParameter;
    private final int dim;
    private final double frobeniusPenaltyScale;

    private final double[] matrixRowMajor;
    private final double[] inverseRowMajor;
    private final double[] gradientRowMajor;
    private final double[] denseGradientParameterOrder;
    private final DenseMatrix64F matrixView;
    private final DenseMatrix64F inverseView;

    private boolean likelihoodKnown = false;
    private double logLikelihood;
    private boolean storedLikelihoodKnown;
    private double storedLogLikelihood;

    public InvertibleMatrixShrinkagePrior(final String name,
                                          final MatrixParameter matrix,
                                          final Parameter lambda) {
        this(name, matrix, lambda, false);
    }

    public InvertibleMatrixShrinkagePrior(final String name,
                                          final MatrixParameter matrix,
                                          final Parameter lambda,
                                          final boolean scaleFrobeniusByDimension) {
        super(name);
        if (matrix == null) {
            throw new IllegalArgumentException("matrix is required");
        }
        if (matrix.getRowDimension() != matrix.getColumnDimension()) {
            throw new IllegalArgumentException("matrix must be square");
        }
        if (matrix.getDimension() != matrix.getRowDimension() * matrix.getColumnDimension()) {
            throw new IllegalArgumentException("matrix must expose all dense entries as parameters");
        }
        if (lambda == null) {
            throw new IllegalArgumentException("lambda is required");
        }
        if (lambda.getDimension() != 1) {
            throw new IllegalArgumentException("lambda must be a scalar parameter");
        }
        this.matrix = matrix;
        this.dim = matrix.getRowDimension();
        this.frobeniusPenaltyScale = scaleFrobeniusByDimension ? dim : 1.0;
        this.invertibleMatrix = matrix instanceof InvertibleMatrixParametrization
                ? (InvertibleMatrixParametrization) matrix
                : null;
        this.nativeGradientProvider = matrix instanceof NativeMatrixGradientProvider
                ? (NativeMatrixGradientProvider) matrix
                : null;
        this.lambda = lambda;
        this.targetParameter = nativeGradientProvider == null
                ? matrix
                : nativeGradientProvider.getNativeParameter();
        if (targetParameter == null) {
            throw new IllegalArgumentException("matrix native parameter is required");
        }
        this.matrixRowMajor = new double[dim * dim];
        this.inverseRowMajor = new double[dim * dim];
        this.gradientRowMajor = new double[dim * dim];
        this.denseGradientParameterOrder = new double[dim * dim];
        this.matrixView = DenseMatrix64F.wrap(dim, dim, matrixRowMajor);
        this.inverseView = DenseMatrix64F.wrap(dim, dim, inverseRowMajor);

        addVariable(matrix);
        addVariable(lambda);
    }

    @Override
    public double getLogLikelihood() {
        if (!likelihoodKnown) {
            logLikelihood = computeLogLikelihood();
            likelihoodKnown = true;
        }
        return logLikelihood;
    }

    private double computeLogLikelihood() {
        final double lam = getLambda();
        fillMatrixRowMajor();
        double frobeniusSquared = 0.0;
        for (final double value : matrixRowMajor) {
            frobeniusSquared += value * value;
        }
        return lam * getLogAbsDeterminant() - 0.5 * lam * frobeniusPenaltyScale * frobeniusSquared;
    }

    @Override
    public double[] getGradientLogDensity() {
        final double lam = getLambda();
        fillMatrixRowMajor();
        fillInverse();

        for (int row = 0; row < dim; ++row) {
            for (int col = 0; col < dim; ++col) {
                final double inverseTranspose = inverseRowMajor[col * dim + row];
                gradientRowMajor[row * dim + col] = lam *
                        (inverseTranspose - frobeniusPenaltyScale * matrixRowMajor[row * dim + col]);
            }
        }

        if (nativeGradientProvider != null) {
            return nativeGradientProvider.pullBackGradientForParameter(targetParameter, gradientRowMajor, dim);
        }

        fillDenseGradientParameterOrder();
        return denseGradientParameterOrder.clone();
    }

    @Override
    public Likelihood getLikelihood() {
        return this;
    }

    @Override
    public Parameter getParameter() {
        return targetParameter;
    }

    @Override
    public int getDimension() {
        return targetParameter.getDimension();
    }

    @Override
    public Model getModel() {
        return this;
    }

    @Override
    public void makeDirty() {
        likelihoodKnown = false;
    }

    @Override
    protected void handleModelChangedEvent(final Model model, final Object object, final int index) {
        likelihoodKnown = false;
    }

    @Override
    protected void handleVariableChangedEvent(final Variable variable, final int index, final Parameter.ChangeType type) {
        likelihoodKnown = false;
    }

    @Override
    protected void storeState() {
        storedLikelihoodKnown = likelihoodKnown;
        storedLogLikelihood = logLikelihood;
    }

    @Override
    protected void restoreState() {
        likelihoodKnown = storedLikelihoodKnown;
        logLikelihood = storedLogLikelihood;
    }

    @Override
    protected void acceptState() {
        // nothing to do
    }

    private double getLambda() {
        final double lam = lambda.getParameterValue(0);
        if (!Double.isFinite(lam) || lam < 0.0) {
            throw new IllegalArgumentException("lambda must be finite and non-negative");
        }
        return lam;
    }

    private void fillMatrixRowMajor() {
        int index = 0;
        for (int row = 0; row < dim; ++row) {
            for (int col = 0; col < dim; ++col) {
                matrixRowMajor[index++] = matrix.getParameterValue(row, col);
            }
        }
    }

    private void fillInverse() {
        if (invertibleMatrix != null) {
            invertibleMatrix.fillInverse(inverseRowMajor);
        } else if (!CommonOps.invert(matrixView, inverseView)) {
            throw new IllegalArgumentException("matrix is singular and cannot be inverted");
        }
    }

    private double getLogAbsDeterminant() {
        if (invertibleMatrix != null) {
            return invertibleMatrix.getLogAbsDeterminant();
        }

        final double determinant = CommonOps.det(matrixView);
        if (!Double.isFinite(determinant) || determinant == 0.0) {
            throw new IllegalArgumentException("matrix determinant is zero or non-finite");
        }
        return Math.log(Math.abs(determinant));
    }

    private void fillDenseGradientParameterOrder() {
        int index = 0;
        for (int col = 0; col < dim; ++col) {
            for (int row = 0; row < dim; ++row) {
                denseGradientParameterOrder[index++] = gradientRowMajor[row * dim + col];
            }
        }
    }
}
