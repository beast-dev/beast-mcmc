package dr.inference.distribution;

import dr.inference.distribution.NormalStatisticsHelpers.IndependentNormalStatisticsProvider;
import dr.inference.model.*;
import dr.inferencexml.distribution.NormalOrthogonalSubspaceDistributionParser;
import dr.math.distributions.NormalDistribution;
import dr.math.matrixAlgebra.Matrix;
import dr.xml.Reportable;
import org.ejml.data.DenseMatrix64F;
import org.ejml.factory.DecompositionFactory;
import org.ejml.interfaces.decomposition.SingularValueDecomposition;
import org.ejml.ops.CommonOps;

import java.util.Arrays;


public class NormalOrthogonalSubspaceDistribution extends AbstractModelLikelihood
        implements NormalStatisticsHelpers.MatrixNormalStatisticsHelper, IndependentNormalStatisticsProvider,
        Reportable {

    private final IndependentNormalStatisticsProvider priorDistribution;
    private final MatrixParameterInterface matrix;
    private Boolean rotationKnown = false;
    private final int nRows;
    private final int nCols;
    private final SingularValueDecomposition svd;
    private final DenseMatrix64F uBuffer;
    private final DenseMatrix64F vBuffer;
    private final DenseMatrix64F wBuffer;
    private final DenseMatrix64F priorVarBuffer;
    private final DenseMatrix64F rotatedPrecision;
    private final double[] meanBuffer;

    private final double[][] precisionArrayBuffer;

    private final DenseMatrix64F kBuffer;

    private final NormalStatisticsHelpers.IndependentToMatrixAdaptor adaptor;

    public NormalOrthogonalSubspaceDistribution(String name, IndependentNormalStatisticsProvider priorDistribution,
                                                MatrixParameterInterface matrix) {
        super(name);
        this.priorDistribution = priorDistribution;
        this.matrix = matrix;
        this.nRows = matrix.getColumnDimension();
        this.nCols = matrix.getRowDimension();
        this.svd = DecompositionFactory.svd(nRows, nCols, true, true, true);
        this.uBuffer = new DenseMatrix64F(nRows, nRows);
        this.vBuffer = new DenseMatrix64F(nRows, nCols);
        this.wBuffer = new DenseMatrix64F(nRows, nRows);
        this.priorVarBuffer = new DenseMatrix64F(nRows, nRows);
        this.rotatedPrecision = new DenseMatrix64F(nRows, nRows);
        this.kBuffer = new DenseMatrix64F(nRows, nRows);
        this.precisionArrayBuffer = new double[nRows][nRows];
        this.meanBuffer = new double[nRows];

        for (int i = 0; i < matrix.getDimension(); i++) {
            if (priorDistribution.getNormalMean(i) != 0) {
                throw new RuntimeException("Not implemented for prior with non-zero mean.");
            }
        }


        addVariable(matrix);

        this.adaptor = new NormalStatisticsHelpers.IndependentToMatrixAdaptor(priorDistribution, nRows, nCols);
    }


    @Override
    public double getScalarPrecision() {
        return priorDistribution.matrixNormalHelper(nRows, nCols).getScalarPrecision();
    }

    @Override
    public double[] precisionMeanProduct(int col) {
        Arrays.fill(meanBuffer, 0);
        return meanBuffer;
    }

    @Override
    public double[][] getColumnPrecision(int col) {
        updateRotatedPrecision(col);

        return precisionArrayBuffer;
    }

    private void updateRotatedPrecision(int col) {  //TODO: cache based on prior variance
        if (!rotationKnown) {
            updateRotation();
        }
        for (int i = 0; i < nRows; i++) {
            priorVarBuffer.set(i, i, 1.0 / adaptor.getColumnPrecisionDiagonal(i, col));
        }

        CommonOps.mult(uBuffer, priorVarBuffer, rotatedPrecision); //rotatedPrecision just a buffer here
        CommonOps.multTransB(rotatedPrecision, uBuffer, kBuffer);
        CommonOps.invert(kBuffer, rotatedPrecision); //TODO: use cholesky to invert?

        for (int i = 0; i < nRows; i++) {
            for (int j = 0; j < nRows; j++) {
                precisionArrayBuffer[i][j] = rotatedPrecision.get(i, j);
            }
        }
    }

    private void updateRotation() {
        double[] matBuffer = matrix.getParameterValues();
        DenseMatrix64F mat = DenseMatrix64F.wrap(nRows, nCols, matBuffer);
        svd.decompose(mat);

        svd.getU(uBuffer, false);
        svd.getV(vBuffer, true);
        svd.getW(wBuffer);

        rotationKnown = true;
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
        // No models
    }

    @Override
    protected void storeState() {
        // do nothing
    }

    @Override
    protected void restoreState() {
        // do nothing
    }

    @Override
    protected void acceptState() {
        // do nothing
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        if (variable == matrix) {
            makeDirty();
        }
    }

    @Override
    public Model getModel() {
        return this;
    }

    @Override
    public double getLogLikelihood() {

        if (!rotationKnown) {
            updateRotation();
        }

        double logLikelihood = 0;
        for (int row = 0; row < nRows; row++) {
            for (int col = 0; col < nCols; col++) {
                double prec = adaptor.getColumnPrecisionDiagonal(row, col);
                double sd = 1.0 / Math.sqrt(prec);
                logLikelihood += NormalDistribution.logPdf(
                        matrix.getParameterValue(col, row), adaptor.getNormalMean(col, row), sd); //assumes matrix transposed
            }
        }

        return logLikelihood;
    }

    @Override
    public void makeDirty() {
        rotationKnown = false;
    }

    @Override
    public double getNormalMean(int dim) {
        return priorDistribution.getNormalMean(dim);
    }

    @Override
    public double getNormalPrecision(int dim) {
        return priorDistribution.getNormalPrecision(dim);
    }

    @Override
    public NormalStatisticsHelpers.MatrixNormalStatisticsHelper matrixNormalHelper(int nRows, int nCols) {
        if (nRows != nRows || nCols != nCols) throw new RuntimeException("Incompatible dimensions.");
        return this;
    }

    @Override
    public String getReport() {
        StringBuilder sb = new StringBuilder(NormalOrthogonalSubspaceDistributionParser.PARSER + " report:\n");
        sb.append("\tlikelihood: " + getLogLikelihood() + "\n");

        sb.append("\trotated precisions:\n");
        for (int i = 0; i < nCols; i++) {
            sb.append("\t\tcolumn " + (i + 1) + ":\n");
            Matrix pMat = new Matrix(this.getColumnPrecision(i));
            sb.append(pMat.toString(3));
            sb.append("\n");
        }

        sb.append("\n\n");
        return sb.toString();
    }
}
