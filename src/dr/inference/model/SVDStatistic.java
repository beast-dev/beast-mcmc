package dr.inference.model;

import dr.math.matrixAlgebra.Vector;
import dr.util.HeapSort;
import dr.xml.*;
import org.ejml.data.DenseMatrix64F;
import org.ejml.factory.DecompositionFactory;
import org.ejml.interfaces.decomposition.SingularValueDecomposition;

/**
 * @author Gabriel Hassler
 */

public class SVDStatistic extends Statistic.Abstract implements VariableListener, Reportable {

    private final MatrixParameterInterface parameter;
    private final double[] singularVals;
    private final double[] Vbuffer;
    private boolean svdKnown = false;

    SVDStatistic(MatrixParameterInterface parameter) {

        this.parameter = parameter;
        this.singularVals = new double[parameter.getColumnDimension()];
        this.Vbuffer = new double[parameter.getRowDimension() * parameter.getColumnDimension()];

        parameter.addParameterListener(this);


    }

    private static String SINGULAR_VALUE = "sv";
    private static String V = "V";

    @Override
    public String getDimensionName(int dim) {
        int k = parameter.getColumnDimension();
        int p = parameter.getRowDimension();

        if (dim < k) {
            return getStatisticName() + "." + SINGULAR_VALUE + (dim + 1);
        } else {
            int row = (dim - k) / p;
            int col = dim - k - row * p;
            return getStatisticName() + "." + V + (row + 1) + (col + 1);
        }
    }


    @Override
    public int getDimension() {
        int k = parameter.getColumnDimension();
        int p = parameter.getRowDimension();

        return k + k * p;
    }

    @Override
    public double getStatisticValue(int dim) {

        if (!svdKnown) {
            recomputeStatistic();
            enforceConstraints();
            svdKnown = true;
        }

        int k = parameter.getColumnDimension();


        if (dim < k) {
            return singularVals[dim];
        } else {
            return Vbuffer[dim - k];
        }

    }

    private void recomputeStatistic() {

        double[] values = parameter.getParameterValues();

        int k = parameter.getColumnDimension();
        int p = parameter.getRowDimension();

        DenseMatrix64F matrix = DenseMatrix64F.wrap(k, p, values);
        DenseMatrix64F buffer = DenseMatrix64F.wrap(k, p, Vbuffer);

        SingularValueDecomposition svd = DecompositionFactory.svd(matrix.numRows, matrix.numCols,
                false, true, true);

        svd.decompose(matrix);

        System.arraycopy(svd.getSingularValues(), 0, singularVals, 0, k);
        svd.getV(buffer, true);

    }

    private void enforceConstraints() {
        int k = parameter.getColumnDimension();
        int p = parameter.getRowDimension();

        for (int i = 1; i < k; i++) {
            if (singularVals[i] > singularVals[i - 1]) {
                reorderSVD();
                break;
            }
        }

        for (int i = 0; i < k; i++) {
            int offset = i * p;
            if (Vbuffer[offset] < 0) {
                for (int j = offset; j < offset + p; j++) {
                    Vbuffer[j] = -Vbuffer[j];
                }
            }
        }

    }

    private void reorderSVD() {
        int k = parameter.getColumnDimension();
        int p = parameter.getRowDimension();
        int[] indices = new int[k];

        HeapSort.sort(singularVals, indices);

        for (int i = 0; i < k; i++) { // want descending order
            indices[i] = k - indices[i] - 1;
        }

        double[] kBuffer = new double[k];
        double[] kpBuffer = new double[k * p];
        for (int i = 0; i < k; i++) {
            kBuffer[i] = singularVals[indices[i]];
            for (int j = 0; j < p; j++) {
                kpBuffer[i * p + j] = Vbuffer[indices[i] * p + j];
            }
        }

        System.arraycopy(kBuffer, 0, singularVals, 0, k);
        System.arraycopy(kpBuffer, 0, Vbuffer, 0, k * p);
    }

    @Override
    public void variableChangedEvent(Variable variable, int index, Variable.ChangeType type) {
        svdKnown = false;
    }

    @Override
    public String getReport() {

        StringBuilder sb = new StringBuilder("svdStatistic Report\n\n");
        sb.append("dimension names: ");

        int n = getDimension();

        double[] values = new double[n];
        for (int i = 0; i < n; i++) {
            sb.append(getDimensionName(i));
            if (i != n - 1) {
                sb.append(" ");
            }
            values[i] = getStatisticValue(i);
        }
        sb.append("\n\n");
        sb.append("values: ");
        sb.append(new Vector(values));
        sb.append("\n\n");
        return sb.toString();

    }

    public static final String SVD_STATISTIC = "svdStatistic";

    public static AbstractXMLObjectParser PARSER = new AbstractXMLObjectParser() {
        @Override
        public Object parseXMLObject(XMLObject xo) {
            MatrixParameterInterface parameter = (MatrixParameterInterface) xo.getChild(MatrixParameterInterface.class);
            return new SVDStatistic(parameter);
        }

        @Override
        public XMLSyntaxRule[] getSyntaxRules() {
            return new XMLSyntaxRule[]{
                    new ElementRule(MatrixParameterInterface.class)
            };
        }

        @Override
        public String getParserDescription() {
            return "Statistic that computes and returns the SVD of a given parameter.";
        }

        @Override
        public Class getReturnType() {
            return SVDStatistic.class;
        }

        @Override
        public String getParserName() {
            return SVD_STATISTIC;
        }
    };


}
