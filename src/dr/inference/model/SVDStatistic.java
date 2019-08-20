package dr.inference.model;

import dr.math.matrixAlgebra.Vector;
import dr.xml.*;
import org.ejml.data.DenseMatrix64F;
import org.ejml.data.Matrix;
import org.ejml.factory.DecompositionFactory;
import org.ejml.interfaces.decomposition.SingularValueDecomposition;

public class SVDStatistic extends Statistic.Abstract implements VariableListener, Reportable {

    private final MatrixParameterInterface parameter;
    private final double[] singularVals;
    private final double[] Vbuffer;
    private boolean svdKnown = false;

    SVDStatistic(MatrixParameterInterface parameter) {

        this.parameter = parameter;
        this.singularVals = new double[parameter.getRowDimension()];
        this.Vbuffer = new double[parameter.getRowDimension() * parameter.getColumnDimension()];

        parameter.addParameterListener(this);


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
        Matrix V = svd.getV(buffer, true);
    }

    @Override
    public void variableChangedEvent(Variable variable, int index, Variable.ChangeType type) {
        svdKnown = false;
    }

    @Override
    public String getReport() {
        StringBuilder sb = new StringBuilder();

        int n = getDimension();
        double[] values = new double[n];
        for (int i = 0; i < n; i++) {
            values[i] = getStatisticValue(i);
        }
        sb.append("values: ");
        sb.append(new Vector(values));
        sb.append("\n\n");
        return sb.toString();

    }

    public static final String SVD_STATISTIC = "svdStatistic";

    public static AbstractXMLObjectParser PARSER = new AbstractXMLObjectParser() {
        @Override
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
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
