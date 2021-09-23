package dr.inference.model;

import dr.xml.*;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;

/**
 * @author Gabriel Hassler
 * @author Marc A. Suchard
 */

public class DeterminantStatistic extends Statistic.Abstract implements VariableListener {

    private final MatrixParameterInterface matrix;
    private final int matrixDim;
    private boolean detKnown = false;
    private double det;

    public DeterminantStatistic(String name, MatrixParameterInterface matrix) {
        super(name);

        this.matrix = matrix;
        this.matrixDim = matrix.getRowDimension();
        matrix.addParameterListener(this);

    }

    @Override
    public int getDimension() {
        return 1;
    }

    @Override
    public double getStatisticValue(int dim) {
        if (!detKnown) {
            double[] values = matrix.getParameterValues();
            DenseMatrix64F M = DenseMatrix64F.wrap(matrixDim, matrixDim, values);
            det = CommonOps.det(M);
            detKnown = true;
        }

        return det;
    }

    @Override
    public void variableChangedEvent(Variable variable, int index, Variable.ChangeType type) {
        detKnown = false;
    }

    public static final AbstractXMLObjectParser PARSER = new AbstractXMLObjectParser() {

        private static final String DETERMINANT_STATISTIC = "determinant";


        @Override
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            MatrixParameterInterface matrix = (MatrixParameterInterface) xo.getChild(MatrixParameterInterface.class);
            if (matrix.getColumnDimension() != matrix.getRowDimension()) {
                throw new XMLParseException("can only calculate determinant for square matrices");
            }
            final String name;
            if (xo.hasId()) {
                name = xo.getId();
            } else if (matrix.getId() != null) {
                name = "determinant." + matrix.getId();
            } else {
                name = "determinant";
            }
            return new DeterminantStatistic(name, matrix);
        }

        @Override
        public XMLSyntaxRule[] getSyntaxRules() {
            return new XMLSyntaxRule[]{
                    new ElementRule(MatrixParameterInterface.class)
            };
        }

        @Override
        public String getParserDescription() {
            return "Statistic that computes the determinant of a matrix";
        }

        @Override
        public Class getReturnType() {
            return DeterminantStatistic.class;
        }

        @Override
        public String getParserName() {
            return DETERMINANT_STATISTIC;
        }
    };


}
