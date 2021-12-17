package dr.inference.model;

import dr.xml.*;

public class MatrixDiagonalLogger extends Statistic.Abstract {
    private final MatrixParameterInterface matrix;

    public MatrixDiagonalLogger(MatrixParameterInterface matrix) {
        this.matrix = matrix;
    }


    @Override
    public int getDimension() {
        return matrix.getColumnDimension();
    }

    @Override
    public double getStatisticValue(int dim) {
        return matrix.getParameterValue(dim, dim);
    }

    @Override
    public String getDimensionName(int dim) {
        return getStatisticName() + "." + matrix.getDimensionName(dim);
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        private static final String MATRIX_DIAGONAL = "matrixDiagonals";


        @Override
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            MatrixParameterInterface matrix = (MatrixParameterInterface) xo.getChild(MatrixParameterInterface.class);
            if (matrix.getColumnDimension() != matrix.getRowDimension()) {
                throw new XMLParseException("Only square matrices can be converted to correlation matrices");
            }

            return new MatrixDiagonalLogger(matrix);
        }

        @Override
        public XMLSyntaxRule[] getSyntaxRules() {
            return new XMLSyntaxRule[]{
                    new ElementRule(MatrixParameterInterface.class)
            };
        }

        @Override
        public String getParserDescription() {
            return "This element returns a statistic that is the diagonals of the associated matrix.";
        }

        @Override
        public Class getReturnType() {
            return MatrixDiagonalLogger.class;
        }

        @Override
        public String getParserName() {
            return MATRIX_DIAGONAL;
        }
    };
}
