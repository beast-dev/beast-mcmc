package dr.inference.operators.factorAnalysis;

import dr.inference.model.MatrixParameterInterface;
import dr.inference.operators.SimpleMCMCOperator;
import dr.math.MathUtils;
import dr.xml.*;

public class ColumnSwapOperator extends SimpleMCMCOperator {

    private final MatrixParameterInterface matrix;
    private static final String COLUMN_SWAP = "columnSwapOperator";

    public ColumnSwapOperator(MatrixParameterInterface matrix, double weight) {
        setWeight(weight);
        this.matrix = matrix;
    }


    @Override
    public String getOperatorName() {
        return COLUMN_SWAP;
    }

    @Override
    public double doOperation() {
        int n = matrix.getRowDimension();
        int p = matrix.getColumnDimension();

        int i1 = MathUtils.nextInt(p);
        int i2 = MathUtils.nextInt(p);
        while (i2 == i1) {
            i2 = MathUtils.nextInt(p);
        }

        for (int i = 0; i < n; i++) {
            double v1 = matrix.getParameterValue(i, i1);
            double v2 = matrix.getParameterValue(i, i2);
            matrix.setParameterValueQuietly(i, i1, v2);
            matrix.setParameterValueQuietly(i, i2, v1);
        }
        matrix.fireParameterChangedEvent();
        return 0;
    }

    private static final String WEIGHT = "weight";

    public static AbstractXMLObjectParser PARSER = new AbstractXMLObjectParser() {
        @Override
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            double weight = xo.getDoubleAttribute(WEIGHT);
            MatrixParameterInterface matrix = (MatrixParameterInterface) xo.getChild(MatrixParameterInterface.class);

            return new ColumnSwapOperator(matrix, weight);
        }

        @Override
        public XMLSyntaxRule[] getSyntaxRules() {
            return new XMLSyntaxRule[]{
                    new ElementRule(MatrixParameterInterface.class),
                    AttributeRule.newDoubleRule(WEIGHT)
            };
        }

        @Override
        public String getParserDescription() {
            return "swaps random columns of a matrix";
        }

        @Override
        public Class getReturnType() {
            return ColumnSwapOperator.class;
        }

        @Override
        public String getParserName() {
            return "columnSwapOperator";
        }
    };
}
