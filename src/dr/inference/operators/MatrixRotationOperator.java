package dr.inference.operators;

import dr.inference.model.MatrixParameterInterface;
import dr.math.MathUtils;
import dr.xml.*;

import static dr.math.MathUtils.nextDouble;

public class MatrixRotationOperator extends AbstractAdaptableOperator {

    private double rotationWindow;
    private final MatrixParameterInterface param;

    private static final double TWO_PI = 2 * Math.PI;

    MatrixRotationOperator(MatrixParameterInterface param, double rotationWindow, double weight) {
        this.param = param;
        this.rotationWindow = rotationWindow;
        setWeight(weight);
    }

    @Override
    public String getOperatorName() {
        return MATRIX_ROTATION_OPERATOR;
    }

    @Override
    public double doOperation() {
        int nRows = param.getRowDimension();
        int row1 = MathUtils.nextInt(nRows);
        int row2 = MathUtils.nextInt(nRows);
        while (row1 == row2) {
            row2 = MathUtils.nextInt(nRows);
        }

        double theta = TWO_PI * rotationWindow * (nextDouble() - 0.5);
        double sinTheta = Math.sin(theta);
        double cosTheta = Math.cos(theta);
        for (int i = 0; i < param.getColumnDimension(); i++) {
            double originalValue1 = param.getParameterValue(row1, i);
            double originalValue2 = param.getParameterValue(row2, i);
            param.setParameterValueQuietly(row1, i, originalValue1 * cosTheta - originalValue2 * sinTheta);
            param.setParameterValueQuietly(row2, i, originalValue1 * sinTheta + originalValue2 * cosTheta);
        }

        param.fireParameterChangedEvent();

        return 0; //symmetric proposals
    }

    @Override
    protected void setAdaptableParameterValue(double value) {
        rotationWindow = value;

    }

    @Override
    protected double getAdaptableParameterValue() {
        return rotationWindow;
    }

    @Override
    public double getRawParameter() {
        return rotationWindow;
    }

    @Override
    public String getAdaptableParameterName() {
        return "Rotation fraction";
    }

    public static final String MATRIX_ROTATION_OPERATOR = "matrixRotationOperator";


    public static final AbstractXMLObjectParser PARSER = new AbstractXMLObjectParser() {

        private static final String ROTATION_WINDOW = "rotationWindow";

        @Override
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            MatrixParameterInterface matParam = (MatrixParameterInterface) xo.getChild(MatrixParameterInterface.class);
            double rotationWindow = xo.getDoubleAttribute(ROTATION_WINDOW);
            double weight = xo.getDoubleAttribute(WEIGHT);
            return new MatrixRotationOperator(matParam, rotationWindow, weight);
        }

        @Override
        public XMLSyntaxRule[] getSyntaxRules() {
            return new XMLSyntaxRule[]{
                    new ElementRule(MatrixParameterInterface.class),
                    AttributeRule.newDoubleRule(ROTATION_WINDOW),
                    AttributeRule.newDoubleRule(WEIGHT)
            };
        }

        @Override
        public String getParserDescription() {
            return null;
        }

        @Override
        public Class getReturnType() {
            return MatrixRotationOperator.class;
        }

        @Override
        public String getParserName() {
            return MATRIX_ROTATION_OPERATOR;
        }
    };
}
