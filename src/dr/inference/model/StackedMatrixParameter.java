package dr.inference.model;

import dr.xml.*;

public class StackedMatrixParameter extends MatrixParameter {
    public StackedMatrixParameter(String name, Parameter[] parameters) {
        super(name, parameters);
    }

    @Override
    public double getParameterValue(int dim) {
        int colDim = dim / rowDimension;
        int rowDim = dim - rowDimension * colDim;

        return getParameterValue(rowDim, colDim);
    }


    public static final String STACKED_MATRIX_PARAMETER = "stackedMatrixParameter";

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {
        @Override
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            final String name = xo.hasId() ? xo.getId() : null;
            Parameter[] params = new Parameter[xo.getChildCount()];
            int dim = -1;
            for (int i = 0; i < xo.getChildCount(); i++) {
                Parameter param = (Parameter) xo.getChild(i);
                if (i == 0) {
                    dim = param.getDimension();
                } else {
                    if (param.getDimension() != dim) {
                        throw new XMLParseException("All parameters must have the same dimension to construct a" +
                                " rectangular matrix");
                    }

                }

                params[i] = param;
            }

            return new StackedMatrixParameter(name, params);
        }

        @Override
        public XMLSyntaxRule[] getSyntaxRules() {
            return new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class, 1, Integer.MAX_VALUE)
            };
        }

        @Override
        public String getParserDescription() {
            return "Matrix parameter with unusual indexing";
        }

        @Override
        public Class getReturnType() {
            return StackedMatrixParameter.class;
        }

        @Override
        public String getParserName() {
            return STACKED_MATRIX_PARAMETER;
        }
    };
}
