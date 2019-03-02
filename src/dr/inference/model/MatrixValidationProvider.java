package dr.inference.model;

import dr.xml.*;

public class MatrixValidationProvider implements CrossValidationProvider {

    final MatrixParameter trueParameter;
    final MatrixParameter inferredParameter;
    final int dimParameter;
    final int[] relevantDimensions;
    final String[] colNames;

    MatrixValidationProvider(MatrixParameter trueParameter, MatrixParameter inferredParameter, String id) {
        this.trueParameter = trueParameter;
        this.inferredParameter = inferredParameter;

        this.dimParameter = trueParameter.getDimension();

        this.relevantDimensions = new int[dimParameter];
        for (int i = 0; i < dimParameter; i++) {
            relevantDimensions[i] = i;
        }

        this.colNames = new String[dimParameter];

        for (int i = 0; i < dimParameter; i++) {

            int row = i / trueParameter.getRowDimension();
            int col = i - row * trueParameter.getRowDimension();

            colNames[i] = id + (row + 1) + (col + 1);

        }


    }

    @Override
    public Parameter getTrueParameter() {
        return trueParameter;
    }

    @Override
    public Parameter getInferredParameter() {
        return inferredParameter;
    }

    @Override
    public int[] getRelevantDimensions() {
        return relevantDimensions;
    }

    @Override
    public String getName(int dim) {
        return colNames[dim];
    }

    public static dr.xml.XMLObjectParser PARSER = new dr.xml.AbstractXMLObjectParser() {

        final static String PARSER_NAME = "matrixValidation";
        final static String TRUE_PARAMETER = "trueParameter";
        final static String INFERRED_PARAMETER = "inferredParameter";

        @Override
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            MatrixParameter trueParameter = (MatrixParameter) xo.getElementFirstChild(TRUE_PARAMETER);
            MatrixParameter inferredParameter = (MatrixParameter) xo.getElementFirstChild(INFERRED_PARAMETER);


            String id = PARSER_NAME;
            if (xo.hasId()) {
                id = xo.getId();
            }

            if (trueParameter.getRowDimension() != inferredParameter.getRowDimension() ||
                    trueParameter.getColumnDimension() != inferredParameter.getColumnDimension()) {

                throw new XMLParseException("The matrix parameters contained in " + TRUE_PARAMETER + " and " + INFERRED_PARAMETER +
                        " must have the same dimensions.");
            }

            MatrixValidationProvider provider = new MatrixValidationProvider(trueParameter, inferredParameter, id);

            CrossValidator crossValidator = new CrossValidator(provider);

            return crossValidator;

        }

        @Override
        public XMLSyntaxRule[] getSyntaxRules() {
            return new XMLSyntaxRule[]{

            };
        }

        @Override
        public String getParserDescription() {
            return null;
        }

        @Override
        public Class getReturnType() {
            return CrossValidator.class;
        }

        @Override
        public String getParserName() {
            return PARSER_NAME;
        }


    };
}
