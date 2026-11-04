package dr.util;

import dr.xml.*;

public class CovarianceToCorrelationTransform extends Transform.MatrixVariateTransform {

    private static final String COVARIANCE_TO_CORRELATION = "covarianceToCorrelation";

    public CovarianceToCorrelationTransform(int dim) {
        super(dim * dim, dim, dim);
    }

    @Override
    public double[] inverse(double[] y, int from, int to, double sum) {
        throw new RuntimeException("inverse transform does not exist");
    }

    @Override
    public double[] gradient(double[] values, int from, int to) {
        throw new RuntimeException("not yet implemented");
    }

    @Override
    public double[] gradientInverse(double[] values, int from, int to) {
        throw new RuntimeException("inverse transform does not exist");
    }

    @Override
    public String getTransformName() {
        return COVARIANCE_TO_CORRELATION;
    }

    @Override
    protected double[] transform(double[] values) {
        double[] newValues = new double[values.length];
        int n = (int) Math.sqrt(values.length);
        for (int i = 0; i < n; i++) {
            int indii = i * n + i;
            newValues[indii] = 1;
            double xii = values[indii];

            for (int j = i + 1; j < n; j++) {
                double xjj = values[j * n + j];
                int indij = j * n + i;
                double value = values[indij] / Math.sqrt(xii * xjj);
                newValues[indij] = value;
                newValues[i * n + j] = value;

            }
        }
        return newValues;
    }

    @Override
    protected double[] inverse(double[] values) {
        throw new RuntimeException("inverse transform does not exist");
    }

    @Override
    protected double getLogJacobian(double[] values) {
        throw new RuntimeException("not yet implemented");
    }

    @Override
    protected double[] getGradientLogJacobianInverse(double[] values) {
        throw new RuntimeException("inverse transform does not exist");
    }

    @Override
    public double[][] computeJacobianMatrixInverse(double[] values) {
        throw new RuntimeException("inverse transform does not exist");
    }

    @Override
    protected boolean isInInteriorDomain(double[] values) {
        throw new RuntimeException("not yet implemented");
    }

    public static final AbstractXMLObjectParser PARSER = new AbstractXMLObjectParser() {

        private static final String DIMENSION = "dimension";


        @Override
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            int dim = xo.getIntegerAttribute(DIMENSION);
            return new CovarianceToCorrelationTransform(dim);
        }

        @Override
        public XMLSyntaxRule[] getSyntaxRules() {
            return new XMLSyntaxRule[]{
                    AttributeRule.newIntegerRule(DIMENSION)
            };
        }

        @Override
        public String getParserDescription() {
            return "Converts a covariance matrix to a correlation matrix";
        }

        @Override
        public Class getReturnType() {
            return CovarianceToCorrelationTransform.class;
        }

        @Override
        public String getParserName() {
            return COVARIANCE_TO_CORRELATION;
        }
    };
}
