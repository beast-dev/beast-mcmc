package dr.util;

import dr.xml.AbstractXMLObjectParser;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

public class BivariateDifferenceTransform extends Transform.MultivariateTransform {
    public static String NAME = "BivariateDifferenceTransform";

    BivariateDifferenceTransform() {
        super(2, 1);
    }

    @Override
    public double[] inverse(double[] values, int from, int to, double sum) {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public double[] gradient(double[] values, int from, int to) {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public double[] gradientInverse(double[] values, int from, int to) {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public String getTransformName() {
        return NAME;
    }

    @Override
    protected double[] transform(double[] values) {
        double[] out = new double[1];
        out[0] = values[0] - values[1];
        return out;
    }

    @Override
    protected double[] inverse(double[] values) {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    protected double getLogJacobian(double[] values) {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    protected double[] getGradientLogJacobianInverse(double[] values) {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public double[][] computeJacobianMatrixInverse(double[] values) {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    protected boolean isInInteriorDomain(double[] values) {
        return true;
    }

    public static final AbstractXMLObjectParser PARSER = new AbstractXMLObjectParser() {

        @Override
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            final String name = xo.hasId() ? xo.getId() : null;

            BivariateDifferenceTransform transform = new BivariateDifferenceTransform();

            return transform;
        }

        @Override
        public XMLSyntaxRule[] getSyntaxRules() {
            return new XMLSyntaxRule[0];
        }

        @Override
        public String getParserDescription() {
            return null;
        }

        @Override
        public Class getReturnType() {
            return Transform.MultivariateTransform.class;
        }

        @Override
        public String getParserName() {
            return NAME;
        }
    };
}
