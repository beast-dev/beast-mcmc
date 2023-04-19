package dr.util;

import dr.inference.model.Parameter;
import dr.inference.model.TransformedMultivariateParameter;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

public class TransformedVectorSumTransform extends Transform.MultivariateTransform {

    private static final String NAME = "transformedVectorSumTransform";
    private final Transform incrementTransform;
    private final int dim;

    public TransformedVectorSumTransform(int dim, Transform incrementTransform) {
        super(dim);
        this.dim = dim;
        this.incrementTransform = incrementTransform;
    }

    @Override
    public double[] inverse(double[] values, int from, int to, double sum) {
        throw new RuntimeException("Not yet implemented.");
    }

    @Override
    public double[] gradient(double[] values, int from, int to) {
        throw new RuntimeException("Not yet implemented.");
    }

    @Override
    public double[] gradientInverse(double[] values, int from, int to) {
        throw new RuntimeException("Not yet implemented.");
    }

    @Override
    public String getTransformName() {
        return NAME;
    }

    @Override
    protected double[] transform(double[] values) {
        double[] fx = new double[values.length];
        fx[0] = values[0];
        for (int i = 1; i < values.length; i++) {
            fx[i] = fx[i-1] + values[i];
        }
        return incrementTransform.inverse(fx, 0, values.length);
    }

    @Override
    protected double[] inverse(double[] values) {
        double[] increments = incrementTransform.transform(values,0, values.length);
        for (int i = 1; i < values.length; i++) {
            increments[i] = increments[i] - increments[i - 1];
        }
        return increments;
    }

    @Override
    protected double getLogJacobian(double[] values) {
        throw new RuntimeException("Not yet implemented.");
    }

    @Override
    protected double[] getGradientLogJacobianInverse(double[] values) {
        throw new RuntimeException("Not yet implemented.");
    }

    @Override
    public double[][] computeJacobianMatrixInverse(double[] values) {
        throw new RuntimeException("Not yet implemented.");
    }

    @Override
    protected boolean isInInteriorDomain(double[] values) {
        throw new RuntimeException("Not yet implemented.");
    }

    public static final AbstractXMLObjectParser PARSER = new AbstractXMLObjectParser() {
        public final static String INCREMENT_TRANSFORM = "incrementTransformType";

        @Override
        public TransformedMultivariateParameter parseXMLObject(XMLObject xo) throws XMLParseException {
            final String name = xo.hasId() ? xo.getId() : null;

            Parameter param = (Parameter) xo.getChild(Parameter.class);

            int dim = param.getDimension();

            double upper = Double.POSITIVE_INFINITY;
            double lower = Double.NEGATIVE_INFINITY;
            if( xo.hasAttribute("upper") && xo.hasAttribute("lower")) {
                upper = xo.getDoubleAttribute("upper");
                lower = xo.getDoubleAttribute("lower");
            }

            Transform incrementTransform = null;
            String ttype = (String) xo.getAttribute(INCREMENT_TRANSFORM);
            if (ttype.equalsIgnoreCase("log")) {
                incrementTransform = Transform.LOG;
            } else if (ttype.equalsIgnoreCase("logit")) {
                incrementTransform = new Transform.ScaledLogitTransform(upper, lower);
            } else if (ttype.equalsIgnoreCase("none")) {
                incrementTransform = new Transform.NoTransform();
            } else {
                throw new RuntimeException("Invalid transform type");
            }

            TransformedVectorSumTransform transform = new TransformedVectorSumTransform(dim, incrementTransform);

            return new TransformedMultivariateParameter(param, transform);
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
            return TransformedMultivariateParameter.class;
        }

        @Override
        public String getParserName() {
            return NAME;
        }
    };
}
