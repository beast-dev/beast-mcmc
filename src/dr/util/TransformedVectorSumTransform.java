package dr.util;

import dr.inference.hmc.GradientWrtIncrement;
import dr.inference.hmc.IncrementTransformType;
import dr.inference.model.Parameter;
import dr.inference.model.TransformedMultivariateParameter;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

public class TransformedVectorSumTransform extends Transform.MultivariateTransform {

    private static final String NAME = "transformedVectorSumTransform";
    private final IncrementTransformType type;
    private final int dim;

    public TransformedVectorSumTransform(int dim, IncrementTransformType type) {
        super(dim);
        this.dim = dim;
        this.type = type;
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
        return type.parameterFromIncrements(values);
    }

    @Override
    protected double[] inverse(double[] values) {
        return type.incrementsFromParameter(values);
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

            String ttype = (String) xo.getAttribute(INCREMENT_TRANSFORM);
            double upper = Double.POSITIVE_INFINITY;
            double lower = Double.NEGATIVE_INFINITY;

            if( xo.hasAttribute("upper") && xo.hasAttribute("lower")) {
                upper = xo.getDoubleAttribute("upper");
                lower = xo.getDoubleAttribute("lower");
            }

            // Get upper and lower from xml
            IncrementTransformType type = IncrementTransformType.factory(ttype, upper, lower);

            TransformedVectorSumTransform transform = new TransformedVectorSumTransform(dim, type);

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
