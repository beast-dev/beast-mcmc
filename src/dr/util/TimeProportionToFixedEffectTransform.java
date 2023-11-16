package dr.util;

import dr.xml.AbstractXMLObjectParser;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

public class TimeProportionToFixedEffectTransform extends Transform.MultivariateTransform {
    public static String NAME = "TimeProportionToFixedEffectTransform";

    // The variables are assumed to be in the following order:
    // 0: the proportion in time for forward transform, the fixed-effect for reverse transform
    // 1: the log-scale rate in the ancestral-model portion of the branch
    // 2: the log-scale rate in the descendant-model portion of the branch
    // The input is these as a length-3 vector, the output is a length-1 vector of the transformed value
    TimeProportionToFixedEffectTransform() {
        super(3, 1);
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
        double propTime = values[0];
        double rateAncestral = Math.exp(values[1]);
        double rateDescendant = Math.exp(values[2]);
        double[] transformed = new double[1];
        transformed[0] = Math.log(propTime * rateDescendant / rateAncestral + (1.0 - propTime));
        return transformed;
    }

    private void verifyOutput(double[] values) {
        double propTime = values[0];
        double rateAncestral = Math.exp(values[1]);
        double rateDescendant = Math.exp(values[2]);
        double[] transformed = new double[1];
        transformed[0] = Math.log(propTime * rateDescendant / rateAncestral + (1.0 - propTime));

        double fe = transformed[0];
        double lb = rateAncestral;
        double ub = rateDescendant;
        if (rateAncestral > rateDescendant) {
            lb = rateDescendant;
            ub = rateAncestral;
        }

        double newRate = rateAncestral * Math.exp(fe);
        assert(newRate >= lb && newRate <= ub);
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
        // Only the proportion is bounded
        return values[0] >= 0.0 && values[0] <= 1.0;
    }

    public static final AbstractXMLObjectParser PARSER = new AbstractXMLObjectParser() {

        @Override
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            final String name = xo.hasId() ? xo.getId() : null;

            TimeProportionToFixedEffectTransform transform = new TimeProportionToFixedEffectTransform();

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
