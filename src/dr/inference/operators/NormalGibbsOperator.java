package dr.inference.operators;

import dr.inference.distribution.MultivariateDistributionLikelihood;
import dr.inference.model.Parameter;
import dr.math.MathUtils;
import dr.xml.*;

public class NormalGibbsOperator extends SimpleMCMCOperator implements GibbsOperator {

    private final Parameter parameter;
    private final Parameter mean;
    private final Parameter prec;

    NormalGibbsOperator(Parameter parameter, Parameter mean, Parameter precision) {
        this.parameter = parameter;
        this.mean = mean;
        this.prec = precision;
    }

    @Override
    public String getOperatorName() {
        return NORMAL_GIBBS_OPERATOR;
    }

    @Override
    public double doOperation() {
        int dim = parameter.getDimension();
        for (int i = 0; i < dim; i++) {
            double draw = MathUtils.nextGaussian() / Math.sqrt(prec.getParameterValue(i)) + mean.getParameterValue(i);
            parameter.setParameterValueQuietly(i, draw);
        }
        parameter.fireParameterChangedEvent();
        return 0;
    }

    public static final String NORMAL_GIBBS_OPERATOR = "normalGibbsOperator";


    public static final AbstractXMLObjectParser PARSER = new AbstractXMLObjectParser() {

        private static final String DATA = "data";
        private static final String PRECISION = "precision";
        private static final String MEAN = "mean";

        @Override
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            XMLObject dxo = xo.getChild(DATA);
            Parameter data = (Parameter) dxo.getChild(Parameter.class);

            XMLObject pxo = xo.getChild(PRECISION);
            Parameter prec = (Parameter) pxo.getChild(Parameter.class);

            XMLObject mxo = xo.getChild(MEAN);
            Parameter mean = (Parameter) mxo.getChild(Parameter.class);

            return new NormalGibbsOperator(data, mean, prec);
        }

        @Override
        public XMLSyntaxRule[] getSyntaxRules() {
            return new XMLSyntaxRule[]{
                    new ElementRule(DATA, new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class)
                    }),
                    new ElementRule(PRECISION, new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class)
                    }),
                    new ElementRule(MEAN, new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class)
                    })
            };
        }

        @Override
        public String getParserDescription() {
            return "Gibbs operator that samples from a normal distribution.";
        }

        @Override
        public Class getReturnType() {
            return NormalGibbsOperator.class;
        }

        @Override
        public String getParserName() {
            return NORMAL_GIBBS_OPERATOR;
        }
    };

}
