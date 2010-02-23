package dr.inference.operators;

import dr.inference.model.Likelihood;
import dr.xml.*;

/**
 * @author Marc Suchard
 */
public class DirtyLikelihoodOperator extends SimpleMCMCOperator implements GibbsOperator {

    public static final String TOUCH_OPERATOR = "dirtyLikelihood";

    public DirtyLikelihoodOperator(Likelihood likelihood, double weight) {
        this.likelihood = likelihood;
        setWeight(weight);
    }

    public String getPerformanceSuggestion() {
        return null;
    }

    public String getOperatorName() {
        return TOUCH_OPERATOR;
    }

    public double doOperation() throws OperatorFailedException {
        likelihood.makeDirty();
        return 0;
    }

    public int getStepCount() {
        return 1;
    }

        public static dr.xml.XMLObjectParser PARSER = new dr.xml.AbstractXMLObjectParser() {

        public String getParserName() {
            return TOUCH_OPERATOR;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            double weight = xo.getDoubleAttribute(WEIGHT);

            Likelihood likelihood = (Likelihood) xo.getChild(Likelihood.class);

            return new DirtyLikelihoodOperator(likelihood,weight);

        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "This element returns a operator that forces the entire model likelihood recomputation";
        }

        public Class getReturnType() {
            return MCMCOperator.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                AttributeRule.newDoubleRule(WEIGHT),
                new ElementRule(Likelihood.class),
        };
    };

    private Likelihood likelihood;
}
