package dr.evomodelxml.treelikelihood.thorneytreelikelihood;

import dr.evomodel.treelikelihood.thorneytreelikelihood.AURCBranchLengthLikelihoodDelegate;
import dr.evomodel.treelikelihood.thorneytreelikelihood.PoissonBranchLengthLikelihoodDelegate;
import dr.inference.model.Parameter;
import dr.xml.*;

public class AURCBranchLengthLikelihoodParser extends AbstractXMLObjectParser {

    public static final String ADDITIVE_UNCORRELATED_RELAXED_CLOCK = "AURCBranchLengthLikelihood";
    public static final String MU = "mu";
    public static final String OMEGA = "omega";

    public String getParserName() {
        return ADDITIVE_UNCORRELATED_RELAXED_CLOCK;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        double scale = xo.getAttribute("scale",1.0);

        Parameter mu = (Parameter) xo.getElementFirstChild(MU);
        Parameter omega = (Parameter) xo.getElementFirstChild(OMEGA);


        return new AURCBranchLengthLikelihoodDelegate(ADDITIVE_UNCORRELATED_RELAXED_CLOCK,mu,omega,scale);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element provides the likelihood of observing a branchlength assuming the rate of evolution follows a gamma distribution.";
    }

    public Class getReturnType() {
        return PoissonBranchLengthLikelihoodDelegate.class;
    }

    public static final XMLSyntaxRule[] rules = {
            new ElementRule(MU, Parameter.class, "mean clock rate  ", false),
            new ElementRule(OMEGA, Parameter.class, "relaxation parameter", false),

            AttributeRule.newDoubleRule("scale",true,"a scale factor to muliply by the rate such as sequence length. default is 1"),

    };

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }
}