package dr.evomodelxml.bigfasttree.thorney;

import dr.evomodel.bigfasttree.thorney.PoissonBranchLengthLikelihoodDelegate;
import dr.xml.*;

public class PoissonBranchLengthLikelihoodParser extends AbstractXMLObjectParser {

    public static final String POISSON_BRANCH_LENGTH_LIKELIHOOD = "poissonBranchLengthLikelihood";

    public String getParserName() {
        return POISSON_BRANCH_LENGTH_LIKELIHOOD;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        double scale = xo.getAttribute("scale",1.0);
        return new PoissonBranchLengthLikelihoodDelegate(POISSON_BRANCH_LENGTH_LIKELIHOOD,scale);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element provides the likelihood of observing a branchlength given a mutation rate.";
    }

    public Class getReturnType() {
        return PoissonBranchLengthLikelihoodDelegate.class;
    }

    public static final XMLSyntaxRule[] rules = {
            AttributeRule.newDoubleRule("scale",true,"a scale factor to multiply by the rate such as sequence length. default is 1"),

    };

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }
}
