package dr.evomodelxml.treelikelihood.thorneytreelikelihood;

import dr.evomodel.treelikelihood.thorneytreelikelihood.BranchLengthProvider;
import dr.evomodel.treelikelihood.thorneytreelikelihood.StrictClockBranchLengthLikelihoodDelegate;
import dr.inference.model.Parameter;
import dr.xml.*;

public class StrictClockBranchLengthLikelihoodParser extends AbstractXMLObjectParser {

    public static final String STRICT_CLOCK_BRANCHLENGTH_LIKELIHOOD = "strictClockBranchLengthLikelihood";
    public static final String DATA_TREE = "dataTree";

    public String getParserName() {
        return STRICT_CLOCK_BRANCHLENGTH_LIKELIHOOD;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        double scale = xo.getAttribute("scale",1.0);

        Parameter mutationRate = (Parameter)xo.getChild(Parameter.class);

        return new StrictClockBranchLengthLikelihoodDelegate(STRICT_CLOCK_BRANCHLENGTH_LIKELIHOOD,mutationRate,scale);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element provides the likelihood of observing a branchlength given a mutation rate.";
    }

    public Class getReturnType() {
        return BranchLengthProvider.class;
    }

    public static final XMLSyntaxRule[] rules = {
            new ElementRule(Parameter.class),
            AttributeRule.newDoubleRule("scale",true,"a scale factor to muliply by the rate such as sequence length. default is 1"),

    };

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }
}
