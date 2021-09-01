package dr.evomodelxml.treelikelihood.thorneytreelikelihood;

import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.treelikelihood.thorneytreelikelihood.PoissonBranchLengthLikelihoodDelegate;
import dr.inference.model.Parameter;
import dr.xml.*;

public class PoissonBranchLengthLikelihoodParser extends AbstractXMLObjectParser {

    public static final String STRICT_CLOCK_BRANCHLENGTH_LIKELIHOOD = "poissonBranchLengthLikelihood";
    public static final String DATA_TREE = "dataTree";

    public String getParserName() {
        return STRICT_CLOCK_BRANCHLENGTH_LIKELIHOOD;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        double scale = xo.getAttribute("scale",1.0);

        BranchRateModel branchRateModel = (BranchRateModel)xo.getChild(BranchRateModel.class);

        return new PoissonBranchLengthLikelihoodDelegate(STRICT_CLOCK_BRANCHLENGTH_LIKELIHOOD,branchRateModel,scale);
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
            new ElementRule(BranchRateModel.class),
            AttributeRule.newDoubleRule("scale",true,"a scale factor to muliply by the rate such as sequence length. default is 1"),

    };

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }
}
