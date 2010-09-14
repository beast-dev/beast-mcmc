package dr.evomodelxml.MSSD;

import dr.evomodel.MSSD.ExponentialBranchLengthTreePrior;
import dr.evomodel.tree.TreeModel;
import dr.xml.*;

/**
 *
 */
public class ExponentialBranchLengthTreePriorParser extends AbstractXMLObjectParser {
    public static final String MODEL_NAME = "exponentialBranchLengthsPrior";
    
    public String getParserName() {
        return MODEL_NAME;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);

        return new ExponentialBranchLengthTreePrior(treeModel);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element represents a tree prior assuming " +
                "exponentially distributed branch lengths.";
    }

    public Class getReturnType() {
        return ExponentialBranchLengthTreePrior.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new ElementRule(TreeModel.class),
    };
}
