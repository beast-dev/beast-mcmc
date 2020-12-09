package dr.evomodelxml.treelikelihood.thorneytreelikelihood;

import dr.evolution.tree.Tree;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.treelikelihood.thorneytreelikelihood.*;
import dr.inference.model.Likelihood;
import dr.xml.*;

public class ConstrainedBranchLengthProviderParser extends AbstractXMLObjectParser {

    public static final String CONSTRAINED_BRANCHLENGTH_PROVIDER = "constrainedBranchLengthProvider";
    public static final String DATA_TREE = "dataTree";

    public String getParserName() {
        return CONSTRAINED_BRANCHLENGTH_PROVIDER;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        double scale = xo.getAttribute("scale",1.0);
        boolean discrete = xo.getAttribute("discrete", true);
        double minBranchLength=xo.getAttribute("minBranchlength", 0.0);


        ConstrainedTreeModel constrainedTree = (ConstrainedTreeModel) xo.getChild(ConstrainedTreeModel.class);
        Tree dataTree = (Tree) xo.getElementFirstChild(DATA_TREE);

        return new ConstrainedTreeBranchLengthProvider(constrainedTree,dataTree,scale,minBranchLength,discrete);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element represents a mapping between branches in a  constrained tree and the number of mutations on the same branch in a fixed tree ";
    }

    public Class getReturnType() {
        return BranchLengthProvider.class;
    }

    public static final XMLSyntaxRule[] rules = {
            new ElementRule(ConstrainedTreeModel.class),
            new ElementRule(DATA_TREE, new XMLSyntaxRule[]{
                    new ElementRule(Tree.class)
            }),
            AttributeRule.newDoubleRule("scale",true,"a scale factor to muliply by the branchlengths in the data tree such as sequence length. default is 1"),
            AttributeRule.newDoubleRule("minBranchlength",true,"minimum branch length to be provided for branches not in data tree. default is 0.0"),
            AttributeRule.newBooleanRule("discrete",true,"should branchlengths be rounded to a discrete number of mutations? default is true")
    };

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }
}
