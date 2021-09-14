//
// Source code recreated from a .class file by IntelliJ IDEA
// (powered by FernFlower decompiler)
//

package dr.evomodelxml.treelikelihood.thorneytreelikelihood;

import dr.evolution.tree.Tree;
import dr.evomodel.treelikelihood.thorneytreelikelihood.BranchLengthProvider;
import dr.evomodel.treelikelihood.thorneytreelikelihood.FixedTreeBranchLengthProvider;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.AttributeRule;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

public class FixedTreeBranchLengthProviderParser extends AbstractXMLObjectParser {
    public static final String FIXED_BRANCHLENGTH_PROVIDER = "fixedBranchLengthProvider";
    public static final String DATA_TREE = "dataTree";
    public static final XMLSyntaxRule[] rules = new XMLSyntaxRule[]{new ElementRule(Tree.class), new ElementRule("dataTree", new XMLSyntaxRule[]{new ElementRule(Tree.class)}), AttributeRule.newDoubleRule("scale", true, "a scale factor to muliply by the branchlengths in the data tree such as sequence length. default is 1"), AttributeRule.newDoubleRule("minBranchlength", true, "minimum branch length to be provided for branches not in data tree. default is 0.0"), AttributeRule.newBooleanRule("discrete", true, "should branchlengths be rounded to a discrete number of mutations? default is true")};

    public FixedTreeBranchLengthProviderParser() {
    }

    public String getParserName() {
        return "fixedBranchLengthProvider";
    }

    public Object parseXMLObject(XMLObject var1) throws XMLParseException {
        double scale =var1.getAttribute("scale", 1.0);
        boolean discrete = var1.getAttribute("discrete", true);
        double minBranchlength = var1.getAttribute("minBranchlength", 0.0);
        Tree fixedTree = (Tree)var1.getChild(Tree.class);
        Tree dataTree = (Tree)var1.getElementFirstChild("dataTree");
        return new FixedTreeBranchLengthProvider(fixedTree, dataTree, scale, minBranchlength, discrete);
    }

    public String getParserDescription() {
        return "This element represents a mapping between branches in a topologically fixed tree and the number of mutations on the same branch in a data tree ";
    }

    public Class getReturnType() {
        return BranchLengthProvider.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }
}
