package dr.evomodelxml.tree;

import dr.evolution.tree.Tree;
import dr.evomodel.tree.NodeHeightsStatistic;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Parameter;
import dr.inference.model.Statistic;
import dr.xml.*;

/**
 */
public class NodeHeightsStatisticParser extends AbstractXMLObjectParser {

    public static final String NODE_HEIGHTS_STATISTIC = "nodeHeightsStatistic";

    public String getParserName() {
        return NODE_HEIGHTS_STATISTIC;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        String name = xo.getAttribute(Statistic.NAME, xo.getId());
        Tree tree = (Tree) xo.getChild(Tree.class);

        Parameter groupSizes = (Parameter) xo.getChild(Parameter.class);

        return new NodeHeightsStatistic(name, tree, groupSizes);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "A statistic that returns the heights of each internal node in increasing order (or groups them by a group size parameter)";
    }

    public Class getReturnType() {
        return NodeHeightsStatistic.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newStringRule(Statistic.NAME, true),
            new ElementRule(TreeModel.class),
            new ElementRule(Parameter.class, true),
    };

}
