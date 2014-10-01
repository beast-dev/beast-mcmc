package dr.evomodelxml.tree;

import dr.evomodel.tree.TreeModel;
import dr.evomodel.tree.TreeShapeStatistic;
import dr.inference.model.Statistic;
import dr.xml.*;

/**
 */
public class TreeShapeStatisticParser extends AbstractXMLObjectParser {

    public static final String TREE_SHAPE_STATISTIC = "treeShapeStatistics";
    public static final String TARGET = "target";

    public String getParserName() {
        return TREE_SHAPE_STATISTIC;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        String name = xo.getAttribute(Statistic.NAME, xo.getId());
        TreeModel target = (TreeModel) xo.getElementFirstChild(TARGET);

        return new TreeShapeStatistic(name, target);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "A statistic that reports a handful of tree shape statistics on the given target tree.";
    }

    public Class getReturnType() {
        return TreeShapeStatistic.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new StringAttributeRule(Statistic.NAME, "A name for this statistic primarily for the purposes of logging", true),
            new ElementRule(TARGET,
                    new XMLSyntaxRule[]{new ElementRule(TreeModel.class)})
    };

}
