package dr.evomodelxml.tree;

import dr.evolution.tree.Tree;
import dr.evomodel.tree.TreeHeightStatistic;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Statistic;
import dr.xml.*;

/**
 * Joseph Heled
 */
public class TreeHeightStatisticParser extends AbstractXMLObjectParser {

    public static final String TREE_HEIGHT_STATISTIC = "treeHeightStatistic";

        public String getParserName() {
            return TREE_HEIGHT_STATISTIC;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            final String name = xo.getAttribute(Statistic.NAME, xo.getId());
            final Tree tree = (Tree) xo.getChild(Tree.class);

            return new TreeHeightStatistic(name, tree);
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "A statistic that returns the height of the tree";
        }

        public Class getReturnType() {
            return TreeHeightStatistic.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                AttributeRule.newStringRule(Statistic.NAME, true),
                new ElementRule(TreeModel.class),
        };
}