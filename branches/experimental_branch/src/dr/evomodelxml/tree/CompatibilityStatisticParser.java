package dr.evomodelxml.tree;

import dr.evolution.tree.Tree;
import dr.evomodel.tree.CompatibilityStatistic;
import dr.inference.model.Statistic;
import dr.xml.*;

/**
 */
public class CompatibilityStatisticParser extends AbstractXMLObjectParser {

    public static final String COMPATIBILITY_STATISTIC = "compatibilityStatistic";
    public static final String COMPATIBLE_WITH = "compatibleWith";

    public String getParserName() {
        return COMPATIBILITY_STATISTIC;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        String name = xo.getAttribute(Statistic.NAME, xo.getId());

        Tree tree1 = (Tree) xo.getChild(Tree.class);

        XMLObject cxo = xo.getChild(COMPATIBLE_WITH);
        Tree tree2 = (Tree) cxo.getChild(Tree.class);

        try {
            return new CompatibilityStatistic(name, tree1, tree2);
        } catch (Tree.MissingTaxonException mte) {
            throw new XMLParseException("Taxon, " + mte + ", in " + getParserName() + "was in the source tree but not the constraints tree.");
        }
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "A statistic that returns true if a pair of trees are compatible";
    }

    public Class getReturnType() {
        return CompatibilityStatistic.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new StringAttributeRule(Statistic.NAME, "A name for this statistic for the purpose of logging", true),
            new ElementRule(Tree.class),
            new ElementRule(COMPATIBLE_WITH, new XMLSyntaxRule[]{
                    new ElementRule(Tree.class)
            }),
            new ElementRule(Tree.class)
    };

}
