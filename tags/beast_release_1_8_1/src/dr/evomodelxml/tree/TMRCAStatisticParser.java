package dr.evomodelxml.tree;

import dr.evolution.tree.Tree;
import dr.evolution.util.Taxa;
import dr.evolution.util.TaxonList;
import dr.evomodel.tree.TMRCAStatistic;
import dr.inference.model.Statistic;
import dr.xml.*;

/**
 * @author Alexei Drummond
 * @author Andrew Rambaut
 */
public class TMRCAStatisticParser extends AbstractXMLObjectParser {

    public static final String TMRCA_STATISTIC = "tmrcaStatistic";
    public static final String MRCA = "mrca";
    // The tmrcaStatistic will represent that age of the parent node of the MRCA, rather than the MRCA itself
    public static final String PARENT = "forParent";
    public static final String STEM = "includeStem";


    public String getParserName() {
        return TMRCA_STATISTIC;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        String name = xo.getAttribute(Statistic.NAME, xo.getId());
        Tree tree = (Tree) xo.getChild(Tree.class);
        TaxonList taxa = (TaxonList) xo.getElementFirstChild(MRCA);
        boolean isRate = xo.getAttribute("rate", false);
        boolean includeStem = false;
        if (xo.hasAttribute(PARENT) && xo.hasAttribute(STEM)) {
             throw new XMLParseException("Please use either " + PARENT + " or " + STEM + "!");
        } else if (xo.hasAttribute(PARENT)) {
             includeStem = xo.getBooleanAttribute(PARENT);
        } else if (xo.hasAttribute(STEM)) {
             includeStem = xo.getBooleanAttribute(STEM);
        }

        try {
            return new TMRCAStatistic(name, tree, taxa, isRate, includeStem);
        } catch (Tree.MissingTaxonException mte) {
            throw new XMLParseException(
                    "Taxon, " + mte + ", in " + getParserName() + "was not found in the tree.");
        }
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "A statistic that has as its value the height of the most recent common ancestor " +
                "of a set of taxa in a given tree";
    }

    public Class getReturnType() {
        return TMRCAStatistic.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(Tree.class),
            new StringAttributeRule("name",
                    "A name for this statistic primarily for the purposes of logging", true),
            AttributeRule.newBooleanRule("rate", true),
            new ElementRule(MRCA,
                    new XMLSyntaxRule[]{new ElementRule(Taxa.class)}),
            new OrRule(
                    new XMLSyntaxRule[]{
                            AttributeRule.newBooleanRule(PARENT, true),
                            AttributeRule.newBooleanRule(STEM, true)
                    }
            )
    };

}
