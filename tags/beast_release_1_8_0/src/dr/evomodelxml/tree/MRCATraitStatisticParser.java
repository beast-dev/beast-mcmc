package dr.evomodelxml.tree;

import dr.evolution.tree.Tree;
import dr.evolution.util.Taxa;
import dr.evolution.util.TaxonList;
import dr.evomodel.tree.MRCATraitStatistic;
import dr.evomodel.tree.TreeModel;
import dr.xml.*;

/**
 */
public class MRCATraitStatisticParser extends AbstractXMLObjectParser {

    public static final String MRCA_TRAIT_STATISTIC = "mrcaTraitStatistic";
    public static final String MRCA = "mrca";
    public static final String NAME = "name";
    public static final String TRAIT = "trait";

    public String getParserName() {
        return MRCA_TRAIT_STATISTIC;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        String name = xo.getAttribute(NAME, xo.getId());
        String trait = xo.getStringAttribute(TRAIT);

        TreeModel tree = (TreeModel) xo.getChild(TreeModel.class);
        TaxonList taxa = (TaxonList) xo.getElementFirstChild(MRCA);

        try {
            return new MRCATraitStatistic(name, trait, tree, taxa);
        } catch (Tree.MissingTaxonException mte) {
            throw new XMLParseException("Taxon, " + mte + ", in " + getParserName() + "was not found in the tree.");
        }
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "A statistic that has as its value the height of the most recent common ancestor of a set of taxa in a given tree";
    }

    public Class getReturnType() {
        return MRCATraitStatistic.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new ElementRule(TreeModel.class),
            new StringAttributeRule("name", "A name for this statistic primarily for the purposes of logging", true),
            new StringAttributeRule("trait", "The name of the trait (can be rate)"),
            AttributeRule.newBooleanRule("rate", true),
            new ElementRule(MRCA,
                    new XMLSyntaxRule[]{new ElementRule(Taxa.class)})
    };

}
