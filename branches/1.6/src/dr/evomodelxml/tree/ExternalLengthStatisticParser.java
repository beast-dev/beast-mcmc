package dr.evomodelxml.tree;

import dr.evolution.tree.Tree;
import dr.evolution.util.Taxa;
import dr.evolution.util.TaxonList;
import dr.evomodel.tree.ExternalLengthStatistic;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Statistic;
import dr.xml.*;

/**
 */
public class ExternalLengthStatisticParser extends AbstractXMLObjectParser {

    public static final String EXTERNAL_LENGTH_STATISTIC = "externalLengthStatistic";

    public String getParserName() {
        return EXTERNAL_LENGTH_STATISTIC;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        String name = xo.getAttribute(Statistic.NAME, xo.getId());

        Tree tree = (Tree) xo.getChild(Tree.class);
        TaxonList taxa = (TaxonList) xo.getChild(Taxa.class);

        try {
            return new ExternalLengthStatistic(name, tree, taxa);
        } catch (Tree.MissingTaxonException mte) {
            throw new XMLParseException("Taxon, " + mte + ", in " + getParserName() + "was not found in the tree.");
        }
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "A statistic that has as its value(s) the length of the external branch length(s) of a set of one or more taxa in a given tree";
    }

    public Class getReturnType() {
        return ExternalLengthStatistic.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new ElementRule(TreeModel.class),
            new StringAttributeRule("name", "A name for this statistic primarily for the purposes of logging", true),
            new ElementRule(Taxa.class)
    };

}
