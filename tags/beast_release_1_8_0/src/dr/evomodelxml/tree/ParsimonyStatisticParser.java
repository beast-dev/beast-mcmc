package dr.evomodelxml.tree;

import dr.evolution.tree.Tree;
import dr.evolution.util.Taxa;
import dr.evolution.util.TaxonList;
import dr.evomodel.tree.ParsimonyStatistic;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Statistic;
import dr.xml.*;

/**
 */
public class ParsimonyStatisticParser extends AbstractXMLObjectParser {

    public static final String PARSIMONY_STATISTIC = "parsimonyStatistic";
    public static final String STATE = "state";

    public String getParserName() {
        return PARSIMONY_STATISTIC;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        String name = xo.getAttribute(Statistic.NAME, xo.getId());
        Tree tree = (Tree) xo.getChild(Tree.class);
        XMLObject cxo = xo.getChild(STATE);
        TaxonList taxa = (TaxonList) cxo.getChild(TaxonList.class);

        try {
            return new ParsimonyStatistic(name, tree, taxa);
        } catch (Tree.MissingTaxonException mte) {
            throw new XMLParseException("Taxon, " + mte + ", in " + getParserName() + "was not found in the tree.");
        }
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "A statistic that has as its value the parsimony tree length of a set of a " +
                "binary state defined by a set of taxa for a given tree";
    }

    public Class getReturnType() {
        return ParsimonyStatistic.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new StringAttributeRule(Statistic.NAME, "A name for this statistic primarily for the purposes of logging", true),
            new ElementRule(TreeModel.class),
            new ElementRule(STATE,
                    new XMLSyntaxRule[]{new ElementRule(Taxa.class)})
    };

}
