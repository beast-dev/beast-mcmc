package dr.evomodelxml.tree;

import dr.evolution.tree.Tree;
import dr.evolution.util.Taxa;
import dr.evolution.util.TaxonList;
import dr.evomodel.tree.ParsimonyStateStatistic;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Statistic;
import dr.xml.*;

/**
 */
public class ParsimonyStateStatisticParser extends AbstractXMLObjectParser {

    public static final String PARSIMONY_STATE_STATISTIC = "parsimonyStateStatistic";
    public static final String STATE = "state";
    public static final String MRCA = "mrca";

    public String getParserName() {
        return PARSIMONY_STATE_STATISTIC;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        String name = xo.getAttribute(Statistic.NAME, xo.getId());

        Tree tree = (Tree) xo.getChild(Tree.class);
        XMLObject cxo = xo.getChild(STATE);
        TaxonList stateTaxa = (TaxonList) cxo.getChild(TaxonList.class);

        cxo = xo.getChild(MRCA);
        TaxonList mrcaTaxa = (TaxonList) cxo.getChild(TaxonList.class);

        try {
            return new ParsimonyStateStatistic(name, tree, stateTaxa, mrcaTaxa);
        } catch (Tree.MissingTaxonException mte) {
            throw new XMLParseException("Taxon, " + mte + ", in " + getParserName() + "was not found in the tree.");
        }
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "A statistic that has as its value the parsimony state reconstruction of a " +
                "binary state defined by a set of taxa at a given MRCA of a tree";
    }

    public Class getReturnType() {
        return ParsimonyStateStatistic.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new StringAttributeRule(Statistic.NAME, "A name for this statistic for the purposes of logging", true),
            new ElementRule(TreeModel.class),
            new ElementRule(STATE,
                    new XMLSyntaxRule[]{new ElementRule(Taxa.class)}),
            new ElementRule(MRCA,
                    new XMLSyntaxRule[]{new ElementRule(Taxa.class)})
    };

}
