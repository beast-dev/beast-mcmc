package dr.evomodelxml.treelikelihood;

import dr.evolution.tree.Tree;
import dr.evolution.util.Taxa;
import dr.evolution.util.TaxonList;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.treelikelihood.AncestralState;
import dr.evomodel.treelikelihood.AncestralStateTreeLikelihood;
import dr.xml.*;

/**
 */
public class AncestralStateParser extends AbstractXMLObjectParser {

    public static final String ANCESTRAL_STATE = "ancestralState";
    public static final String NAME = "name";
    public static final String MRCA = "mrca";

    public String getParserName() {
        return ANCESTRAL_STATE;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        String name = xo.getAttribute(NAME, xo.getId());
        Tree tree = (Tree) xo.getChild(Tree.class);
        TaxonList taxa = (TaxonList) xo.getElementFirstChild(MRCA);
        AncestralStateTreeLikelihood ancestralTreeLikelihood = (AncestralStateTreeLikelihood) xo.getChild(AncestralStateTreeLikelihood.class);
        try {
            return new AncestralState(name, ancestralTreeLikelihood, tree, taxa);
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
        return AncestralState.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new StringAttributeRule("name", "A name for this statistic primarily for the purposes of logging", true),
            new ElementRule(TreeModel.class),
            new ElementRule(AncestralStateTreeLikelihood.class),
            new ElementRule(MRCA,
                    new XMLSyntaxRule[]{new ElementRule(Taxa.class)})
    };
}
