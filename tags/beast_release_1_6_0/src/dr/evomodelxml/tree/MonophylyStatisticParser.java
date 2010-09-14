package dr.evomodelxml.tree;

import dr.evolution.tree.Tree;
import dr.evolution.util.Taxa;
import dr.evolution.util.Taxon;
import dr.evolution.util.TaxonList;
import dr.evomodel.tree.MonophylyStatistic;
import dr.inference.model.Statistic;
import dr.xml.*;

/**
 */
public class MonophylyStatisticParser extends AbstractXMLObjectParser {

    public static final String MONOPHYLY_STATISTIC = "monophylyStatistic";
    public static final String MRCA = "mrca";
    public static final String IGNORE = "ignore";

    public String getParserName() {
        return MONOPHYLY_STATISTIC;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        String name = xo.getAttribute(Statistic.NAME, xo.getId());

        Tree tree = (Tree) xo.getChild(Tree.class);

        XMLObject cxo = xo.getChild(MRCA);
        TaxonList taxa = (TaxonList) cxo.getChild(TaxonList.class);
        if (taxa == null) {
            Taxa taxa1 = new Taxa();
            for (int i = 0; i < cxo.getChildCount(); i++) {
                Object ccxo = cxo.getChild(i);
                if (ccxo instanceof Taxon) {
                    taxa1.addTaxon((Taxon) ccxo);
                }
            }
            taxa = taxa1;
        }

        TaxonList ignore = null;
        if (xo.hasChildNamed(IGNORE)) {
            cxo = xo.getChild(IGNORE);
            ignore = (TaxonList) cxo.getChild(TaxonList.class);
            if (ignore == null) {
                Taxa taxa1 = new Taxa();
                for (int i = 0; i < cxo.getChildCount(); i++) {
                    Object ccxo = cxo.getChild(i);
                    if (ccxo instanceof Taxon) {
                        taxa1.addTaxon((Taxon) ccxo);
                    }
                }
                ignore = taxa1;
            }
        }

        try {
            return new MonophylyStatistic(name, tree, taxa, ignore);
        } catch (Tree.MissingTaxonException mte) {
            throw new XMLParseException("Taxon, " + mte + ", in " + getParserName() + "was not found in the tree.");
        }
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "A statistic that returns true if a given set of taxa are monophyletic for a given tree";
    }

    public Class getReturnType() {
        return MonophylyStatistic.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new StringAttributeRule(Statistic.NAME, "A name for this statistic for the purpose of logging", true),
            // Any tree will do, no need to insist on a Tree Model
            new ElementRule(Tree.class),
            new ElementRule(MRCA, new XMLSyntaxRule[]{
                    new XORRule(
                            new ElementRule(Taxon.class, 1, Integer.MAX_VALUE),
                            new ElementRule(Taxa.class)
                    )
            }),
            new ElementRule(IGNORE, new XMLSyntaxRule[]{
                    new XORRule(
                            new ElementRule(Taxon.class, 1, Integer.MAX_VALUE),
                            new ElementRule(Taxa.class)
                    )
            }, "An optional list of taxa to ignore from the test of monophyly", true)
    };

}
