package dr.evomodelxml.tree;

import dr.evolution.tree.MutableTreeModel;
import dr.evolution.tree.NodeRef;
import dr.evolution.util.Taxa;
import dr.evolution.util.Taxon;
import dr.evomodel.continuous.AncestralTaxonInTree;
import dr.evomodel.tree.DuplicationTreeModel;
import dr.inference.model.Parameter;
import dr.xml.*;

import java.util.List;


/**
 * @author Xiang Ji
 */
public class DuplicationTreeModelParser extends AbstractXMLObjectParser {

    private static final String DUPLICATION_TREE_MODEL = "duplicationTreeModel";
    private static final String DUPLICATION = "duplication";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        MutableTreeModel tree = (MutableTreeModel) xo.getChild(MutableTreeModel.class);
        List<AncestralTaxonInTree> duplications = AncestralTraitTreeModelParser.parseAllAncestors(tree, xo, DUPLICATION);

        int index = tree.getExternalNodeCount();
        for (AncestralTaxonInTree duplication : duplications) {
            duplication.setIndex(index);
            duplication.setNode(new NodeRef() {
                @Override
                public int getNumber() {
                    return 0;
                }

                @Override
                public void setNumber(int n) {

                }
            });
        }

        return new DuplicationTreeModel(xo.getId(), tree, duplications);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules =
        new XMLSyntaxRule[]{
            new ElementRule(MutableTreeModel.class),
            new ElementRule(DUPLICATION, new XMLSyntaxRule[] {
                new ElementRule(Taxon.class),
                new ElementRule(Parameter.class, "Branch length towards duplication node."),
                new ElementRule(MonophylyStatisticParser.MRCA, Taxa.class),
            }, 1, 1),
        };

    @Override
    public String getParserDescription() {
        return null;
    }

    @Override
    public Class getReturnType() {
        return DuplicationTreeModel.class;
    }

    @Override
    public String getParserName() {
        return DUPLICATION_TREE_MODEL;
    }
}
