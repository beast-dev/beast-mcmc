package dr.evomodelxml.treelikelihood.thorneytreelikelihood;


import dr.evolution.util.TaxonList;
import dr.evomodel.treelikelihood.thorneytreelikelihood.ConstrainedTreeModel;
import dr.evomodel.treelikelihood.thorneytreelikelihood.SubtreeRootHeightStatistic;
import dr.inference.model.Statistic;
import dr.xml.*;

public class SubtreeRootHeightStatisticParser extends AbstractXMLObjectParser {
    public final static String SUBTREE_ROOT_HEIGHT_STATISTIC = "subtreeRootHeightStatistic";
    public static final String ABSOLUTE = "absolute";
    public static final String MRCA = "mrca";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        ConstrainedTreeModel tree = (ConstrainedTreeModel) xo.getChild(ConstrainedTreeModel.class);

        TaxonList taxa = null;

        if (xo.hasChildNamed(MRCA)) {
            taxa = (TaxonList) xo.getElementFirstChild(MRCA);
        }
        boolean isAbsolute = xo.getAttribute(ABSOLUTE, false);
        String name = xo.getAttribute(Statistic.NAME, xo.getId());
        return new SubtreeRootHeightStatistic(name,tree,taxa,isAbsolute);
    }

    /**
     * @return an array of syntax rules required by this element.
     * Order is not important.
     */
    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{
                new ElementRule(ConstrainedTreeModel.class),
        };    }

    @Override
    public String getParserDescription() {
        return "A parser for logging the heights of subtrees in a constrained tree model";
    }

    @Override
    public Class getReturnType() {
        return SubtreeRootHeightStatistic.class;
    }

    /**
     * @return Parser name, which is identical to name of xml element parsed by it.
     */
    @Override
    public String getParserName() {
        return SUBTREE_ROOT_HEIGHT_STATISTIC;
    }
}