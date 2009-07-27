package dr.evomodel.continuous;

import dr.evolution.tree.Tree;
import dr.evolution.tree.NodeRef;
import dr.xml.*;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.tree.TreeStatistic;
import dr.inference.model.Statistic;

/**
 * @author Marc Suchard
 * @author Philippe Lemey
 */
public class TreeDispersionStatistic extends Statistic.Abstract implements TreeStatistic {

    public static final String TREE_DISPERSION_STATISTIC = "treeDispersionStatistic";
    public static final String BOOLEAN_OPTION = "booleanOption";

    public TreeDispersionStatistic(String name, TreeModel tree, SampledMultivariateTraitLikelihood traitLikelihood,
                                   boolean genericOption) {
        super(name);
        this.tree = tree;
        this.traitLikelihood = traitLikelihood;
        this.genericOption = genericOption;
    }

    public void setTree(Tree tree) {
        this.tree = (TreeModel) tree;
    }

    public Tree getTree() {
        return tree;
    }

    public int getDimension() {
        return 1;
    }

    /**
     * @return whatever Philippe wants
     */
    public double getStatisticValue(int dim) {

        String traitName = traitLikelihood.getTraitName();

        for (int i = 0; i < tree.getNodeCount(); i++) {
            NodeRef node = tree.getNode(i);

            double[] trait = tree.getMultivariateNodeTrait(node, traitName);

            if (node != tree.getRoot()) {

                double branchLength = tree.getBranchLength(node);
                // or
                double diffusionRate = traitLikelihood.getRescaledBranchLength(node);

                if (genericOption)
                    ;
            }
        }
        return 666;
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return TREE_DISPERSION_STATISTIC;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            String name = xo.getAttribute(NAME, xo.getId());
            TreeModel tree = (TreeModel) xo.getChild(Tree.class);

            boolean option = xo.getAttribute(BOOLEAN_OPTION,false); // Default value is false

            SampledMultivariateTraitLikelihood traitLikelihood = (SampledMultivariateTraitLikelihood)
                    xo.getChild(SampledMultivariateTraitLikelihood.class);

            return new TreeDispersionStatistic(name, tree, traitLikelihood, option);
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "A statistic that returns the average of the branch rates";
        }

        public Class getReturnType() {
            return TreeStatistic.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                AttributeRule.newStringRule(NAME, true),
                AttributeRule.newBooleanRule(BOOLEAN_OPTION,true),
                new ElementRule(TreeModel.class),
                new ElementRule(SampledMultivariateTraitLikelihood.class),
        };
    };

    private TreeModel tree = null;
    private boolean genericOption;
    private SampledMultivariateTraitLikelihood traitLikelihood;
}
