package dr.evomodel.continuous;

import dr.evolution.tree.Tree;
import dr.evolution.tree.NodeRef;
import dr.xml.*;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.tree.TreeStatistic;
import dr.inference.model.Statistic;
import dr.geo.math.SphericalPolarCoordinates;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Marc Suchard
 * @author Philippe Lemey
 * @author Andrew Rambaut
 */
public class DiffusionRateStatistic extends Statistic.Abstract {

    public static final String DIFFUSION_RATE_STATISTIC = "diffusionRateStatistic";
    public static final String TREE_DISPERSION_STATISTIC = "treeDispersionStatistic";
    public static final String BOOLEAN_OPTION = "greatCircleDistance";

    public DiffusionRateStatistic(String name, TreeModel tree, List<AbstractMultivariateTraitLikelihood> traitLikelihoods,
                                   boolean genericOption) {
        super(name);
        this.traitLikelihoods = traitLikelihoods;
        this.genericOption = genericOption;
    }

    public int getDimension() {
        return 1;
    }

    /**
     * @return whatever Philippe wants
     */
    public double getStatisticValue(int dim) {

        String traitName = traitLikelihoods.get(0).getTraitName();
        double treelength = 0;
        double treeDistance = 0;

        for (AbstractMultivariateTraitLikelihood traitLikelihood : traitLikelihoods) {
            TreeModel tree = traitLikelihood.getTreeModel();

            for (int i = 0; i < tree.getNodeCount(); i++) {
                NodeRef node = tree.getNode(i);
                double[] trait = traitLikelihood.getTraitForNode(tree, node, traitName);

                if (node != tree.getRoot()) {

                    double[] parentTrait = traitLikelihood.getTraitForNode(tree, tree.getParent(node), traitName);
                    treelength += tree.getBranchLength(node);

                    if (genericOption) { // Great Circle distance
                        SphericalPolarCoordinates coord1 = new SphericalPolarCoordinates(trait[0], trait[1]);
                        SphericalPolarCoordinates coord2 = new SphericalPolarCoordinates(parentTrait[0], parentTrait[1]);
                        treeDistance += coord1.distance(coord2);
                    } else {
                        treeDistance += getNativeDistance(trait, parentTrait);
                    }

                }
            }
        }
        return treeDistance / treelength;
    }

    private double getNativeDistance(double[] location1, double[] location2) {
        return Math.sqrt(Math.pow((location2[0] - location1[0]), 2.0) + Math.pow((location2[1] - location1[1]), 2.0));
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return DIFFUSION_RATE_STATISTIC;
        }

        @Override
        public String[] getParserNames() {
            return new String[] { getParserName(), TREE_DISPERSION_STATISTIC };
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            String name = xo.getAttribute(NAME, xo.getId());
            TreeModel tree = (TreeModel) xo.getChild(Tree.class);

            boolean option = xo.getAttribute(BOOLEAN_OPTION, false); // Default value is false

            List<AbstractMultivariateTraitLikelihood> traitLikelihoods = new ArrayList<AbstractMultivariateTraitLikelihood>();

            for (int i = 0; i < xo.getChildCount(); i++) {
                if (xo.getChild(i) instanceof AbstractMultivariateTraitLikelihood) {
                     traitLikelihoods.add((AbstractMultivariateTraitLikelihood) xo.getChild(i));
                }
            }

            return new DiffusionRateStatistic(name, tree, traitLikelihoods, option);
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
                AttributeRule.newBooleanRule(BOOLEAN_OPTION, true),
                new ElementRule(AbstractMultivariateTraitLikelihood.class, 1, Integer.MAX_VALUE),
        };
    };

    private boolean genericOption;
    private List<AbstractMultivariateTraitLikelihood> traitLikelihoods;
}
