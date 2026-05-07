
package dr.evomodelxml.treedatalikelihood;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.evolution.tree.TreeTraitProvider;
import dr.evolution.util.Taxon;
import dr.evomodel.treedatalikelihood.*;
import dr.evomodel.treedatalikelihood.discrete.beastBasedDiscreteTreeLikelihood.PreOrderMessageProvider;
import dr.evomodel.treedatalikelihood.preorder.AbstractBeagleGradientDelegate;
import dr.evomodel.treedatalikelihood.preorder.DiscretePreOrderType;
import dr.evomodel.treedatalikelihood.preorder.ProcessSimulationDelegate;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.ElementRule;
import dr.xml.Reportable;
import dr.xml.XMLObject;
import dr.xml.XMLObjectParser;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

import java.util.Arrays;

/**
 * Reportable XML helper for dumping all cached discrete pre-order partials.
 */
public class DiscretePreOrderReportParser extends AbstractXMLObjectParser {

    public static final String PARSER_NAME = "discretePreOrderReport";

    @Override
    public String getParserName() {
        return PARSER_NAME;
    }

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        TreeDataLikelihood likelihood = (TreeDataLikelihood) xo.getChild(TreeDataLikelihood.class);
        if (likelihood == null) {
            throw new XMLParseException("Expected a treeDataLikelihood child in " + PARSER_NAME);
        }
        return new DiscretePreOrderReport(likelihood);
    }

    @Override
    public String getParserDescription() {
        return "Reports branch-start and branch-end pre-order partials for a discrete TreeDataLikelihood.";
    }

    @Override
    public Class getReturnType() {
        return DiscretePreOrderReport.class;
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{
                new ElementRule(TreeDataLikelihood.class)
        };
    }

    public static final XMLObjectParser PARSER = new DiscretePreOrderReportParser();

    public static final class DiscretePreOrderReport implements Reportable {

        private final TreeDataLikelihood treeDataLikelihood;
        private final PreOrderMessageProvider discreteDelegate;
        private final TreeTrait treeTrait;

        private DiscretePreOrderReport(TreeDataLikelihood treeDataLikelihood) {
            this.treeDataLikelihood = treeDataLikelihood;

            DataLikelihoodDelegate delegate = treeDataLikelihood.getDataLikelihoodDelegate();
            if (delegate instanceof PreOrderMessageProvider) {
                this.discreteDelegate = (PreOrderMessageProvider) delegate;
                this.treeTrait = null;
            } else if (delegate instanceof BeagleDataLikelihoodDelegate) {

                BeagleDataLikelihoodDelegate likelihoodDelegate = (BeagleDataLikelihoodDelegate) delegate;
                String name = "preorder.partials";

                if (treeDataLikelihood.getTreeTrait(name) == null) {

                    ProcessSimulationDelegate preOrderDelegate = new AbstractBeagleGradientDelegate(
                            "test", treeDataLikelihood.getTree(), likelihoodDelegate) {
                        @Override
                        protected int getGradientLength() {
                            return 0;
                        }

                        @Override
                        protected void getNodeDerivatives(Tree tree, double[] first, double[] second) {

                        }

                        @Override
                        protected void constructTraits(Helper treeTraitHelper) {

                            treeTraitHelper.addTrait(new TreeTrait.DA() {
                                @Override
                                public String getTraitName() {
                                    return name;
                                }

                                @Override
                                public Intent getIntent() {
                                    return Intent.WHOLE_TREE;
                                }

                                @Override
                                public double[] getTrait(Tree tree, NodeRef node) {
                                    return getGradient(node);
                                }
                            });
                        }
                    };

                    TreeTraitProvider traitProvider = new ProcessSimulation(treeDataLikelihood, preOrderDelegate);
                    treeDataLikelihood.addTraits(traitProvider.getTreeTraits());

                    this.discreteDelegate = (PreOrderMessageProvider) preOrderDelegate;
                    this.treeTrait = treeDataLikelihood.getTreeTrait("preorder.partials");

                } else {
                    throw new RuntimeException("Unknown error");
                }
            } else {
                throw new IllegalArgumentException("Unknown likelihood");
            }
        }

        @Override
        public String getReport() {

            treeDataLikelihood.getLogLikelihood();
            Tree tree = treeDataLikelihood.getTree();

            if (treeTrait != null) {
                treeTrait.getTrait(tree, null);
            }

            int categoryCount = discreteDelegate.getCategoryCount();
            int patternCount = discreteDelegate.getPatternCount();
            int stateCount = discreteDelegate.getStateCount();

            double[] allStart = new double[stateCount * categoryCount * patternCount];
            double[] allEnd = new double[stateCount * categoryCount * patternCount];

            double[] start = new double[stateCount];
            double[] end = new double[stateCount];

            StringBuilder sb = new StringBuilder();
            sb.append("Discrete pre-order partials\n");
            sb.append("  treeDataLikelihood = ").append(treeDataLikelihood.getId()).append('\n');
            sb.append("  nodes = ").append(tree.getNodeCount()).append('\n');
            sb.append("  patterns = ").append(patternCount).append('\n');
            sb.append("  categories = ").append(categoryCount).append('\n');
            sb.append("  states = ").append(stateCount).append('\n');

            for (int i = 0; i < tree.getNodeCount(); i++) {
                NodeRef node = tree.getNode(i);
                int nodeNumber = node.getNumber();
                sb.append("node ").append(nodeNumber);
                if (tree.isRoot(node)) {
                    sb.append(" root");
                }
                if (tree.isExternal(node)) {
                    Taxon taxon = tree.getNodeTaxon(node);
                    sb.append(" taxon=").append(taxon == null ? "null" : taxon.getId());
                }
                sb.append('\n');

                discreteDelegate.getPreorderPartials(nodeNumber, DiscretePreOrderType.TOP, allStart);
                discreteDelegate.getPreorderPartials(nodeNumber, DiscretePreOrderType.BOTTOM, allEnd);

                for (int c = 0; c < categoryCount; c++) {
                    for (int p = 0; p < patternCount; p++) {

                        int offset = (c * patternCount + p) * stateCount;
                        System.arraycopy(allStart, offset, start, 0, stateCount);
                        System.arraycopy(allEnd, offset, end, 0, stateCount);

                        sb.append("  category ").append(c)
                                .append(" pattern ").append(p)
                                .append(" start=").append(Arrays.toString(start))
                                .append(" end=").append(Arrays.toString(end))
                                .append('\n');
                    }
                }
            }

            return sb.toString();
        }
    }
}

