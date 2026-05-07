
package dr.evomodelxml.treedatalikelihood;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxon;
import dr.evomodel.treedatalikelihood.DataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.DiscreteDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
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

        private DiscretePreOrderReport(TreeDataLikelihood treeDataLikelihood) {
            this.treeDataLikelihood = treeDataLikelihood;
        }

        @Override
        public String getReport() {
            DataLikelihoodDelegate delegate = treeDataLikelihood.getDataLikelihoodDelegate();
            if (!(delegate instanceof DiscreteDataLikelihoodDelegate)) {
                throw new IllegalStateException(PARSER_NAME + " requires a DiscreteDataLikelihoodDelegate");
            }

            DiscreteDataLikelihoodDelegate discreteDelegate = (DiscreteDataLikelihoodDelegate) delegate;
            if (!discreteDelegate.hasPreOrder()) {
                throw new IllegalStateException(PARSER_NAME + " requires usePreOrder=\"true\"");
            }

            treeDataLikelihood.getLogLikelihood();

            Tree tree = treeDataLikelihood.getTree();
            int categoryCount = discreteDelegate.getCategoryCount();
            int patternCount = discreteDelegate.getPatternCount();
            int stateCount = discreteDelegate.getStateCount();
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

                for (int c = 0; c < categoryCount; c++) {
                    for (int p = 0; p < patternCount; p++) {
                        discreteDelegate.getPreOrderBranchTopInto(nodeNumber, c, p, start);
                        discreteDelegate.getPreOrderBranchBottomInto(nodeNumber, c, p, end);
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

