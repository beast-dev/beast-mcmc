package dr.evomodelxml.branchratemodel;

import dr.evolution.tree.NodeRef;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.tree.TreeModel;
import dr.inference.loggers.LogColumn;
import dr.inference.loggers.Loggable;
import dr.inference.loggers.NumberColumn;
import dr.xml.*;

/**
 * @author Alexander Fisher
 */

public class LogRatesFromBranchRateModel implements Loggable {

    private static final String BRANCH_RATES = "ratesFromBranchRateModel";

    private BranchRateModel model;
    private TreeModel tree;
    private int numBranches;

    private LogRatesFromBranchRateModel(BranchRateModel model, TreeModel tree){
        this.model = model;
        this.tree = tree;
        this.numBranches = tree.getNodeCount() - 1;
    }

    @Override
    public LogColumn[] getColumns() {
        LogColumn[] logs = new LogColumn[numBranches];
//        double[] branchRates = new double[numBranches];
        NodeRef node;
        for (int i = 0; i < numBranches; i++) {
            node = tree.getNode(i);
            logs[i] = new RateColumn(getName(i), i, node) ;
//                    model.getBranchRate(tree, node);
        }
        return logs;
    }


    private String getName(int dim) {
        return "rate." + (dim + 1);
    }

    private class RateColumn extends NumberColumn {
        private final int dim;
        private final NodeRef node;
        public RateColumn(String label, int dim, NodeRef node) {
            super(label);
            this.dim = dim;
            this.node = node;
        }

        public double getDoubleValue() {
            return model.getBranchRate(tree, node); }
    }


    // **************************************************************
    // XMLObjectParser
    // **************************************************************

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return BRANCH_RATES;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            BranchRateModel model = (BranchRateModel) xo.getChild(BranchRateModel.class);
            TreeModel tree = (TreeModel) xo.getChild(TreeModel.class);

            return new LogRatesFromBranchRateModel(model, tree);
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                new ElementRule(BranchRateModel.class),
                new ElementRule(TreeModel.class),

        };

        public String getParserDescription() {
            return "Logs rates using getBranchRate from a branchRateModel";
        }

        public Class getReturnType() {
            return LogRatesFromBranchRateModel.class;
        }
    };
}

