package dr.evomodel.branchratemodel;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.tree.TreeParameterModel;
import dr.evomodel.tree.randomlocalmodel.RandomLocalTreeVariable;
import dr.evomodelxml.branchratemodel.RelaxedDriftModelParser;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;

/**
 * Created by IntelliJ IDEA.
 * User: mandevgill
 * Date: 7/25/14
 * Time: 4:27 PM
 * To change this template use File | Settings | File Templates.
 */


public class RelaxedDriftModel extends AbstractBranchRateModel
        implements RandomLocalTreeVariable {

    public RelaxedDriftModel(TreeModel treeModel,
                             Parameter rateIndicatorParameter,
                             Parameter ratesParameter,
                             Parameter driftRates) {

        super(RelaxedDriftModelParser.RELAXED_DRIFT);

        rates = new TreeParameterModel(treeModel, ratesParameter, true);
        ratesParameter.addBounds(new Parameter.DefaultBounds(Double.MAX_VALUE, -Double.MAX_VALUE, ratesParameter.getDimension()));
        indicators = new TreeParameterModel(treeModel, rateIndicatorParameter, true);
        rateIndicatorParameter.addBounds(new Parameter.DefaultBounds(1, -1, rateIndicatorParameter.getDimension()));

        for (int i = 0; i < rateIndicatorParameter.getDimension(); i++) {
            rateIndicatorParameter.setParameterValue(i, 0.0);
        }
        for (int i = 0; i < ratesParameter.getDimension(); i++) {
            ratesParameter.setParameterValue(i, 0.0);
        }

        addModel(treeModel);
        this.treeModel = treeModel;
        addModel(indicators);
        addModel(rates);

        if (driftRates != null) {
            this.driftRates = driftRates;
            driftRates.setDimension(ratesParameter.getDimension());
        } else {
            driftRates = null;
        }

        branchRates = new double[treeModel.getNodeCount()];

        // Logger.getLogger("dr.evomodel").info("  indicator parameter name is '" + rateIndicatorParameter.getId() + "'");

        calculateBranchRates(treeModel);
        //recalculateScaleFactor();
    }


    /**
     * @param tree the tree
     * @param node the node to retrieve the variable of
     * @return the raw real-valued variable at this node
     */
    public final double getVariable(Tree tree, NodeRef node) {
        return rates.getNodeValue(tree, node);
    }

    /**
     * @param tree the tree
     * @param node the node
     * @return true of the variable at this node is included in function, thus representing a change in the
     * function looking down the tree.
     */
    public final boolean isVariableSelected(Tree tree, NodeRef node) {
        return indicators.getNodeValue(tree, node) != 0;
    }

    public void handleModelChangedEvent(Model model, Object object, int index) {
        recalculationNeeded = true;
        fireModelChanged();
    }

    protected final void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        recalculationNeeded = true;
        fireModelChanged();
    }

    protected void storeState() {
    }

    protected void restoreState() {
        calculateBranchRates(treeModel);
        // recalculateScaleFactor();
    }

    protected void acceptState() {
    }

    public double getBranchRate(final Tree tree, final NodeRef node) {
        if (recalculationNeeded) {
            calculateBranchRates(treeModel);
            // recalculateScaleFactor();
            recalculationNeeded = false;
        }
        return branchRates[node.getNumber()];
    }

    private void calculateBranchRates(TreeModel tree) {
        branchRates[tree.getRoot().getNumber()] = getVariable(tree, tree.getRoot());
        if (driftRates != null) {
            driftRates.setParameterValue(tree.getRoot().getNumber(), getVariable(tree, tree.getRoot()));
        }

        cbr(tree, tree.getRoot(), branchRates[tree.getRoot().getNumber()]);
    }

    /**
     * This is a recursive function that does the work of
     * calculating the unscaled branch rates across the tree
     * taking into account the indicator variables.
     *
     * @param tree the tree
     * @param node the node
     * @param rate the rate of the parent node
     */
    private void cbr(TreeModel tree, NodeRef node, double rate) {

        NodeRef childNode0 = tree.getChild(node, 0);
        NodeRef childNode1 = tree.getChild(node, 1);
        int nodeNumber0 = childNode0.getNumber();
        int nodeNumber1 = childNode1.getNumber();
        double nodeIndicator = indicators.getNodeValue(tree, node);
        if (indicators.getNodeValue(tree, node) != 0) {
            //     System.err.println("nodeIndicator: " + nodeIndicator);
        }
        if (nodeIndicator < 0) {
            //  System.err.println("child0 change");
            branchRates[nodeNumber0] = rate + getVariable(tree, childNode0);
            // branchRates[nodeNumber0] = getVariable(tree, childNode0);
            branchRates[nodeNumber1] = rate;

            if (driftRates != null) {
                driftRates.setParameterValue(nodeNumber0, rate + getVariable(tree, childNode0));
                driftRates.setParameterValue(nodeNumber1, rate);
            }

        } else if (nodeIndicator > 0) {
            //  System.err.println("child1 change");
            branchRates[nodeNumber0] = rate;
            branchRates[nodeNumber1] = rate + getVariable(tree, childNode1);
            // branchRates[nodeNumber1] = getVariable(tree, childNode1);

            if (driftRates != null) {
                driftRates.setParameterValue(nodeNumber0, rate);
                driftRates.setParameterValue(nodeNumber1, rate + getVariable(tree, childNode1));
            }

        } else {
            // System.err.println("NO CHANGES!!!");
            branchRates[nodeNumber0] = rate;
            branchRates[nodeNumber1] = rate;

            if (driftRates != null) {
                driftRates.setParameterValue(nodeNumber0, rate);
                driftRates.setParameterValue(nodeNumber1, rate);
            }
        }

        if (tree.getChildCount(childNode0) > 0) {
            cbr(tree, childNode0, branchRates[nodeNumber0]);
        }
        if (tree.getChildCount(childNode1) > 0) {
            cbr(tree, childNode1, branchRates[nodeNumber1]);
        }


        /*
        if (indicators.getNodeValue(tree, childNode0) > 0.5 && indicators.getNodeValue(tree, childNode1) < 0.5) {
            branchRates[nodeNumber0] = rate + getVariable(tree, childNode0);
           // branchRates[nodeNumber0] = getVariable(tree, childNode0);
            branchRates[nodeNumber1] = rate;
        } else if (indicators.getNodeValue(tree, childNode0) < 0.5 && indicators.getNodeValue(tree, childNode1) > 0.5) {
            branchRates[nodeNumber0] = rate;
            branchRates[nodeNumber1] = rate + getVariable(tree, childNode1);
           // branchRates[nodeNumber1] = getVariable(tree, childNode1);
        } else {
            branchRates[nodeNumber0] = rate;
            branchRates[nodeNumber1] = rate;
        }

        if (tree.getChildCount(childNode0) > 0) {
            cbr(tree, childNode0, branchRates[nodeNumber0]);
        }
        if (tree.getChildCount(childNode1) > 0) {
            cbr(tree, childNode1, branchRates[nodeNumber1]);
        }
        */
    }


    // the tree model
    private TreeModel treeModel;

    // the unscaled rates of each branch, taking into account the indicators
    private double[] branchRates;

    private TreeParameterModel indicators;
    private TreeParameterModel rates;
    private Parameter driftRates;

    boolean recalculationNeeded = true;
}