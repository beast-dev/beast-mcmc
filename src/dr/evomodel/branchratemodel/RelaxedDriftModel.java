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
                             boolean randInheritance,
                             boolean randNumChanges) {

        super(RelaxedDriftModelParser.RELAXED_DRIFT);


        //indicators = new TreeParameterModel(treeModel, rateIndicatorParameter, true);
        rates = new TreeParameterModel(treeModel, ratesParameter, true);
        randomInheritance = randInheritance;
        randomNumChanges = randNumChanges;

        // rateIndicatorParameter.addBounds(new Parameter.DefaultBounds(1, 0, rateIndicatorParameter.getDimension()));
        ratesParameter.addBounds(new Parameter.DefaultBounds(Double.MAX_VALUE, -Double.MAX_VALUE, ratesParameter.getDimension()));

        if (rateIndicatorParameter != null) {
            indicators = new TreeParameterModel(treeModel, rateIndicatorParameter, true);
            rateIndicatorParameter.addBounds(new Parameter.DefaultBounds(1, 0, rateIndicatorParameter.getDimension()));
            addModel(indicators);
            for (int i = 0; i < rateIndicatorParameter.getDimension(); i++) {
                rateIndicatorParameter.setParameterValue(i, 0.0);
            }
        }

        for (int i = 0; i < ratesParameter.getDimension(); i++) {
            ratesParameter.setParameterValue(i, 0.0);
        }

        // rootDrift = rootDriftParameter;

        addModel(treeModel);
        this.treeModel = treeModel;

        // addModel(indicators);
        addModel(rates);


        unscaledBranchRates = new double[treeModel.getNodeCount()];

        // Logger.getLogger("dr.evomodel").info("  indicator parameter name is '" + rateIndicatorParameter.getId() + "'");

        recalculateScaleFactor();
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
        return indicators.getNodeValue(tree, node) > 0.5;
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
        recalculateScaleFactor();
    }

    protected void acceptState() {
    }

    public double getBranchRate(final Tree tree, final NodeRef node) {
        if (recalculationNeeded) {
            recalculateScaleFactor();
            recalculationNeeded = false;
        }
        return unscaledBranchRates[node.getNumber()] * scaleFactor;
    }

    private void calculateUnscaledBranchRates(TreeModel tree) {
        unscaledBranchRates[tree.getRoot().getNumber()] = getVariable(tree, tree.getRoot());
        cubr(tree, tree.getRoot(), unscaledBranchRates[tree.getRoot().getNumber()]);
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
    private void cubr(TreeModel tree, NodeRef node, double rate) {

        NodeRef childNode0 = tree.getChild(node, 0);
        NodeRef childNode1 = tree.getChild(node, 1);
        int nodeNumber0 = childNode0.getNumber();
        int nodeNumber1 = childNode1.getNumber();


        if (randomInheritance) {


            if (randomNumChanges) {
                if (indicators.getNodeValue(tree, childNode0) > 0.5 && indicators.getNodeValue(tree, childNode1) < 0.5) {
                    unscaledBranchRates[nodeNumber0] = rate + getVariable(tree, childNode0);
                    unscaledBranchRates[nodeNumber1] = rate;
                } else if (indicators.getNodeValue(tree, childNode0) < 0.5 && indicators.getNodeValue(tree, childNode1) > 0.5) {
                    unscaledBranchRates[nodeNumber0] = rate;
                    unscaledBranchRates[nodeNumber1] = rate + getVariable(tree, childNode1);
                } else {
                    unscaledBranchRates[nodeNumber0] = rate;
                    unscaledBranchRates[nodeNumber1] = rate;
                }

            } else {
                if (indicators.getNodeValue(tree, node) > 0.5) {
                    unscaledBranchRates[nodeNumber0] = rate + getVariable(tree, childNode0);
                    unscaledBranchRates[nodeNumber1] = rate;
                } else {
                    unscaledBranchRates[nodeNumber0] = rate;
                    unscaledBranchRates[nodeNumber1] = rate + getVariable(tree, childNode1);
                }

            }

        } else {
            //  System.err.println("I get heeeeeere");
            unscaledBranchRates[nodeNumber0] = rate + getVariable(tree, childNode0);
            unscaledBranchRates[nodeNumber1] = rate;
        }

        if (tree.getChildCount(childNode0) > 0) {
            cubr(tree, childNode0, unscaledBranchRates[nodeNumber0]);
        }
        if (tree.getChildCount(childNode1) > 0) {
            cubr(tree, childNode1, unscaledBranchRates[nodeNumber1]);
        }


    }

    private void recalculateScaleFactor() {
        calculateUnscaledBranchRates(treeModel);
        scaleFactor = 1;
    }

    // AR - as TreeParameterModels are now loggable, the indicator parameter should be logged
    // directly.
//    private static String[] attributeLabel = {"changed"};
//
//    public String[] getNodeAttributeLabel() {
//        return attributeLabel;
//    }
//
//    public String[] getAttributeForNode(Tree tree, NodeRef node) {
//
//        if (tree.isRoot(node)) {
//            return new String[]{"false"};
//        }
//
//        return new String[]{(isVariableSelected((TreeModel) tree, node) ? "true" : "false")};
//    }

    // the scale factor necessary to maintain the mean rate
    private double scaleFactor;

    // the tree model
    private TreeModel treeModel;

    // true if the rate variables are treated as relative
    // to the parent rate rather than absolute rates
    private boolean randomInheritance = false;
    private boolean randomNumChanges = false;

    // the unscaled rates of each branch, taking into account the indicators
    private double[] unscaledBranchRates;

    // the mean rate across all the tree, if null then mean rate is scaled to 1.0
    private Parameter meanRateParameter;

    private TreeParameterModel indicators;
    private TreeParameterModel rates;
    private Parameter rootDrift;
    private double myRandom;

    boolean recalculationNeeded = true;
}