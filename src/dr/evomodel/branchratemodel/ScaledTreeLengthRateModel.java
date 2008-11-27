package dr.evomodel.branchratemodel;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.ParameterChangeType;
import dr.xml.*;

import java.util.logging.Logger;

/**
 * Package: ScaledTreeLengthRateModel
 * Description:
 * <p/>
 * <p/>
 * Created by
 * Alexander V. Alekseyenko (alexander.alekseyenko@gmail.com)
 * Date: Mar 18, 2008
 * Time: 3:58:43 PM
 */
public class ScaledTreeLengthRateModel extends AbstractModel implements BranchRateModel {
    private Parameter totalLength;
    protected Tree treeModel;
    private double storedRateFactor;
    private boolean currentFactorKnown;
    private double rateFactor;

    ScaledTreeLengthRateModel(TreeModel treeModel, Parameter totalLength) {
        super("ScaledTreeLengthRateModel");
        this.totalLength = totalLength;
        this.treeModel = treeModel;
        currentFactorKnown = false;
        addModel(treeModel);
        addParameter(totalLength);
    }

    public double getBranchRate(Tree tree, NodeRef node) {
        if (!currentFactorKnown)
            updateCurrentLength();
        if (tree == treeModel) {
            return rateFactor;
        }
        return 0; // This is an error, we are referenced through a different Tree!!!
    }

    public double getTotalLength() {
        return totalLength.getParameterValue(0);
    }

    protected void updateCurrentLength() {
        double currentLength = 0;
        NodeRef root = treeModel.getRoot();
        for (int i = 0; i < treeModel.getNodeCount(); ++i) {
            NodeRef node = treeModel.getNode(i);
            if (node != root) {
                currentLength += treeModel.getBranchLength(node);
            }
        }
        rateFactor = totalLength.getParameterValue(0) / currentLength;
        currentFactorKnown = true;
    }

    protected void handleModelChangedEvent(Model model, Object object, int index) {
        if (model == treeModel) {
            currentFactorKnown = false;
        }
    }

    /**
     * This method is called whenever a parameter is changed.
     * <p/>
     * It is strongly recommended that the model component sets a "dirty" flag and does no
     * further calculations. Recalculation is typically done when the model component is asked for
     * some information that requires them. This mechanism is 'lazy' so that this method
     * can be safely called multiple times with minimal computational cost.
     */
    protected final void handleParameterChangedEvent(Parameter parameter, int index, ParameterChangeType type) {
        if (parameter == totalLength) {
            currentFactorKnown = false;
        }
    }

    /**
     * Additional state information, outside of the sub-model is stored by this call.
     */
    protected void storeState() {
        storedRateFactor = rateFactor;
    }

    /**
     * After this call the model is guaranteed to have returned its extra state information to
     * the values coinciding with the last storeState call.
     * Sub-models are handled automatically and do not need to be considered in this method.
     */
    protected void restoreState() {
        rateFactor = storedRateFactor;
    }

    /**
     * This call specifies that the current state is accept. Most models will not need to do anything.
     * Sub-models are handled automatically and do not need to be considered in this method.
     */
    protected void acceptState() {
    }

    public String getBranchAttributeLabel() {
        return "rate";
    }

    public String getAttributeForBranch(Tree tree, NodeRef node) {
        return Double.toString(getBranchRate(tree, node));
    }

    public static final String MODEL_NAME = "scaledTreeLengthModel";
    public static final String SCALING_FACTOR = "scalingFactor";
    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************
    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return MODEL_NAME;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            TreeModel tree = (TreeModel) xo.getChild(TreeModel.class);
            Parameter totalLength = (Parameter) xo.getElementFirstChild(SCALING_FACTOR);
            if (totalLength == null) {
                totalLength = new Parameter.Default(1, 1.0);
            }
            Logger.getLogger("dr.evomodel.branchratemodel").info("\n ---------------------------------\nCreating ScaledTreeLengthRateModel model.");
            Logger.getLogger("dr.evomodel.branchratemodel").info("\tTotal tree length will be scaled to " + totalLength.getParameterValue(0) + ".");
            Logger.getLogger("dr.evomodel.branchratemodel").info("\tIf you publish results using this rate model, please reference Alekseyenko, Lee and Suchard (in submision).\n---------------------------------\n");

            return new ScaledTreeLengthRateModel(tree, totalLength);
        }

        public String getParserDescription() {
            return
                    "This element returns a branch rate model that scales the total length of the tree to specified valued (default=1.0).";
        }

        public Class getReturnType() {
            return ScaledTreeLengthRateModel.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                new ElementRule(SCALING_FACTOR,
                        new XMLSyntaxRule[]{new ElementRule(Parameter.class)}, true),
                new ElementRule(TreeModel.class)
        };
    };
}
