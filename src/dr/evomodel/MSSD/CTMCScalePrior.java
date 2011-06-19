package dr.evomodel.MSSD;

import dr.evolution.tree.Tree;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.AbstractModelLikelihood;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;

/**
 * @author Alexander V. Alekseyenko (alexander.alekseyenko@gmail.com)
 * @author Marc A. Suchard
 * 
 * Date: Aug 22, 2008
 * Time: 3:26:57 PM
 */
public class CTMCScalePrior extends AbstractModelLikelihood {
    final private Parameter ctmcScale;
    final private TreeModel treeModel;
    private double treeLength;
    private boolean treeLengthKnown;

    public CTMCScalePrior(String name, Parameter ctmcScale, TreeModel treeModel) {
        super(name);
        this.ctmcScale = ctmcScale;
        this.treeModel = treeModel;
        addModel(treeModel);
        treeLengthKnown = false;
    }

    private void updateTreeLength() {
        treeLength = Tree.Utils.getTreeLength(treeModel, treeModel.getRoot());
    }

    protected void handleModelChangedEvent(Model model, Object object, int index) {
        if (model == treeModel) {
            treeLengthKnown = false;
        }
    }

    protected final void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
    }

    protected void storeState() {
    }

    protected void restoreState() {
        treeLengthKnown = false;
    }

    protected void acceptState() {
    }

    public Model getModel() {
        return this;
    }

    public double getLogLikelihood() {
        double ab = ctmcScale.getParameterValue(0);
//        if (!treeLengthKnown) {
//            updateTreeLength();
//            treeLengthKnown = true;
//        }
        double totalTreeTime = Tree.Utils.getTreeLength(treeModel, treeModel.getRoot());
        return -0.5 * Math.log(ab) - ab * totalTreeTime; // TODO Change to treeLength and confirm results
    }

    public void makeDirty() {
        treeLengthKnown = false;
    }
}
