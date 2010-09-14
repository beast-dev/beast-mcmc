package dr.evomodel.MSSD;

import dr.evolution.tree.Tree;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.AbstractModelLikelihood;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;

/**
 * Package: CTMCScalePrior
 * Description:
 * <p/>
 * <p/>
 * Created by
 * Alexander V. Alekseyenko (alexander.alekseyenko@gmail.com)
 * Date: Aug 22, 2008
 * Time: 3:26:57 PM
 */
public class CTMCScalePrior extends AbstractModelLikelihood {
    Parameter ctmcScale;
    TreeModel treeModel;

    public CTMCScalePrior(String name, Parameter ctmcScale, TreeModel treeModel) {
        super(name);
        this.ctmcScale = ctmcScale;
        this.treeModel = treeModel;
    }

    protected void handleModelChangedEvent(Model model, Object object, int index) {
    }

    protected final void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
    }

    protected void storeState() {
    }

    protected void restoreState() {
    }

    protected void acceptState() {
    }

    public Model getModel() {
        return this;
    }

    public double getLogLikelihood() {
        double ab = ctmcScale.getParameterValue(0);
        double totalTreeTime = Tree.Utils.getTreeLength(treeModel, treeModel.getRoot());
        return -0.5 * Math.log(ab) - ab * totalTreeTime;
    }

    public void makeDirty() {
    }
}
