package dr.evomodel.branchratemodel;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxon;
import dr.evomodel.branchratemodel.AbstractBranchRateModel;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Statistic;
import dr.inference.model.Variable;
import dr.math.MathUtils;

import java.util.Iterator;
import java.util.List;

/**
 * @author Andrew Rambaut
 * @version $Id$
 */
public class AttributeBranchRateModel extends AbstractBranchRateModel {

    public AttributeBranchRateModel(final TreeModel treeModel, final String rateAttributeName) {
        super(ATTRIBUTE_BRANCH_RATE_MODEL);

        this.treeModel = treeModel;
        this.rateAttributeName = rateAttributeName;

        addModel(treeModel);
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
        fireModelChanged();
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        // nothing to do
    }

    @Override
    protected void storeState() {
        // nothing to do
    }

    @Override
    protected void restoreState() {
        // nothing to do
    }

    @Override
    protected void acceptState() {
        // nothing to do
    }

    @Override
    public double getBranchRate(Tree tree, NodeRef node) {
        Object value = tree.getNodeAttribute(node, rateAttributeName);
        return Double.parseDouble((String)value);
    }

    @Override
    public String getTraitName() {
        return rateAttributeName;
    }

    public static final String ATTRIBUTE_BRANCH_RATE_MODEL = "attributeBranchRateModel";

    private final TreeModel treeModel;
    private final String rateAttributeName;

}
