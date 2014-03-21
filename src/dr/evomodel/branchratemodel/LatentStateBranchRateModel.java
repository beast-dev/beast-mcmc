package dr.evomodel.branchratemodel;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.tree.TreeParameterModel;
import dr.inference.model.AbstractModelLikelihood;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;

/**
 * ${CLASS_NAME}
 *
 * @author Andrew Rambaut
 * @version $Id$
 *
 * $HeadURL$
 *
 * $LastChangedBy$
 * $LastChangedDate$
 * $LastChangedRevision$
 */
public class LatentStateBranchRateModel extends AbstractModelLikelihood implements BranchRateModel {

    private final TreeModel tree;
    private final Parameter latentTransitionRateParameter;
    private final Parameter latentStateProportionParameter;
    private final TreeParameterModel latentStateProportions;


    /**
     * @param name Model Name
     */
    public LatentStateBranchRateModel(String name,
                                      TreeModel treeModel,
                                      Parameter latentTransitionRateParameter) {
        super(name);

        this.tree = treeModel;
        addModel(tree);

        this.latentTransitionRateParameter = latentTransitionRateParameter;
        addVariable(latentTransitionRateParameter);

        this.latentStateProportionParameter = new Parameter.Default(0.5);
        this.latentStateProportions = new TreeParameterModel(tree, latentStateProportionParameter, false);

    }

    @Override
    public double getBranchRate(Tree tree, NodeRef node) {
        double latentProportion = latentStateProportions.getNodeValue(tree, node);
        double length = tree.getBranchLength(node);

        double rate = 1.0;

        // double rate = calculateRateHere(length, latentProportion);

        return rate;
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
        if (model == tree) {

        }

    }

    @Override
    protected void storeState() {

    }

    @Override
    protected void restoreState() {

    }

    @Override
    protected void acceptState() {

    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {

    }

    @Override
    public Model getModel() {
        return null;
    }

    @Override
    public double getLogLikelihood() {
        return 0;
    }

    @Override
    public void makeDirty() {

    }

    @Override
    public String getTraitName() {
        return null;
    }

    @Override
    public Intent getIntent() {
        return null;
    }

    @Override
    public Class getTraitClass() {
        return null;
    }

    @Override
    public Double getTrait(Tree tree, NodeRef node) {
        return null;
    }

    @Override
    public String getTraitString(Tree tree, NodeRef node) {
        return null;
    }

    @Override
    public boolean getLoggable() {
        return false;
    }

    @Override
    public TreeTrait[] getTreeTraits() {
        return new TreeTrait[0];
    }

    @Override
    public TreeTrait getTreeTrait(String key) {
        return null;
    }
}
