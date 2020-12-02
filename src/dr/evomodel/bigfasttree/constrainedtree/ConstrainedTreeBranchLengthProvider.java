package dr.evomodel.bigfasttree.constrainedtree;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.bigfasttree.BranchLengthProvider;
import dr.evomodel.bigfasttree.constrainedtree.CladeNodeModel;
import dr.evomodel.bigfasttree.constrainedtree.CladeRef;
import dr.inference.model.*;



//TODO cache branch lengths
public class ConstrainedTreeBranchLengthProvider  extends AbstractModel implements BranchLengthProvider {
    public static final String CONSTRAINED_TREE_BRANCHLENGTH_PROVIDER = "ConstrainedTreeBranchLengthProvider";
    public ConstrainedTreeBranchLengthProvider(CladeNodeModel cladeNodeModel){
        super(CONSTRAINED_TREE_BRANCHLENGTH_PROVIDER);

        this.cladeNodeModel=cladeNodeModel;
        addModel(cladeNodeModel);
        // set up branch lengths
        Tree dataTree = cladeNodeModel.getCladeTree();
        branchLengths = new double[dataTree.getNodeCount()];

        for (int i = 0; i < dataTree.getNodeCount(); i++) {
            double x = dataTree.getBranchLength(dataTree.getNode(i));
            branchLengths[i] = x;
        }
    }


    @Override
    public double getBranchLength(Tree tree, NodeRef node) {
        assert tree== cladeNodeModel.getTreeModel();
        if(tree.isExternal(node)){
            return branchLengths[node.getNumber()];
        }
        CladeRef clade = cladeNodeModel.getClade(node);
        if (cladeNodeModel.getRootNode(clade) == node) {
            return branchLengths[tree.getExternalNodeCount() -1 + clade.getNumber()];
            //get internalNde().getNumber
        }
        return 0d;
    }
    /**
     * Adds a listener that is notified when the this model changes.
     *
     * @param listener
     */
    @Override
    public void addModelListener(ModelListener listener) {

    }

    /**
     * Remove a listener previously addeed by addModelListener
     *
     * @param listener
     */
    @Override
    public void removeModelListener(ModelListener listener) {

    }

    /**
     * @return whether this model is in a valid state
     */
    @Override
    public boolean isValidState() {
        return false;
    }

    /**
     * @return the total number of sub-models
     */
    @Override
    public int getModelCount() {
        return 0;
    }



    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {

    }

    /**
     * This method is called whenever a parameter is changed.
     * <p/>
     * It is strongly recommended that the model component sets a "dirty" flag and does no
     * further calculations. Recalculation is typically done when the model component is asked for
     * some information that requires them. This mechanism is 'lazy' so that this method
     * can be safely called multiple times with minimal computational cost.
     *
     * @param variable
     * @param index
     * @param type
     */
    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {

    }


    /**
     * Additional state information, outside of the sub-model is stored by this call.
     */
    @Override
    protected void storeState() {

    }

    /**
     * After this call the model is guaranteed to have returned its extra state information to
     * the values coinciding with the last storeState call.
     * Sub-models are handled automatically and do not need to be considered in this method.
     */
    @Override
    protected void restoreState() {

    }

    /**
     * This call specifies that the current state is accept. Most models will not need to do anything.
     * Sub-models are handled automatically and do not need to be considered in this method.
     */
    @Override
    protected void acceptState() {

    }

    /**
     * is the model being listened to by another or by a likelihood?
     *
     * @return
     */
    @Override
    public boolean isUsed() {
        return false;
    }

    /**
     * @return the id as a string.
     */
    @Override
    public String getId() {
        return null;
    }

    /**
     * set the id as a string.
     *
     * @param id
     */
    @Override
    public void setId(String id) {

    }

    private final double[] branchLengths;
    private final CladeNodeModel cladeNodeModel;
}
