package dr.evomodel.continuous;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.tree.TreeChangedEvent;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.treedatalikelihood.continuous.cdi.ContinuousDiffusionIntegrator;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;

import java.util.BitSet;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * @author Marc A. Suchard
 */
public class RestrictedPartialsModel extends AbstractModel {

    private final TreeModel treeModel;
    private final List<RestrictedPartials> restrictedPartialsList;
    private final ContinuousDiffusionIntegrator cdi;

    private boolean updateTreeMapping = true;
    private boolean updateRestrictedPartials = true;

    private final int startingPartialCount;
    private final int startingMatrixCount;

    private int partialsCount;
    private int matrixCount;
    private final int sparePartialIndex;

    /**
     * @param name Model Name
     */
    public RestrictedPartialsModel(String name, List<RestrictedPartials> restrictedPartialsList,
                                   int startingPartialBufferCount, int startingMatrixBufferCount) {
        super(name);

        this.treeModel = validateTreeModel(restrictedPartialsList);
        this.startingPartialCount = startingPartialBufferCount;
        this.startingMatrixCount = startingMatrixBufferCount;

        partialsCount = 0;
        matrixCount = 0;

        for (RestrictedPartials partials : restrictedPartialsList) {
            partials.setIndex(partialsCount + startingPartialBufferCount);
            addRestrictedPartials(partials);
            ++partialsCount;
        }
        sparePartialIndex = partialsCount + startingPartialBufferCount;
        ++partialsCount;

        setupClamps();

        this.restrictedPartialsList = restrictedPartialsList;
        this.cdi = null; // likelihoodDelegate.getIntegrator();

    }

    public int getExtraPartialBufferCount() { return partialsCount; }

    public int getExtraMatrixBufferCount() {  return matrixCount; }

    public void updatePartialRestrictions() {

        if (updateTreeMapping) {
            setupClamps();
            updateTreeMapping = false;
        }

        if (updateRestrictedPartials) {
            setupRestrictedPartials();
            updateRestrictedPartials = false;
        }
    }

    public boolean hasAnyPartialRestrictions() { return anyClamps; }

    private TreeModel validateTreeModel(List<RestrictedPartials> restrictedPartialsList) {
        TreeModel treeModel = restrictedPartialsList.get(0).getTreeModel();
        for (RestrictedPartials restrictedPartials : restrictedPartialsList) {
            if (restrictedPartials.getTreeModel() != treeModel) {
                throw new IllegalArgumentException("All tree models must be the same");
            }
        }
        return treeModel;
    }

    private void addRestrictedPartials(RestrictedPartials nodeClamp) {
        if (clampList == null) {
            clampList = new HashMap<BitSet, RestrictedPartials>();
        }

        clampList.put(nodeClamp.getTipBitSet(), nodeClamp);
        addModel(nodeClamp);

        System.err.println("Added a CLAMP!");
    }

    private void setupRestrictedPartials() {
        for (RestrictedPartials clamp : clampList.values()) {
            final int clampIndex = clamp.getIndex();
            final double[] partials = clamp.getRestrictedPartials();

            ContinuousDiffusionIntegrator cdi = null;
            cdi.setPostOrderPartial(clampIndex, partials); // TODO Double buffering?
        }

    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {

        if (model == treeModel) {
            if (object instanceof TreeChangedEvent) {
                TreeChangedEvent event = (TreeChangedEvent) object;
                if (event.isTreeChanged()) {
                    updateTreeMapping = true;
                } else if (event.isNodeChanged()) {
                    updateTreeMapping = true;
                }
            }
        } else if (model instanceof RestrictedPartials) {
            updateRestrictedPartials = true;
        } else {
            throw new RuntimeException("Unknown model");
        }

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
        updateTreeMapping = true; // TODO cache?  Caching is still not working, see IMTL.restoreState()
        updateRestrictedPartials = true;
    }

    /**
     * This call specifies that the current state is accept. Most models will not need to do anything.
     * Sub-models are handled automatically and do not need to be considered in this method.
     */
    @Override
    protected void acceptState() {

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

    public void setupClamps() {
        if (nodeToClampMap == null) {
            nodeToClampMap = new HashMap<NodeRef, RestrictedPartials>();
        }
        nodeToClampMap.clear();

        recursiveSetupClamp(treeModel, treeModel.getRoot(), new BitSet());

        anyClamps = (nodeToClampMap.size() > 0);
    }

    private void recursiveSetupClamp(Tree tree, NodeRef node, BitSet tips) {

        if (tree.isExternal(node)) {
            tips.set(node.getNumber());
        } else {
            for (int i = 0; i < tree.getChildCount(node); i++) {
                NodeRef child = tree.getChild(node, i);

                BitSet childTips = new BitSet();
                recursiveSetupClamp(tree, child, childTips);
                tips.or(childTips);
            }

            if (clampList.containsKey(tips)) {
                RestrictedPartials partials = clampList.get(tips);
                partials.setNode(node);
                nodeToClampMap.put(node, partials);
            }
        }
    }

    protected Map<BitSet, RestrictedPartials> clampList = null;
    protected Map<NodeRef, RestrictedPartials> nodeToClampMap = null;

    protected boolean anyClamps = false;

}
