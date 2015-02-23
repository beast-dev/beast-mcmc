package dr.evomodel.tree;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeNodeFilter;
import dr.evolution.util.TaxonList;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Variable;

import java.util.BitSet;
import java.util.HashSet;
import java.util.Set;

/**
 * @author Marc A. Suchard
 * @author Andrew Rambaut
 */
public class BackboneNodeFilter extends AbstractModel implements TreeNodeFilter {

    public BackboneNodeFilter(String name, Tree tree, TaxonList taxonList, boolean excludeClade, boolean includeStem) {
        super(name);

        this.tree = tree;
        this.excludeClade = excludeClade;
        this.includeStem = includeStem;

        backboneSet = new HashSet<NodeRef>();

        try {
            tipSet = Tree.Utils.getTipsBitSetForTaxa(this.tree, taxonList);
        } catch (Tree.MissingTaxonException e) {
            e.printStackTrace();
        }

        if (tree instanceof TreeModel) {
            addModel((TreeModel) tree);
        }

        backboneSetKnown = false;
    }

    public boolean includeNode(Tree tree, NodeRef node) {

        if (!backboneSetKnown) {
            computeBackBoneMap();
            backboneSetKnown = true;
        }
        return backboneSet.contains(node);
    }

    public void computeBackBoneMap() {
        backboneSet.clear();
        recursivelyComputeBackBoneMap(tree, tree.getRoot());
    }

    private boolean recursivelyComputeBackBoneMap(Tree tree, NodeRef node) {
       boolean onBackBone = false;

        if (tree.isExternal(node)) {
            if (tipSet.get(node.getNumber())) {
                onBackBone = true;
            }
        } else {
            for (int i = 0; i < tree.getChildCount(node); i++) {
                NodeRef child = tree.getChild(node, i);
                if (recursivelyComputeBackBoneMap(tree, child)) {
                    // if any of the desendents are back bone then this node is too
                    onBackBone = true;
                }
            }
        }

        if (onBackBone) {
            recursivelyPruneStemAndClade(tree, node, includeStem, excludeClade);
            return true;
        }
        return false;
    }

    private void recursivelyPruneStemAndClade(Tree tree, NodeRef node,
                                              boolean includeStem, boolean excludeClade) {

        if (!tree.isExternal(node) && !excludeClade) {
            for (int i = 0; i < tree.getChildCount(node); i++) {
                NodeRef child = tree.getChild(node, i);
                recursivelyPruneStemAndClade(tree, child, true, true);
            }
        }

        if (includeStem && !backboneSet.contains(node)) {
            backboneSet.add(node);
        }
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {

        if (model == tree) {
            backboneSetKnown = false;
        }

//            if (model == tree) {
//                if (object instanceof TreeModel.TreeChangedEvent) {
//                    TreeModel.TreeChangedEvent event = (TreeModel.TreeChangedEvent) object;
//                    if (event.isHeightChanged()) {
//                        // Do nothing
//                        return;
//                    }
//                    if (event.isNodeParameterChanged()) {
//                        // Do nothing
//                        return;
//                    } // else a topology change (?)
//                    backboneSetKnown = false;
//                    return;
//                }
//            }
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Variable.ChangeType type) {
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

    private boolean backboneSetKnown;
    private final Set<NodeRef> backboneSet;
    private BitSet tipSet;

    private final Tree tree;
    private final boolean includeStem;
    private final boolean excludeClade;
}