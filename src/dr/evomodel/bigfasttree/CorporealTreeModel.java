package dr.evomodel.bigfasttree;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;

/**
 * A tree model that does not allow for the usual edits so that it's unlikely to be
 * edited outside of the ghost tree it shadows.
 *
 * @author JT McCrone
 */
public class CorporealTreeModel extends BigFastTreeModel {
    public static final String CORPOREAL_TREE_MODEL = "corporealTreeModel";

    public CorporealTreeModel(String name, Tree tree) {
        super(name, tree);
    }
    public CorporealTreeModel(Tree tree) {
        this(CORPOREAL_TREE_MODEL, tree);
    }

    // *****************************************************************
    // Interface MutableTree
    // *****************************************************************

    //The tree should only be edited by the ghosttree model. The usual edit methods throw errors
    /**
     * Set a new node as root node.
     */
    public void setRoot(NodeRef newRoot) {
        throw new UnsupportedOperationException("Function not available in CorporealTreeModel");
    }
    protected void makeRoot(NodeRef newRoot) {
        super.setRoot(newRoot);
    }

    public void addChild(NodeRef p, NodeRef c) {
        throw new UnsupportedOperationException("Function not available in CorporealTreeModel");
    }
    protected void adoptChild(NodeRef p, NodeRef c) {
        super.addChild(p, c);
    }

    public void removeChild(NodeRef p, NodeRef c) {
        throw new UnsupportedOperationException("Function not available in CorporealTreeModel");
    }
    protected void disownChild(NodeRef p, NodeRef c) {
        super.removeChild(p, c);
    }

    public boolean beginTreeEdit() {
        return super.beginTreeEdit();
    }

    public void endTreeEdit() {
        // and cleanup
        super.endTreeEdit();
    }

    public void setNodeHeight(NodeRef n, double height) {
        throw new UnsupportedOperationException("Function not available in CorporealTreeModel");
    }

    public void setNodeHeightQuietly(NodeRef n, double height) {
        throw new UnsupportedOperationException("Function not available in CorporealTreeModel");
    }
    protected void adjustNodeHeight(NodeRef n, double height) {
        super.setNodeHeight(n, height);
    }

    protected void adjustNodeHeightQuietly(NodeRef n, double height) {
        super.setNodeHeightQuietly(n, height);
    }

    public void setNodeRate(NodeRef n, double rate) {
        throw new UnsupportedOperationException("Function not available in CorporealTreeModel");
    }

    public void setNodeTrait(NodeRef n, String name, double value) {
        throw new UnsupportedOperationException("Function not available in CorporealTreeModel");
    }

    public void setMultivariateTrait(NodeRef n, String name, double[] value) {
        throw new UnsupportedOperationException("Function not available in CorporealTreeModel");
    }


}
