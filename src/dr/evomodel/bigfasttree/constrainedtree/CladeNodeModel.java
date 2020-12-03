package dr.evomodel.bigfasttree.constrainedtree;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeUtils;
import dr.evolution.util.Taxon;
import dr.evomodel.tree.TreeChangedEvent;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Variable;

import java.util.*;

/**
 * A  model class maintains a map between nodes an resolved tree and internal nodes in a static, compatible, potentially
 * unresolved tree.
 * <p>
 * NB: It is currently assumed that external nodes
 * never leave the clade. The model will only update internal nodes that may move about the tree while respecting the
 * clades.
 *
 * @author JT McCrone
 */


public class CladeNodeModel extends AbstractModel {
    public static final String CLADE_NODE_MODEL = "cladeNodeModel";

    public CladeNodeModel(String name, Tree cladeTree, TreeModel treeModel) throws TreeUtils.MissingTaxonException {
        super(name);
        if (treeModel.getExternalNodeCount() != cladeTree.getExternalNodeCount()) {
            System.err.println("Tree model and clade tree do not have the same number of taxa. Free taxon are not" +
                    "implemented. If you would like free taxa you can construct a ghost tree and pass it's" +
                    " corporeal tree here");
            // now find missing taxon
        }

        for (int i = 0; i < treeModel.getTaxonCount(); i++) {
            String id = treeModel.getTaxonId(i);
            if (cladeTree.getTaxonIndex(id) == -1) {
                throw new TreeUtils.MissingTaxonException(treeModel.getTaxon(i));
            }
        }
        this.cladeTree = cladeTree;
        this.treeModel = treeModel;
        // need store restore to trigger here so add as dummy likelihood
        // but also need to listen to the tree for changes

        addModel(treeModel);
        nodeCladeMap = new int[treeModel.getNodeCount()];
        storedNodeCladeMap = new int[treeModel.getNodeCount()];
        cladeRoots = new int[cladeTree.getInternalNodeCount()];
        storedCladeRoots = new int[cladeTree.getInternalNodeCount()];
        clades = new CladeRef[cladeTree.getInternalNodeCount()];

        cladesKnown = false;
        storedCladesKnown = false;

        updatedNode = new boolean[treeModel.getNodeCount()];
        storedUpdatedNode = new boolean[treeModel.getNodeCount()];

        Arrays.fill(updatedNode, true);
        Arrays.fill(storedUpdatedNode, true);

        // Get clades in cladeTree
        Set<BitSet> cladesInCladeTree = getClades(cladeTree, treeModel);

        Set<BitSet> cladesInTreeModel = getClades(treeModel, treeModel);

        boolean allFound = true;

        for (BitSet clade :
                cladesInCladeTree) {
            if (!cladesInTreeModel.contains(clade)) {
                allFound = false;
                break;
            }
        }
        if (!allFound) {
            throw new IllegalArgumentException("Tree model is not compatible with constraints tree");
        }


        // find matching clades in tree model
        setUpClades(cladeTree.getRoot(),null);
        updateClades();
        rootClade = clades[cladeTree.getRoot().getNumber()-cladeTree.getExternalNodeCount()];

    }

    public CladeNodeModel(Tree cladeTree, TreeModel treeModel) throws TreeUtils.MissingTaxonException {
        this(CLADE_NODE_MODEL, cladeTree, treeModel);
    }

    public TreeModel getTreeModel() {
        return treeModel;
    }

    public Tree getCladeTree() {
        return cladeTree;
    }

    public CladeRef getRootClade() {
        return rootClade;
    }

    public int getCladeCount() {
        return clades.length;
    }

    public int getChildCount(CladeRef clade) {
        Clade c = (Clade) clade;
        return c.getChildCount();
    }

    public CladeRef getChild(CladeRef clade, int i) {
        Clade c = (Clade) clade;
        return c.getChild(i);
    }

    public CladeRef getParent(CladeRef clade) {
        Clade c = (Clade) clade;
        return c.getParent();
    }

    public NodeRef getRootNode(CladeRef c) {
        if (!cladesKnown) {
            updateClades();
        }
        return treeModel.getNode(cladeRoots[c.getNumber()]);
    }

    public void setRootNode(CladeRef c, NodeRef n) {
//        cladeRoots[c.getNumber()] = n.getNumber();
        throw new UnsupportedOperationException("no longer in use");
    }

    private void setCladeRootNode(CladeRef c, NodeRef n) {
        cladeRoots[c.getNumber()] = n.getNumber();
    }

    public CladeRef getClade(NodeRef nodeRef) {
        if (!cladesKnown) {
            updateClades();
        }
        return getClade(nodeCladeMap[(nodeRef.getNumber())]);
    }

    public CladeRef getClade(int i) {
        return clades[i];
    }

    public List<NodeRef> getCladeRoots() {
        if (!cladesKnown) {
            updateClades();
        }
        List<NodeRef> roots = new ArrayList<>();
        for (int i = 0; i < clades.length; i++) {
            roots.add(getRootNode(clades[i]));
        }
        return roots;
    }

    private void setClade(NodeRef nodeRef, CladeRef cladeRef) {
        nodeCladeMap[nodeRef.getNumber()] = cladeRef.getNumber();
    }


    private void updateClades() {
        CladeRef rootClade = updateClades(treeModel.getRoot());
        setCladeRootNode(rootClade, treeModel.getRoot());
        cladesKnown = true;
        //TODO set make this the root clade
    }

    private CladeRef updateClades(NodeRef nodeRef) {

        if (updatedNode[nodeRef.getNumber()]) {
            if (!treeModel.isExternal(nodeRef)) {
                //assumes bifurcating tree
                CladeRef child0Clade = updateClades(treeModel.getChild(nodeRef, 0));
                CladeRef child1Clade = updateClades(treeModel.getChild(nodeRef, 1));

                if (child0Clade == child1Clade) {
                    //same clade as kids
                    setClade(nodeRef, child0Clade);
                } else {//Either one child clade is the parent of the other or they are siblings
                    if (child0Clade == getParent(child1Clade)) {
                        setClade(nodeRef, child0Clade);
                        setCladeRootNode(child1Clade, treeModel.getChild(nodeRef, 1));
                    } else if (child1Clade == getParent(child0Clade)) {
                        setClade(nodeRef, child1Clade);
                        setCladeRootNode(child0Clade, treeModel.getChild(nodeRef, 0));
                    } else {   // else this belongs to the parent of both clades and both kids are the child clade root
                        assert getParent(child0Clade) == getParent(child1Clade);
                        setClade(nodeRef, getParent(child0Clade));
                        setCladeRootNode(child0Clade, treeModel.getChild(nodeRef, 0));
                        setCladeRootNode(child1Clade, treeModel.getChild(nodeRef, 1));
                    }
                }
            }else{
                Taxon taxon = treeModel.getNodeTaxon(nodeRef);
                NodeRef nodeInCladeTree =  cladeTree.getExternalNode(cladeTree.getTaxonIndex(taxon.getId()));
                setClade(nodeRef,clades[cladeTree.getParent(nodeInCladeTree).getNumber()-cladeTree.getExternalNodeCount()]);
            }
        }
        updatedNode[nodeRef.getNumber()] = false;
        return clades[nodeCladeMap[(nodeRef.getNumber())]];
    }

    private void updateNode(NodeRef node) {
        updatedNode[node.getNumber()] = true;
        cladesKnown = false;
        NodeRef parent = treeModel.getParent(node);
        if (parent != null && !updatedNode[parent.getNumber()]) {
            updateNode(parent);
        }
    }

    private void updateAllNodes() {
        Arrays.fill(updatedNode, true);
        cladesKnown = false;
    }
    // *****************************************************************************************************************
    // ********                                 Setup helpers                                                   ********
    // *****************************************************************************************************************

    /**
     * Compiles the set of clades from the constraints tree using the numbers from the target tree.
     *
     * @param queryTree
     * @param referenceTree
     * @return
     */
    private Set<BitSet> getClades(Tree queryTree, Tree referenceTree) {
        Set<BitSet> clades = new HashSet<>();
        getClades(queryTree, queryTree.getRoot(), referenceTree, clades);
        return clades;
    }

    /**
     * recursive version
     *
     * @param queryTree
     * @param node
     * @param referenceTree
     * @param clades
     * @return
     */
    private BitSet getClades(Tree queryTree, NodeRef node, Tree referenceTree, Set<BitSet> clades) {
        BitSet clade = new BitSet();
        if (queryTree.isExternal(node)) {
            String taxonId = queryTree.getNodeTaxon(node).getId();
            clade.set(referenceTree.getTaxonIndex(taxonId));
        } else {
            for (int i = 0; i < queryTree.getChildCount(node); i++) {
                NodeRef child = queryTree.getChild(node, i);
                clade.or(getClades(queryTree, child, referenceTree, clades));
            }
            clades.add(clade);
        }
        return clade;
    }

    //pre-order traversal that makes a clade for each internal node of the clade defining tree
    private Clade setUpClades(NodeRef nodeRef,Clade parentClade) {
        int internalNodeIndex = nodeRef.getNumber() - cladeTree.getExternalNodeCount();
        Clade clade = new Clade(internalNodeIndex, parentClade);
        clades[internalNodeIndex] = clade;

        for (int i = 0; i < cladeTree.getChildCount(nodeRef); i++) {
            NodeRef child = cladeTree.getChild(nodeRef, i);
            if(!cladeTree.isExternal(child)){
                clade.addChild(setUpClades(child,clade));
            }
        }
        return clade;
    }


    // *********************************************************************************************
    // *****                             Model implementation                                *******
    // *********************************************************************************************

    public void makeDirty(){
        updateAllNodes();
    }
    @Override
    protected void handleModelChangedEvent(Model model, Object object, int i) {
        fireModelChanged();

        if (model == treeModel) {
            if (object instanceof TreeChangedEvent) {

                if (((TreeChangedEvent) object).isNodeChanged()) {
                    // If a node event occurs the node and its two child nodes
                    // are flagged for updating (this will result in everything
                    // above being updated as well. Node events occur when a node
                    // is added to a branch, removed from a branch or its height or
                    // rate changes.
                    NodeRef node = ((TreeChangedEvent) object).getNode();
                    //external nodes can't change clades
                    if(!treeModel.isExternal(node)){
                        updateNode(node);
                    }
                } else if (((TreeChangedEvent) object).isTreeChanged()) {
                    // Full tree events result in a complete updating of the tree likelihood
                    // This event type is now used for EmpiricalTreeDistributions.
//                    System.err.println("Full tree update event - these events currently aren't used\n" +
//                            "so either this is in error or a new feature is using them so remove this message.");
                    updateAllNodes();
                }
            }
        } else {
            throw new RuntimeException("Unknown componentChangedEvent");
        }
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int i, Variable.ChangeType changeType) {

    }

    @Override
    protected void storeState() {
        System.arraycopy(cladeRoots, 0, storedCladeRoots, 0, cladeRoots.length);
        System.arraycopy(nodeCladeMap, 0, storedNodeCladeMap, 0, nodeCladeMap.length);
        System.arraycopy(updatedNode, 0, storedUpdatedNode, 0, updatedNode.length);
        storedCladesKnown = cladesKnown;
    }

    @Override
    protected void restoreState() {
        int[] tmp = storedCladeRoots;
        storedCladeRoots = cladeRoots;
        cladeRoots = tmp;

        int[] tmp2 = storedNodeCladeMap;
        storedNodeCladeMap = nodeCladeMap;
        nodeCladeMap = tmp2;


        boolean[] tmp3 = storedUpdatedNode;
        storedUpdatedNode = updatedNode;
        updatedNode = tmp3;

        cladesKnown = storedCladesKnown;
    }

    @Override
    protected void acceptState() {

    }

    // Private Stuff

    private final CladeRef[] clades;
    private int[] nodeCladeMap;
    private final CladeRef rootClade;
    private int[] storedNodeCladeMap;
    private final TreeModel treeModel;
    private int[] cladeRoots;
    private int[] storedCladeRoots;
    private boolean[] updatedNode;
    private boolean[] storedUpdatedNode;
    private boolean cladesKnown;
    private boolean storedCladesKnown;
    private Tree cladeTree;


    /**
     * A private helper clade class.
     */
    class Clade implements CladeRef {
        //TODO figure out what in interface and what just here.
        public final static String CLADE = "clade";

        public Clade(int number, Clade parent) {
            this.number = number;
            this.parent = parent;
        }

        public void addChild(Clade clade) {
            this.children.add(clade);
        }


        public int getNumber() {
            return number;
        }

        public int getChildCount() {
            return children.size();
        }

        public Clade getChild(int i) {
            return children.get(i);
        }

        public Clade getParent() {
            return parent;
        }


        public final List<Clade> children = new ArrayList<>();
        private final Clade parent;
        private final int number;
    }
}