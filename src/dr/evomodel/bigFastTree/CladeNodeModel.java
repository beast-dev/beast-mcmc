package dr.evomodel.bigFastTree;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeUtils;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Variable;

import java.util.*;
import java.util.stream.Collectors;

/**
 * A  model class that watches a tree model.
 * Things it does
 * 1) It keeps track of pre-specified nodes in the tree,
 * provides access to those nodes for efficient operation and statistics on them.
 *
 * @author JT McCrone
 */

public class CladeNodeModel  extends AbstractModel {
    public static final String CLADE_NODE_MODEL = "cladeNodeModel";

    public CladeNodeModel(String name, Tree cladeTree, TreeModel treeModel) throws TreeUtils.MissingTaxonException {
        super(name);
        if(treeModel.getExternalNodeCount()!=cladeTree.getExternalNodeCount()){
            System.err.println("Tree model and clade tree do not have the same number of taxa. Free taxon are not" +
                    "implemented yet.");
            // now find missing taxon
        }
        for (int i = 0; i < treeModel.getTaxonCount(); i++) {
            String id = treeModel.getTaxonId(i);
            if (cladeTree.getTaxonIndex(id) == -1) {
                throw new TreeUtils.MissingTaxonException(treeModel.getTaxon(i));
            }
        }

        this.treeModel = treeModel;
        treeModel.addModel(this);
//        addModel(treeModel); could make smart enough to update its self instead of having operators do it

        nodeCladeMap = new Clade[treeModel.getNodeCount()];
        cladeRoots = new int[cladeTree.getInternalNodeCount()];
        storedCladeRoots = new int[cladeTree.getInternalNodeCount()];

        // Get clades in cladeTree
        Set<BitSet> cladesInCladeTree = getClades(cladeTree,treeModel);

        Set<BitSet> cladesInTreeModel = getClades(treeModel, treeModel);

        boolean allFound = true;

        for (BitSet clade :
                cladesInCladeTree) {
            if(!cladesInTreeModel.contains(clade)){
                allFound=false;
                break;
            }
        }
        if(!allFound){
            throw new IllegalArgumentException("Tree model is not compatible with constraints tree");
        }
        // find matching clades in tree model
        setupClades(cladesInCladeTree,treeModel);
    }

    public CladeNodeModel(Tree cladeTree, TreeModel treeModel) throws TreeUtils.MissingTaxonException {
        this(CLADE_NODE_MODEL,cladeTree,treeModel);
    }

    public TreeModel getTreeModel() {
        return treeModel;
    }

    public int getCladeCount(){
        return cladeList.size();
    }


    public int getChildCount(CladeRef clade) {
        Clade c = (Clade) clade;
        return  c.getChildCount();
    }

    public CladeRef getChild(CladeRef clade, int i) {
        Clade c = (Clade) clade;
        return c.getChild(i);
    }
    public CladeRef getParent(CladeRef clade) {
        Clade c = (Clade) clade;
        return c.getParent();
    }

    public int getNodeCount(CladeRef c){
        Clade clade = (Clade) c;
        return clade.nodes.size() + clade.getChildCount();
    }

    public List<NodeRef> getNodes(CladeRef c){
        Clade clade = (Clade) c;
        return clade.getNodes().stream().map(treeModel::getNode).collect(Collectors.toList());
    }
    public NodeRef getRootNode(CladeRef c) {
        return treeModel.getNode( cladeRoots[c.getNumber()] );
    }
    public void setRootNode(CladeRef c,NodeRef n) {
        cladeRoots[c.getNumber()] = n.getNumber();
    }

    public NodeRef getNode(CladeRef c, int i) {
        Clade clade = (Clade) c;
         return treeModel.getNode((clade.getNode(i)));
    }

    public CladeRef getClade(NodeRef nodeRef) {
        return nodeCladeMap[(nodeRef.getNumber())];
    }
    public CladeRef getClade(int i) {
        return cladeList.get(i);
    }
    public List<NodeRef> getCladeRoots() {
        return cladeList.stream().map(this::getRootNode).collect(Collectors.toList());
    }
    //TODO set up isCompatible
    // each clade does a dfs for it's tips - in the bit and it's childs' root
    // flag a clade as updated if one of it's nodes updates.


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
        getClades(queryTree,queryTree.getRoot(),referenceTree,clades) ;
        return clades;
    }

    /**
     * recursive version
     * @param queryTree
     * @param node
     * @param referenceTree
     * @param clades
     * @return
     */
    private BitSet getClades(Tree queryTree, NodeRef node,Tree referenceTree,Set<BitSet> clades){
        BitSet clade = new BitSet();
        if (queryTree.isExternal(node)) {
            String taxonId = queryTree.getNodeTaxon(node).getId();
            clade.set(referenceTree.getTaxonIndex(taxonId));
        } else {
            for (int i = 0; i < queryTree.getChildCount(node); i++) {
                NodeRef child = queryTree.getChild(node, i);
                clade.or(getClades(queryTree, child, referenceTree,clades));
            }
            clades.add(clade);
        }
        return clade;
    }

    private void setupClades(Set<BitSet> targetClades, Tree tree){
        Map<NodeRef,BitSet> nodeCladeMap = getBitSetNodeMap(tree);

        setupClades(targetClades,tree,tree.getRoot(),nodeCladeMap,null);
    }
    private void setupClades(Set<BitSet> targetClades,Tree tree, NodeRef node, Map<NodeRef,BitSet> refBitSetMap,Clade parentClade) {
        Clade clade;
        BitSet cladeSet = refBitSetMap.get(node);

        if(targetClades.contains(cladeSet)) {
            clade = new Clade(cladeSet,cladeList.size(),parentClade);
            if(node!=tree.getRoot()){
                parentClade.addChild(clade);
            }
            cladeRoots[cladeList.size()]=node.getNumber();
            cladeList.add(clade);
            targetClades.remove(cladeSet);
        }else {
            clade = parentClade;
        }
        clade.addNode(node);
        nodeCladeMap[node.getNumber()]=clade;

        for (int i = 0; i < tree.getChildCount(node); i++) {
            NodeRef child = tree.getChild(node,i);
            setupClades(targetClades,tree,child,refBitSetMap,clade);
        }

    }

    /**
     * Gets a HashMap of clade bitsets to nodes in tree. This is useful for comparing the topology of trees
     * @param tree the tree for which clades are being defined
     * @return A HashMap with a BitSet of descendent taxa as the key and a node as value
     */
    private HashMap<NodeRef,BitSet> getBitSetNodeMap(Tree  tree) {
        HashMap<NodeRef,BitSet> map = new HashMap<>();
        getBitSetNodeMap(tree,tree.getRoot(),map);
        return map;
    }

    /**
     *  A private recursive function used by getBitSetNodeMap
     *  This is modeled after the addClades in CladeSet and getClades in compatibility statistic
     * @param tree the tree for which clades are being defined
     * @param node current node
     * @param map map that is being appended to
     */
    private BitSet getBitSetNodeMap( Tree tree, NodeRef node, HashMap<NodeRef,BitSet> map) {
        BitSet bits = new BitSet();
        if (tree.isExternal(node)) {
            String taxonId = tree.getNodeTaxon(node).getId();
            bits.set(tree.getTaxonIndex(taxonId));


        } else {
            for (int i = 0; i < tree.getChildCount(node); i++) {
                NodeRef node1 = tree.getChild(node, i);
                bits.or(getBitSetNodeMap(tree, node1, map));
            }
        }
        map.put(node,bits);
        return bits;
    }




    // *********************************************************************************************
    // *****                             Model implementation                                *******
    // *********************************************************************************************



    @Override
    protected void handleModelChangedEvent(Model model, Object object, int i) {
        fireModelChanged();
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int i, Variable.ChangeType changeType) {

    }

    @Override
    protected void storeState() {
        System.arraycopy(cladeRoots,0,storedCladeRoots,0,cladeRoots.length);
    }

    @Override
    protected void restoreState() {
        int[] tmp = storedCladeRoots;
        storedCladeRoots = cladeRoots;
        cladeRoots= tmp;
    }

    @Override
    protected void acceptState() {

    }

    // Private Stuff

    private final List<Clade> cladeList = new ArrayList<>();
    private final Clade[] nodeCladeMap;
    private final TreeModel treeModel;
    private  int[] cladeRoots;
    private  int[] storedCladeRoots;

    /**
     * A private helper clade class.
     */
    class Clade  implements CladeRef {
        //TODO figure out what in interface and what just here.
        public final static String CLADE = "clade";

        public Clade(BitSet bits, int number, Clade parent){
            this.bits=bits;
            this.number = number;
            this.parent = parent;
        }

        public void addChild(Clade clade) {
            this.children.add(clade);
        }

        public void addNode(NodeRef nodeRef) {
            this.nodes.add(nodeRef.getNumber());
        }

        public int getNumber() {
            return number;
        }

        public List<Integer> getNodes() {
            List<Integer> allNodes = new ArrayList<>(nodes);
            for (Clade child : children) {
                allNodes.add(cladeRoots[child.getNumber()]);
            }
            return allNodes;
        }

        public int getNode(int i) {
            if(i<nodes.size()){
                return nodes.get(i);
            }else if (i<nodes.size()+children.size()){
                return cladeRoots[children.get(i-nodes.size()).getNumber()];
            }else{
                throw new IllegalArgumentException("index " + i +"out of bounds. nodes in the clade " + getSize());
            }
        }

        public int getSize(){
            return nodes.size() + children.size();
        }
        public int getChildCount(){
            return children.size();
        }

        public Clade getChild(int i) {
            return children.get(i);
        }

        public Clade getParent() {
            return parent;
        }

        public BitSet getBits() {
            return bits;
        }

        private final List<Integer> nodes = new ArrayList<>();
        public final List<Clade> children = new ArrayList<>();
        private final Clade parent;
        private final int number;
        private final BitSet bits;
    }
}