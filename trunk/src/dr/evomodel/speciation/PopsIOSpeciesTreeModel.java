package dr.evomodel.speciation;

import java.util.*;

import dr.evolution.tree.*;
import dr.evomodel.tree.TreeLogger;
import dr.evomodelxml.speciation.PopsIOSpeciesTreeModelParser;
import dr.inference.loggers.LogColumn;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Variable;
import dr.inference.operators.OperatorFailedException;
import dr.inference.operators.Scalable;
import dr.util.AlloppMisc;
import dr.util.Author;
import dr.util.Citation;
import jebl.util.FixedBitSet;
import dr.evolution.util.Taxon;
import dr.math.MathUtils;

/**
 * User: Graham
 * Date: 10/05/12
 */

/*
nodes[] and root are the fundamental data here. They are used for
the SlidableTree implementation, not just supporting the MCMC move,
but also the calculations dealing with coalescences (compatibility
and likelihood).

oldtree stores the state so it can be restored after a rejected move. It
only stores taxa at tips, topology and node times. The other data (unions, lineages,
coalescences) are reconstructed from that.

stree is used to implement Tree, which is used for logging and for
calculation of the tree prior (eg Yule model). stree is a copy of
topology and node times  from nodes[].

 */
public class PopsIOSpeciesTreeModel extends AbstractModel implements SlidableTree, Tree, Scalable, TreeLogger.LogUpon {
    private PopsIOSpeciesBindings piosb;
    private PriorComponent[] priorComponents;
    private PopsIONode[] nodes;
    private PopsIONode root;
    private SimpleTree oldtree;
    private SimpleTree stree;

    @Override
    public boolean logNow(long state) {
        // for debugging, set logEvery=0 in XML
        if (state == 40) {
            System.out.println("DEBUGGING: PopsIOSpeciesTreeModel.logNow(), state == 40");
        }
        return state % 10000 == 0;
    }


    public static class PriorComponent {
        private double weight;
        private double alpha;
        private double beta;

        // inv gamma pdf is parameterized as  b^a/Gamma(a)  x^(-a-1)  exp(-b/x)
        // mean is b/(a-1) if a>1, var is  b^2/((a-1)^2 (a-2)) if a>2.

        public PriorComponent(double weight, double alpha, double beta) {
            this.weight = weight;
            this.alpha = alpha;
            this.beta = beta;
        }
    }


    /*
     *  The parent, child, height, taxon fileds implement the basic binary tree.
     *  The other fields are 'working' fields.
     *
     *  nodeNumber is required to implement NodeRef, needed for replaceSlidableRoot()
     *  (possibly one can manage this another way?)
     *
     *  union, coalheights, nlineages are used for calculations: determining compatibility
     *  with gene trees and calculating likelihood.
     */
    private class PopsIONode extends AlloppNode.Abstract implements AlloppNode, NodeRef {
        private PopsIONode parent;
        private PopsIONode child[];
        private double height;
        private Taxon taxon;
        private FixedBitSet union;
        private ArrayList<Double> coalheights;
        private int nlineages;
        private int nodeNumber;

        // dud constuctor
        PopsIONode(int nn) {
            parent = null;
            child = new PopsIONode[0];
            height = -1.0;
            coalheights = new ArrayList<Double>();
            taxon = new Taxon("");
            union = null;
            nodeNumber = nn;
        }


        public String asText(int indentlen) {
            StringBuilder s = new StringBuilder();
            Formatter formatter = new Formatter(s, Locale.US);
            if (child.length == 0) {
                formatter.format("%s ", taxon.getId());
            } else {
                formatter.format("%s ", "+");
            }
            while (s.length() < 20-indentlen) {
                formatter.format("%s", " ");
            }
            formatter.format("%s ", AlloppMisc.nonnegIn8Chars(height));
            formatter.format("%20s ", AlloppMisc.FixedBitSetasText(union));
            formatter.format("%3d  ", nlineages);
            for (int c = 0; c < coalheights.size(); c++) {
                formatter.format(AlloppMisc.nonnegIn8Chars(coalheights.get(c)) + ",");
            }
            return s.toString();
        }


        @Override
        public int getNumber() {
            return nodeNumber;
        }

        @Override
        public void setNumber(int n) {
            nodeNumber = n;
        }

        @Override
        public int nofChildren() {
            return child.length;
        }

        @Override
        public AlloppNode getChild(int ch) {
            return child[ch];
        }

        @Override
        public AlloppNode getAnc() {
            return parent;
        }

        @Override
        public Taxon getTaxon() {
            return taxon;
        }

        @Override
        public double getHeight() {
            return height;
        }

        @Override
        public FixedBitSet getUnion() {
            return union;
        }

        @Override
        public void setChild(int ch, AlloppNode newchild) {
            if (child == null) {
                child = new PopsIONode[2];
            }
            child[ch] = (PopsIONode)newchild;
            newchild.setAnc(this);
        }

        @Override
        public void setAnc(AlloppNode anc) {
            parent = (PopsIONode)anc;

        }

        @Override
        public void setTaxon(String name) {
            this.taxon = new Taxon(name);

        }

        @Override
        public void setHeight(double height) {
            this.height = height;

        }

        @Override
        public void setUnion(FixedBitSet union) {
            this.union = union;

        }

        @Override
        public void addChildren(AlloppNode c0, AlloppNode c1) {
            child = new PopsIONode[2];
            child[0] = (PopsIONode)c0;
            child[0].setAnc(this);
            child[1] = (PopsIONode)c1;
            child[1].setAnc(this);
        }

    }


    public PopsIOSpeciesTreeModel(PopsIOSpeciesBindings piosb, PriorComponent[] priorComponents) {
        super(PopsIOSpeciesTreeModelParser.PIO_SPECIES_TREE);
        this.piosb = piosb;
        this.priorComponents = priorComponents;
        PopsIOSpeciesBindings.SpInfo[] species = piosb.getSpecies();
        int nTaxa = species.length;
        int nNodes = 2 * nTaxa - 1;
        nodes = new PopsIONode[nNodes];
        for (int n = 0; n < nNodes; n++) {
            nodes[n] = new PopsIONode(n);
        }
        ArrayList<Integer> tojoin = new ArrayList<Integer>(nTaxa);
        for (int n = 0; n < nTaxa; n++) {
            nodes[n].setTaxon(species[n].name);
            nodes[n].setHeight(0.0);
            nodes[n].setUnion(piosb.tipUnionFromTaxon(nodes[n].getTaxon()));
            tojoin.add(n);
        }
        double rate = 1.0;
        double treeheight = 0.0;
        for (int i = 0; i < nTaxa-1; i++) {
            int numtojoin = tojoin.size();
            int j = MathUtils.nextInt(numtojoin);
            Integer child0 = tojoin.get(j);
            tojoin.remove(j);
            int k = MathUtils.nextInt(numtojoin-1);
            Integer child1 = tojoin.get(k);
            tojoin.remove(k);
            nodes[nTaxa+i].addChildren(nodes[child0], nodes[child1]);
            nodes[nTaxa+i].setHeight(treeheight + randomnodeheight(numtojoin*rate));
            treeheight = nodes[nTaxa+i].getHeight();
            tojoin.add(nTaxa+i);
        }
        root = nodes[nodes.length - 1];

        double scale = 0.99 * piosb.initialMinGeneNodeHeight() / root.height;
        scaleAllHeights(scale);
        root.fillinUnionsInSubtree(piosb.getSpecies().length);

        stree = makeSimpleTree();
    }




    public List<Citation> getCitations() {
        List<Citation> citations = new ArrayList<Citation>();
        citations.add(new Citation(
                new Author[]{
                        new Author("Graham", "Jones")
                },
                "WORKING TITLE: A multi-species coalescent model with population parameters integrated out",
                "??",  // journal
                Citation.Status.IN_PREPARATION
        ));
        return citations;
    }




    public String toString() {
        int ngt = piosb.numberOfGeneTrees();
        String nl = System.getProperty("line.separator");
        String s = nl + pioTreeAsText() + nl;
        for (int g = 0; g < ngt; g++) {
            s += "Gene tree " + g + nl;
            s += piosb.genetreeAsText(g) + nl;
        }
        s += nl;
        return s;
    }



    public LogColumn[] getColumns() {
        LogColumn[] columns = new LogColumn[1];
        columns[0] = new LogColumn.Default("    species-tree and gene trees", this);
        return columns;
    }



    private int scaleAllHeights(double scale) {
        for (int nn = 0; nn < nodes.length; nn++) {
            nodes[nn].height *= scale;
        }
        return nodes.length;
    }


    public String pioTreeAsText() {
        String header = "topology             height         union         nlin coalheights" + System.getProperty("line.separator");

        String s = "";
        Stack<Integer> x = new Stack<Integer>();
        return header + subtreeAsText(root, s, x, 0, "");
    }

    /*
      * Called from PopsIOSpeciesBindings to check if a node in a gene tree
      * is compatible with the network.
      */
    public boolean coalescenceIsCompatible(double height, FixedBitSet union) {
        PopsIONode node = (PopsIONode) root.nodeOfUnionInSubtree(union);
        return (node.height <= height);
    }


    /*
    * Called from PopsIOSpeciesBindings to remove coalescent information
    * from branches of mullabtree. Required before call to recordCoalescence
    */
    public void clearCoalescences() {
        clearSubtreeCoalescences(root);
    }


    /*
      * Called from PopsIOSpeciesBindings to add a node from a gene tree
      * to its branch in mullabtree.
      */
    public void recordCoalescence(double height, FixedBitSet union) {
        PopsIONode node = (PopsIONode) root.nodeOfUnionInSubtree(union);
        assert (node.height <= height);
        while (node.parent != null  &&  node.parent.height <= height) {
            node = node.parent;
        }
        node.coalheights.add(height);    }


    public void sortCoalescences() {
        for (PopsIONode node : nodes) {
            Collections.sort(node.coalheights);
        }
    }



    /*
      * Records the number of gene lineages at nodes of mullabtree.
      */
    public void recordLineageCounts() {
        recordSubtreeLineageCounts(root);
    }


    public void fixupAfterNodeSlide() {
        stree = makeSimpleTree();
    }

    /*
      * Calculates the log-likelihood for a single gene tree in the network
      *
      * Requires that clearCoalescences(), recordCoalescence(), recordLineageCounts()
      * called to fill tree in nodes[] with information about gene tree coalescences first.
      */
    /*
     * The formula comes from my note at http://www.indriid.com/goteborg/2011-09-23-simple-pop-model.pdf
     * See branchLL() for more.
     */
    public double geneTreeInSpeciesTreeLogLikelihood() {
        return geneTreeInPopsIOSubtreeLogLikelihood(root);
    }


    private String subtreeAsText(PopsIONode node, String s, Stack<Integer> x, int depth, String b) {
        Integer[] y = x.toArray(new Integer[x.size()]);
        StringBuffer indent = new StringBuffer();
        for (int i = 0; i < depth; i++) {
            indent.append("  ");
        }
        for (int i = 0; i < y.length; i++) {
            indent.replace(2*y[i], 2*y[i]+1, "|");
        }
        if (b.length() > 0) {
            indent.replace(indent.length()-b.length(), indent.length(), b);
        }
        s += indent;
        s += node.asText(indent.length());
        s += System.getProperty("line.separator");
        String subs = "";
        if (node.child.length > 0) {
            x.push(depth);
            subs += subtreeAsText(node.child[0], "", x, depth+1, "-");
            x.pop();
            subs += subtreeAsText(node.child[1], "", x, depth+1, "`-");
        }
        return s + subs;
    }


    private double geneTreeInPopsIOSubtreeLogLikelihood(PopsIONode node) {
        double loglike = 0.0;
        for (int i = 0; i < node.child.length; i++) {
            loglike += geneTreeInPopsIOSubtreeLogLikelihood(node.child[i]);
        }
        loglike += branchLLInPopsIOtree(node);
        return loglike;
    }


    private double branchLLInPopsIOtree(PopsIONode node) {
        double loglike = 0.0;
        double t[];
        if (node.parent == null) {
            t = new double[node.coalheights.size() + 2];
            t[0] = node.height;
            t[t.length - 1] = piosb.maxGeneTreeHeight();
            for (int i = 0; i < node.coalheights.size(); i++) {
                t[i + 1] = node.coalheights.get(i);
            }
            loglike += branchLL(t, node.nlineages);
        } else {
            t = new double[node.coalheights.size() + 2];
            t[0] = node.height;
            t[t.length - 1] = node.parent.height;
            for (int i = 0; i < node.coalheights.size(); i++) {
                t[i + 1] = node.coalheights.get(i);
            }
            loglike += branchLL(t, node.nlineages);
        }
        return loglike;
    }

    /*
    * For one branch with tipward time t[0], rootward time t[k+1], k-1 coalescent times t[1]...t[k],
    * and n lineages at tipward end, set
    *      x = sum from i=0 to k of
    *           ((n-i) choose 2) * (t[i+1]-t[i])
    * Then sum over j (j is component index) of
    *     weight[j]  *  b[j]^a[i]  *  (b[j] + x)^-(a[j]+k+1)  *  GAMMA(a[j]+k+1)  /  GAMMA(a[j])
    *
    *     G(z+1) = zG(z)
    *     GAMMA(a[j]+k+1) = (a[j]+k)GAMMA(a[j]+k)
    *                 = (a[j]+k)(a[j]+k-1)GAMMA(a[j]+k-1)
    *                 = ...
    *                 = (a[j]+k)(a[j]+k-1)...a[j]GAMMA(a[j])
    *    so GAMMA(a[j]+k+1)  /  GAMMA(a[j]  = (a[j]+k)(a[j]+k-1)...a[j]
    */
    private double branchLL(double t[], int n) {
        double lhood = 0.0;
        double x = 0.0;
        int k = t.length - 2;
        for (int i = 0; i <= k; i++) {
            x += (t[i+1] - t[i]) * 0.5*(n-i)*(n-i-1);
        }
        for (int j = 0; j < priorComponents.length; j++) {
            double w = priorComponents[j].weight;
            double a = priorComponents[j].alpha;
            double b = priorComponents[j].beta;
            double G = 1.0;
            for (int i = 0; i <= k; i++) {
                G *= (a+i);
            }
            lhood += w * Math.pow(a,b) * Math.pow(b+x, -(a+k+1)) * G;
        }
        return Math.log(lhood);
    }


    private SimpleTree makeSimpleTree() {
        SimpleNode[] snodes = new SimpleNode[nodes.length];
        for (int n = 0; n < nodes.length; n++) {
            snodes[n] = new SimpleNode();
            snodes[n].setTaxon(null);
        }
        makesimplesubtree(snodes, 0, root);
        return new SimpleTree(snodes[nodes.length-1]);
    }


    // for makeSimpleTree()
    private int makesimplesubtree(SimpleNode[] snodes, int nextsn, PopsIONode pionode) {
        if (pionode.child.length == 0) {
            Taxon tx = new Taxon(pionode.taxon.getId());
            if (nextsn >= snodes.length) {
                System.out.println("BUG: makesimplesubtree()");
            }
            snodes[nextsn].setTaxon(tx);
        } else {
            nextsn = makesimplesubtree(snodes, nextsn, pionode.child[0]);
            int subtree0 = nextsn-1;
            nextsn = makesimplesubtree(snodes, nextsn, pionode.child[1]);
            int subtree1 = nextsn-1;
            snodes[nextsn].addChild(snodes[subtree0]);
            snodes[nextsn].addChild(snodes[subtree1]);
        }
        snodes[nextsn].setHeight(pionode.height);
        return nextsn+1;
    }



    private void clearSubtreeCoalescences(PopsIONode node) {
        for (int i = 0; i < node.child.length; i++) {
            clearSubtreeCoalescences(node.child[i]);
        }
        node.coalheights.clear();
    }


    private void recordSubtreeLineageCounts(PopsIONode node) {
        if (node.child.length == 0) {
            int spIndex = piosb.speciesId2index(node.getTaxon().getId());
            node.nlineages = piosb.nLineages(spIndex);
        } else {
            node.nlineages = 0;
            for (int i = 0; i < node.child.length; i++) {
                recordSubtreeLineageCounts(node.child[i]);
                node.nlineages += node.child[i].nlineages - node.child[i].coalheights.size();
            }
        }
    }

    private double randomnodeheight(double rate) {
        return MathUtils.nextExponential(rate) + 1e-6/rate;
        // 1e-6/rate to avoid very tiny heights
    }



    /*******************************************************************************/
    // for Scalable.
    /*******************************************************************************/
    @Override
    public int scale(double factor, int nDims) throws OperatorFailedException {
        int n = scaleAllHeights(factor);
        stree = makeSimpleTree();
        return n;
    }


    @Override
    public String getName() {
        return PopsIOSpeciesTreeModelParser.PIO_SPECIES_TREE;

    }
    /*******************************************************************************/
    // for AbstractModel.
    /*******************************************************************************/
    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
        fireModelChanged();
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Variable.ChangeType type) {
    }

    @Override
    protected void storeState() {
        oldtree = makeSimpleTree();
    }

    @Override
    protected void restoreState() {
        SimpleNode root = (SimpleNode) oldtree.getRoot();
        int nextn = 0;
        nodes = new PopsIONode[oldtree.getNodeCount()];
        for (int nn = 0; nn < nodes.length; nn++) {
            nodes[nn] = new PopsIONode(nn);
        }
        nextn = AlloppNode.Abstract.simpletree2piotree(piosb, nodes, nextn, root);
        this.root = nodes[nextn - 1];
        this.root.fillinUnionsInSubtree(piosb.getSpecies().length);
        stree = makeSimpleTree();
    }

    @Override
    protected void acceptState() {
    }



    /*******************************************************************************/
    // for SlidableTree.
    /*******************************************************************************/
    @Override
    public NodeRef getSlidableRoot() {
        if (root.parent != null) {
            System.out.println("BUG: getSlidableRoot");
        }
        assert root.parent == null;
        return root;
    }

    @Override
    public void replaceSlidableRoot(NodeRef newroot) {
        root = (PopsIONode) newroot;
        root.parent = null;
    }

    @Override
    public int getSlidableNodeCount() {
        return nodes.length;
    }

    @Override
    public Taxon getSlidableNodeTaxon(NodeRef node) {
        assert node == nodes[node.getNumber()];
        return ((PopsIONode)node).getTaxon();
    }


    @Override
    public double getSlidableNodeHeight(NodeRef node) {
        assert node == nodes[node.getNumber()];
        return ((PopsIONode)node).getHeight();
    }

    @Override
    public void setSlidableNodeHeight(NodeRef node, double height) {
        assert node == nodes[node.getNumber()];
        ((PopsIONode)node).height = height;
    }

    @Override
    public boolean isExternalSlidable(NodeRef node) {
        assert node == nodes[node.getNumber()];
        return (((PopsIONode)node).child.length == 0);
    }

    @Override
    public NodeRef getSlidableChild(NodeRef node, int j) {
        assert node == nodes[node.getNumber()];
        return ((PopsIONode)node).child[j];
    }

    @Override
    public void replaceSlidableChildren(NodeRef node, NodeRef lft, NodeRef rgt) {
        assert node == nodes[node.getNumber()];
        assert lft == nodes[lft.getNumber()];
        assert rgt == nodes[rgt.getNumber()];
        assert ((PopsIONode)node).child.length == 2;
        ((PopsIONode)node).child[0] = (PopsIONode)lft;
        ((PopsIONode)node).child[1] = (PopsIONode)rgt;
        ((PopsIONode)lft).parent = (PopsIONode)node;
        ((PopsIONode)rgt).parent = (PopsIONode)node;
    }



    /*******************************************************************************/
    // For Tree
    /*******************************************************************************/

    @Override
    public NodeRef getRoot() {
        return stree.getRoot();
    }

    @Override
    public int getNodeCount() {
        return stree.getNodeCount();
    }

    @Override
    public NodeRef getNode(int i) {
        return stree.getNode(i);
    }

    @Override
    public NodeRef getInternalNode(int i) {
        return stree.getInternalNode(i);
    }

    @Override
    public NodeRef getExternalNode(int i) {
        return stree.getExternalNode(i);
    }

    @Override
    public int getExternalNodeCount() {
        return stree.getExternalNodeCount();
    }

    @Override
    public int getInternalNodeCount() {
        return stree.getInternalNodeCount();
    }

    @Override
    public Taxon getNodeTaxon(NodeRef node) {
        return stree.getNodeTaxon(node);
    }

    @Override
    public boolean hasNodeHeights() {
        return stree.hasNodeHeights();
    }

    @Override
    public double getNodeHeight(NodeRef node) {
        return stree.getNodeHeight(node);
    }

    @Override
    public boolean hasBranchLengths() {
        return stree.hasBranchLengths();
    }

    @Override
    public double getBranchLength(NodeRef node) {
        return stree.getBranchLength(node);
    }

    @Override
    public double getNodeRate(NodeRef node) {
        return stree.getNodeRate(node);
    }

    @Override
    public Object getNodeAttribute(NodeRef node, String name) {
        return stree.getNodeAttribute(node, name);
    }

    @Override
    public Iterator getNodeAttributeNames(NodeRef node) {
        return stree.getNodeAttributeNames(node);
    }

    @Override
    public boolean isExternal(NodeRef node) {
        return stree.isExternal(node);
    }

    @Override
    public boolean isRoot(NodeRef node) {
        return stree.isRoot(node);
    }

    @Override
    public int getChildCount(NodeRef node) {
        return stree.getChildCount(node);
    }

    @Override
    public NodeRef getChild(NodeRef node, int j) {
        return stree.getChild(node, j);
    }

    @Override
    public NodeRef getParent(NodeRef node) {
        return stree.getParent(node);
    }

    @Override
    public Tree getCopy() {
        return stree.getCopy();
    }

    @Override
    public void setAttribute(String name, Object value) {
        stree.setAttribute(name, value);
    }

    @Override
    public Object getAttribute(String name) {
        return stree.getAttribute(name);
    }

    @Override
    public Iterator<String> getAttributeNames() {
        return stree.getAttributeNames();
    }

    @Override
    public int getTaxonCount() {
        return stree.getTaxonCount();
    }

    @Override
    public Taxon getTaxon(int taxonIndex) {
        return stree.getTaxon(taxonIndex);
    }

    @Override
    public String getTaxonId(int taxonIndex) {
        return stree.getTaxonId(taxonIndex);
    }

    @Override
    public int getTaxonIndex(String id) {
        return stree.getTaxonIndex(id);
    }

    @Override
    public int getTaxonIndex(Taxon taxon) {
        return stree.getTaxonIndex(taxon);
    }

    @Override
    public List<Taxon> asList() {
        return stree.asList();
    }

    @Override
    public Object getTaxonAttribute(int taxonIndex, String name) {
        return stree.getTaxonAttribute(taxonIndex, name);
    }

    @Override
    public Iterator<Taxon> iterator() {
        return stree.iterator();
    }

    @Override
    public Type getUnits() {
        return stree.getUnits();
    }

    @Override
    public void setUnits(Type units) {
        stree.setUnits(units);
    }


}

