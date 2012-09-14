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
    private PopsIONode[] pionodes;
    private int rootn;
    private PopsIONode[] oldpionodes;
    private int oldrootn;
    private SimpleTree stree;

    @Override
    public boolean logNow(long state) {
        // for debugging, set logEvery=0 in XML:
        //    <!--  species tree log file.  -->
        //    <logTree id="pioTreeFileLog" logEvery="0" fileName="C:\Users\....

        /*if (state == 40) {
            System.out.println("DEBUGGING: PopsIOSpeciesTreeModel.logNow(), state == 40");
        }  */

        if (state <= 100) {
            return true;
        }
        if (state <= 10000) {
            return (state % 100) == 0;
        }
        return (state % 10000) == 0;
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
        private int anc;
        private int lft;
        private int rgt;
        private double height;
        private Taxon taxon;
        private FixedBitSet union;
        private ArrayList<Double> coalheights;
        private int nlineages;
        private int nodeNumber;

        // dud constuctor
        PopsIONode(int nn) {
            anc = -1;
            lft = -1;
            rgt = -1;
            height = -1.0;
            coalheights = new ArrayList<Double>();
            taxon = new Taxon("");
            union = null;
            nodeNumber = nn;
        }



        // copy constructor
        public PopsIONode(PopsIONode node) {
            anc = node.anc;
            lft = node.lft;
            rgt = node.rgt;
            nodeNumber = node.nodeNumber;
            copyNonTopologyFields(node);
        }




        private void copyNonTopologyFields(PopsIONode node) {
            height = node.height;
            taxon = new Taxon(node.taxon.getId());
            nlineages = node.nlineages;
            if (node.union == null) {
                union = null;
            } else {
                union = new FixedBitSet(node.union);
            }
            coalheights = new ArrayList<Double>();
            for (int i = 0; i < node.coalheights.size(); i++) {
                coalheights.add(node.coalheights.get(i));
            }
        }




        public String asText(int indentlen) {
            StringBuilder s = new StringBuilder();
            Formatter formatter = new Formatter(s, Locale.US);
            if (lft < 0) {
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
            return (lft < 0) ? 0 : 2;
        }

        @Override
        public AlloppNode getChild(int ch) {
            return ch==0 ? pionodes[lft] : pionodes[rgt];
        }

        @Override
        public AlloppNode getAnc() {
            return pionodes[anc];
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
            int newch = ((PopsIONode)newchild).nodeNumber;
            if (ch == 0) {
                lft = newch;
            } else {
                rgt = newch;
            }
        }

        @Override
        public void setAnc(AlloppNode anc) {
            this.anc = ((PopsIONode)anc).nodeNumber;

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
            lft = ((PopsIONode)c0).nodeNumber;
            pionodes[lft].anc = nodeNumber;
            rgt = ((PopsIONode)c1).nodeNumber;
            pionodes[rgt].anc = nodeNumber;
        }

    }


    public PopsIOSpeciesTreeModel(PopsIOSpeciesBindings piosb, PriorComponent[] priorComponents) {
        super(PopsIOSpeciesTreeModelParser.PIO_SPECIES_TREE);
        this.piosb = piosb;
        this.priorComponents = priorComponents;
        PopsIOSpeciesBindings.SpInfo[] species = piosb.getSpecies();
        int nTaxa = species.length;
        int nNodes = 2 * nTaxa - 1;
        pionodes = new PopsIONode[nNodes];
        for (int n = 0; n < nNodes; n++) {
            pionodes[n] = new PopsIONode(n);
        }
        ArrayList<Integer> tojoin = new ArrayList<Integer>(nTaxa);
        for (int n = 0; n < nTaxa; n++) {
            pionodes[n].setTaxon(species[n].name);
            pionodes[n].setHeight(0.0);
            pionodes[n].setUnion(piosb.tipUnionFromTaxon(pionodes[n].getTaxon()));
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
            pionodes[nTaxa+i].addChildren(pionodes[child0],pionodes[child1]);
            pionodes[nTaxa+i].setHeight(treeheight + randomnodeheight(numtojoin*rate));
            treeheight = pionodes[nTaxa+i].getHeight();
            tojoin.add(nTaxa+i);
        }
        rootn = pionodes.length - 1;

        double scale = 0.99 * piosb.initialMinGeneNodeHeight() / pionodes[rootn].height;
        scaleAllHeights(scale);
        pionodes[rootn].fillinUnionsInSubtree(piosb.getSpecies().length);

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
        for (int nn = 0; nn < pionodes.length; nn++) {
            pionodes[nn].height *= scale;
        }
        return pionodes.length;
    }


    public String pioTreeAsText() {
        String header = "topology             height         union         nlin coalheights" + System.getProperty("line.separator");

        String s = "";
        Stack<Integer> x = new Stack<Integer>();
        return header + subtreeAsText(pionodes[rootn], s, x, 0, "");
    }

    /*
      * Called from PopsIOSpeciesBindings to check if a node in a gene tree
      * is compatible with the network.
      */
    public boolean coalescenceIsCompatible(double height, FixedBitSet union) {
        PopsIONode node = (PopsIONode) pionodes[rootn].nodeOfUnionInSubtree(union);
        return (node.height <= height);
    }


    /*
    * Called from PopsIOSpeciesBindings to remove coalescent information
    * from branches of mullabtree. Required before call to recordCoalescence
    */
    public void clearCoalescences() {
        clearSubtreeCoalescences(pionodes[rootn]);
    }


    /*
      * Called from PopsIOSpeciesBindings to add a node from a gene tree
      * to its branch in mullabtree.
      */
    public void recordCoalescence(double height, FixedBitSet union) {
        PopsIONode node = (PopsIONode) pionodes[rootn].nodeOfUnionInSubtree(union);
        assert (node.height <= height);
        while (node.anc >= 0  &&  pionodes[node.anc].height <= height) {
            node = pionodes[node.anc];
        }
        node.coalheights.add(height);

    }


    public void sortCoalescences() {
        for (PopsIONode node : pionodes) {
            Collections.sort(node.coalheights);
        }
    }



    /*
      * Records the number of gene lineages at nodes of mullabtree.
      */
    public void recordLineageCounts() {
        recordSubtreeLineageCounts(pionodes[rootn]);
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
        return geneTreeInPopsIOSubtreeLogLikelihood(pionodes[rootn]);
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
        if (node.lft >= 0) {
            x.push(depth);
            subs += subtreeAsText(pionodes[node.lft], "", x, depth+1, "-");
            x.pop();
            subs += subtreeAsText(pionodes[node.rgt], "", x, depth+1, "`-");
        }
        return s + subs;
    }


    private double geneTreeInPopsIOSubtreeLogLikelihood(PopsIONode node) {
        double loglike = 0.0;
        if (node.lft >= 0) {
            loglike += geneTreeInPopsIOSubtreeLogLikelihood(pionodes[node.lft]);
            loglike += geneTreeInPopsIOSubtreeLogLikelihood(pionodes[node.rgt]);
        }
        loglike += branchLLInPopsIOtree(node);
        return loglike;
    }


    private double branchLLInPopsIOtree(PopsIONode node) {
        double loglike = 0.0;
        double t[];
        if (node.anc < 0) {
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
            t[t.length - 1] = pionodes[node.anc].height;
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
        SimpleNode[] snodes = new SimpleNode[pionodes.length];
        for (int n = 0; n < pionodes.length; n++) {
            snodes[n] = new SimpleNode();
            snodes[n].setTaxon(null);
        }
        makesimplesubtree(snodes, 0, pionodes[rootn]);
        return new SimpleTree(snodes[pionodes.length-1]);
    }


    // for makeSimpleTree()
    private int makesimplesubtree(SimpleNode[] snodes, int nextsn, PopsIONode pionode) {
        if (pionode.lft < 0) {
            Taxon tx = new Taxon(pionode.taxon.getId());
            if (nextsn >= snodes.length) {
                System.out.println("BUG: makesimplesubtree()");
            }
            snodes[nextsn].setTaxon(tx);
        } else {
            nextsn = makesimplesubtree(snodes, nextsn, pionodes[pionode.lft]);
            int subtree0 = nextsn-1;
            nextsn = makesimplesubtree(snodes, nextsn, pionodes[pionode.rgt]);
            int subtree1 = nextsn-1;
            snodes[nextsn].addChild(snodes[subtree0]);
            snodes[nextsn].addChild(snodes[subtree1]);
        }
        snodes[nextsn].setHeight(pionode.height);
        return nextsn+1;
    }



    private void clearSubtreeCoalescences(PopsIONode node) {
        if (node.lft >= 0) {
            clearSubtreeCoalescences(pionodes[node.lft]);
            clearSubtreeCoalescences(pionodes[node.rgt]);
        }
        if (node == null) {
            System.out.println("BUG");
        }
        if (node.coalheights == null) {
            System.out.println("BUG");
        }
        node.coalheights.clear();
    }


    private void recordSubtreeLineageCounts(PopsIONode node) {
        if (node.lft < 0) {
            int spIndex = piosb.speciesId2index(node.getTaxon().getId());
            node.nlineages = piosb.nLineages(spIndex);
        } else {
            node.nlineages = 0;
            recordSubtreeLineageCounts(pionodes[node.lft]);
            node.nlineages += pionodes[node.lft].nlineages - pionodes[node.lft].coalheights.size();
            recordSubtreeLineageCounts(pionodes[node.rgt]);
            node.nlineages += pionodes[node.rgt].nlineages - pionodes[node.rgt].coalheights.size();
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
        oldpionodes = new PopsIONode[pionodes.length];
        for (int n = 0; n < oldpionodes.length; n++) {
            oldpionodes[n] = new PopsIONode(pionodes[n]);
        }
        oldrootn = rootn;
    }

    @Override
    protected void restoreState() {
        pionodes = new PopsIONode[oldpionodes.length];
        for (int n = 0; n < pionodes.length; n++) {
            pionodes[n] = new PopsIONode(oldpionodes[n]);
        }
        rootn = oldrootn;
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
        assert pionodes[rootn].anc < 0;
        return pionodes[rootn];

    }

    @Override
    public void replaceSlidableRoot(NodeRef newroot) {
        rootn = newroot.getNumber();
        pionodes[rootn].anc = -1;
    }

    @Override
    public int getSlidableNodeCount() {
        return pionodes.length;
    }

    @Override
    public Taxon getSlidableNodeTaxon(NodeRef node) {
        assert node == pionodes[node.getNumber()];
        return ((PopsIONode)node).getTaxon();
    }


    @Override
    public double getSlidableNodeHeight(NodeRef node) {
        assert node == pionodes[node.getNumber()];
        return ((PopsIONode)node).getHeight();
    }

    @Override
    public void setSlidableNodeHeight(NodeRef node, double height) {
        assert node == pionodes[node.getNumber()];
        ((PopsIONode)node).height = height;
    }

    @Override
    public boolean isExternalSlidable(NodeRef node) {
        return (pionodes[node.getNumber()].lft < 0);
    }

    @Override
    public NodeRef getSlidableChild(NodeRef node, int j) {
        int n = node.getNumber();
        return j == 0 ? pionodes[ pionodes[n].lft ] : pionodes[ pionodes[n].rgt ];
    }

    @Override
    public void replaceSlidableChildren(NodeRef node, NodeRef lft, NodeRef rgt) {
        int nn = node.getNumber();
        int lftn = lft.getNumber();
        int rgtn = rgt.getNumber();
        assert pionodes[nn].lft >= 0;
        pionodes[nn].lft = lftn;
        pionodes[nn].rgt = rgtn;
        pionodes[lftn].anc = pionodes[nn].nodeNumber;
        pionodes[rgtn].anc = pionodes[nn].nodeNumber;
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

