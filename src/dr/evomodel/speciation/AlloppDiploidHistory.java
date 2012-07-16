package dr.evomodel.speciation;

import java.util.ArrayList;
import java.util.List;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.SimpleNode;
import dr.evolution.tree.SimpleTree;
import dr.evolution.tree.SlidableTree;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxon;
import dr.util.AlloppMisc;
import jebl.util.FixedBitSet;

/**
 * AlloppDiploidHistory represents part of the network before hybridizations.
 * It is basically a tree with some tips representing diploid species at present time
 * and others representing the points at which hybridization occurs.
 * 
 * @author Graham Jones
 *         Date: 13/03/2012
 */


/*
 * This a tree with some tips representing diploid species at present time
 * and others representing the points at which hybridization (to form a
 * tetraploid) occurs. The latter are in pairs and have times before present,
 * and I call them `hydridization tips' or `hyb-tips'.
 * 
 * The purpose of this class is to extract the part of the network before
 * hybridizations (=`diploid history') in a form (the SimpleTree simptree) 
 * that can be subjected to Mau-type moves. 
 * Then, after a move, this class provides ways of re-making the 
 * network, as a set of HybHistory's, one per tetraploid subtree, and a SimpleTree
 * which is the diploid tree on extant tips only, with no hyb-tips.
 * 
 * The sequence goes
 * 
 * 1. Set of AlloppLeggedTree's representing network (AlloppSpeciesNetworkModel)
 * 2. The diploid history as a tree in DipHistNode[] nodes
 * 3. The SimpleTree TESTtree for testing only
 * 4. Move
 * 5. Extract the HybHistory's
 * 6. Obtain the diploid tree on extant tips only
 * 
 * 1 to 3 are done by constructor
 * 5 is extractHybHistory(), 6 is ditreeFromDipHist()
 * 
 * 4 is done elsewhere. 
 * 5 and 6 are used by AlloppSpeciesNetworkModel to replace diploid history.
 */



public class AlloppDiploidHistory implements SlidableTree {
	private DipHistNode[] nodes;
	private int rootn;
	private int nextn;
	private AlloppSpeciesBindings apsp;
	
	// for replacing diploid history. Collects
	public class HybHistory {
		public double hybheight;
		public double splitheight;
		public AlloppTreeLeg legs[];

		public HybHistory(AlloppTreeLeg[] legs, double hybheight, double splitheight) {
			this.hybheight = hybheight;
			this.splitheight = splitheight;
			this.legs = legs;
		}
			
	}
	
	
	/*
	 * parent, child[] implement the tree topology.
	 * 
	 * height is the node height; can be > 0 for tips.
	 * 
	 * union is a spseq-union. At diploid tips they are used normally. 
	 * Hyb-tips have unions derived from the tetratree and leg index (0 or 1). 
	 * Unions are taken towards the root.
	 * They therefore are `the same thing' as footUnion's in AlloppTreeLeg's.
	 * 
	 * For a hyb-tip, tettree specifies the index of the tetratree whose
	 * root comes from it. It is not used for other nodes.   
	 */
	private class DipHistNode extends AlloppNode.Abstract implements AlloppNode, NodeRef {
		private DipHistNode parent;
		private DipHistNode child[];
		private double height;
		private Taxon taxon;
		private FixedBitSet union;
		private int tettree;
		private int leg;
		private int nodeNumber;
		private boolean del;  // for marking nodes that don't go into ditree on extant diploids
		
		// dud constuctor
		DipHistNode(int nn) {
			parent = null;
			child = new DipHistNode[0];
			height = -1.0;
			taxon = new Taxon("");
			union = null;
			tettree = -1;
			leg = -1;
			nodeNumber = nn;
			del = false;
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
		public double getHeight() {
			return height;
		}
		
		@Override
		public Taxon getTaxon() {
			return taxon;
		}
		
		@Override
		public FixedBitSet getUnion() {
			return union;
		}

		@Override
		public void setChild(int ch, AlloppNode newchild) {
			if (child == null) {
				child = new DipHistNode[2];
			}
			child[ch] = (DipHistNode)newchild;
		}

		@Override
		public void setAnc(AlloppNode anc) {
			parent = (DipHistNode)anc;
			
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
			child = new DipHistNode[2];
			child[0] = (DipHistNode)c0;
			child[0].setAnc(this);
			child[1] = (DipHistNode)c1;
			child[1].setAnc(this);
			}
			
		

		@Override
		public int nofChildren() {
			return child.length;
		}

		@Override
		public int getNumber() {
			return nodeNumber;
		}

		@Override
		public void setNumber(int n) {
			nodeNumber = n;
			
		}



	}
	
	
	/*
	 * this is for testing, debugging. It is filled from
	 * the SimpleTree simptree, not directly from DipHistNode's
	 */
	public class DipHistTESTINGNode {
		public int num;
		public int ch0;
		public int ch1;
		public double height;
		public String name;
		public FixedBitSet union;
		public int tettree;
		public int leg;
		
		DipHistTESTINGNode(int num, int ch0, int ch1, double height, String name, FixedBitSet union, int tettree, int leg) {
			this.num = num;
			this.ch0 = ch0;
			this.ch1 = ch1;
			this.height = height;
			this.name = name;
			this.union = union;
			this.tettree = tettree;
			this.leg = leg;
		}
		
	}
	
	


	// for extracting the HybHistory's
	private class HybTipPair {
		public int n0;
		public int n1;
		HybTipPair(int n0, int n1) {
			this.n0 = n0;
			this.n1 = n1;
		}
	}

	
	
	/*
	 * Constructor makes AlloppDiploidHistory from collection of trees in network. 
	 * apsp needed for speciesseqEmptyUnion(), speciesseqToTipUnion()
	 * 
	 */
	AlloppDiploidHistory(AlloppLeggedTree[][] trees, AlloppSpeciesBindings apsp) {
		this.apsp = apsp;
		// Count tips for the tree to be made
		int ntips = 0;
		AlloppLeggedTree[] ditrees = trees[AlloppSpeciesNetworkModel.DITREES];
		AlloppLeggedTree[]  tetratrees = trees[AlloppSpeciesNetworkModel.TETRATREES];
		if (ditrees.length > 0) {
			ntips += ditrees[0].getTaxonCount(); // never more than one diploid tree
		}
		for (int i = 0; i < tetratrees.length; i++) {
			ntips += 2;
		}
		// Make array of dud nodes
		nodes = new DipHistNode[2 * ntips - 1];
		for (int i = 0; i < nodes.length; i++) {
			nodes[i] = new DipHistNode(i);
		}
		// Copy homoploid trees into array
		nextn = 0;
		// Ditree if it exists
        // grjtodo I don't think diroot==null case is sustainable.
		DipHistNode diroot = null;
        if (ditrees.length > 0) {
            SimpleNode root = (SimpleNode) ditrees[0].getRoot();
            nextn = AlloppNode.Abstract.simpletree2allopptree(apsp, nodes, nextn, root, false, 0);
            rootn = nextn - 1;
            diroot = nodes[rootn];
        }
		if (diroot != null &&  diroot.height <= 0.0) {
			System.err.println("AlloppDiploidHistory constructor: bug");
		}
        assert diroot != null;
        diroot.fillinUnionsInSubtree(apsp.numberOfSpSeqs());

		assert (diroot != null  ||  tetratrees.length == 1);

		List<AlloppLegLink> leglinks = new ArrayList<AlloppLegLink>();
		for (int i = 0; i < tetratrees.length; i++) {
			// Tetratrees. Two hyb-tips each.
			addhybtiptodiphisttree(apsp, tetratrees, i, 0, "tt" + i + "leg0");
			int r0 = nextn - 1;
			addhybtiptodiphisttree(apsp, tetratrees, i, 1, "tt" + i + "leg1");
			int r1 = nextn - 1;
			nextn = AlloppNode.Abstract.collectLegLinksAndMakeLegJoin(leglinks, tetratrees[i], diroot, nodes, nextn, r0, r1);
		}
		if (diroot == null) {
			rootn = nextn - 1;
		}
		AlloppNode.Abstract.convertLegLinks(nodes, nextn, leglinks);
		assert diphistOK();
		nodes[rootn].fillinUnionsInSubtree(apsp.numberOfSpSeqs());
		makesimpletree();
	}
	
	
	
	public AlloppDiploidHistory(AlloppSpeciesNetworkModel asnm) {
		this(asnm.getAllHomoploidTrees(), asnm.getSpeciesBindings());
	}
	
	

	
	
	
	/*
	 * For testing
	 */
	public String diphistTreeAsUniqueNewick() {
		SimpleTree TESTtree = makesimpletree();
		return Tree.Utils.uniqueNewick(TESTtree, TESTtree.getRoot());
	}
	
	/*
	 * For testing
	 */	
	public String diphistTreeAsText() {
		SimpleTree TESTtree = makesimpletree();
		String s = AlloppMisc.SimpleTreeAsTextualNodeList(TESTtree);
		return s;
		//return simptree.toString();
	} 
	
	
	/*
	 * For testing
	 */	
	public DipHistTESTINGNode[] diphistTreeAsNodeList() {
		SimpleTree TESTtree = makesimpletree();
		int nnodes = TESTtree.getNodeCount();
		DipHistTESTINGNode[] tnodes = new DipHistTESTINGNode[nnodes];
		for (int n = 0;  n < nnodes;  ++n) {
			NodeRef node = TESTtree.getNode(n);
			int num = node.getNumber();
			
			int ch0 = -1; 
			int ch1 = -1;
			double height = -1.0;
			String name = "*";
			FixedBitSet union = apsp.speciesseqEmptyUnion();
			int tettree = -1;
			int leg = -1;
			
			if (TESTtree.getChildCount(node) > 0) {
				assert(TESTtree.getChildCount(node) == 2);
				ch0 = TESTtree.getChild(node, 0).getNumber();
				ch1 = TESTtree.getChild(node, 1).getNumber();
			}
			height = TESTtree.getNodeHeight(node);

			Taxon tx = TESTtree.getNodeTaxon(node);
			if (tx != null) {
				String txid = tx.getId();
				if (txid != null && txid.length() > 0) {
					name = txid;
				}
			}
			
			Object x = TESTtree.getNodeAttribute(node, "union");
			if (x != null) {
				union = (FixedBitSet)x;
			}
			x = TESTtree.getNodeAttribute(node, "tettree");
			if (x != null) {
				tettree = (Integer)x;
			}			
			x = TESTtree.getNodeAttribute(node, "leg");
			if (x != null) {
				leg = (Integer)x;
			}
			DipHistTESTINGNode tnode = new DipHistTESTINGNode(num, ch0, ch1, height, name, union, tettree, leg);	
			tnodes[n] = tnode;
		}
		return tnodes;
	}
	
	
	/*
	 * For testing
	 */
	public boolean diphistOK() {
		boolean bad = false;
		for (int i = 0; i < nodes.length; i++) {
			if (nodes[i].child.length > 0) {
				bad = bad || (nodes[i].child[0].parent != nodes[i]);
				bad = bad || (nodes[i].child[1].parent != nodes[i]);
			}
		}
		bad = bad || (nodes[rootn].parent != null);
		return !bad;
	}
	
	
	
	/*
	 * For testing only (2012-05-04).
	 */
	public int nofTettrees() {
		int ttreeattribcount = 0;
		for (int n = 0; n < nodes.length; ++n) {
			if (nodes[n].tettree >= 0) {
				ttreeattribcount ++;
			}
		}
		assert(ttreeattribcount > 0);
		assert(ttreeattribcount % 2 == 0);
		return ttreeattribcount/2;
	}
	
	
	public boolean nodeIsDiploidTip(NodeRef node) {
		return (nodes[node.getNumber()].tettree < 0);
	}

	
	public FixedBitSet getNodeUnion(NodeRef node) {
		return nodes[node.getNumber()].getUnion();
	}
	
	 
	/*
	 *  Extracts the HybHistory for one tetratree, Marks nodes at and after
	 *  hybridization for deletion.
	 *  
	 */
	public HybHistory extractHybHistory(int tt) {
	
		// recalc unions with empty unions at hyb-tips and possibly new topology after move
		for (int i = 0;  i < nodes.length; i++) {
			if (nodes[i].child.length == 0  &&  nodes[i].tettree < 0) {
				
			} else {
				nodes[i].setUnion(apsp.speciesseqEmptyUnion());
			}
		}
		nodes[rootn].fillinUnionsInSubtree(apsp.numberOfSpSeqs());
		
		// get the HybHistory and mark nodes for deletion
		// joined legs are a complication
		HybTipPair htp = getHybTipPair(tt);
		DipHistNode node0 = nodes[htp.n0];
		DipHistNode node1 = nodes[htp.n1];
		DipHistNode anc0 = node0.parent;
		DipHistNode anc1 = node1.parent;
		node0.del = true;
		node1.del = true;
		anc0.del = true;
		anc1.del = true;
		
		AlloppTreeLeg legs[];
		double splitheight = -1.0;
		if (anc0 == anc1) {
			legs = new AlloppTreeLeg[1];
			assert(anc0.parent != null);
			DipHistNode ancanc = anc0.parent;
			assert(ancanc.union != null);
			ancanc.del = true;
			legs[0] = new AlloppTreeLeg(ancanc.height);
			legs[0].footUnion = ancanc.union;
			splitheight = anc0.height;
		} else {
			legs = new AlloppTreeLeg[2];
			legs[0] = new AlloppTreeLeg(anc0.height);
			legs[1] = new AlloppTreeLeg(anc1.height);
			legs[0].footUnion = anc0.union;
			legs[1].footUnion = anc1.union;
		}
		double hybhgt = node0.height;
		assert(hybhgt == node1.height);		
		HybHistory hybtip = new HybHistory(legs, hybhgt, splitheight);
		return hybtip;
	}
		
	
	/*
	 * Gets the diploid tree on extant tips only
	 */
	public SimpleTree ditreeFromDipHist() {
		int ndnodes = 0;
		for (int n = 0; n < nodes.length; n++) {
			if (!nodes[n].del) {
				ndnodes++;
			}
		}
		SimpleNode[] dnodes = new SimpleNode[ndnodes];
		ditreeSubtree(dnodes, 0, nodes[rootn]);
		return new SimpleTree(dnodes[ndnodes-1]);
	}
	
	
	
	
	/*
	 * **************************************************************************
	 *                      PRIVATE methods
	 * **************************************************************************
	 */

	
	
	// For main constructor
	private void addhybtiptodiphisttree(AlloppSpeciesBindings apsp, AlloppLeggedTree[] tetratrees, int tt, int leg, String name) {
		DipHistNode[] tetnodes;
		int tetrootn;
		int tetnextn;
		double hybheight;
		
		// make a temporary tettree in order to find the union at root. 
		tetnodes = new DipHistNode[tetratrees[tt].getNodeCount()];
		for (int i = 0; i < tetnodes.length; i++) {
			tetnodes[i] = new DipHistNode(i);
		}
		tetnextn = 0;
		hybheight = tetratrees[tt].getHybridHeight();
		SimpleNode root = (SimpleNode)tetratrees[tt].getRoot();
		tetnextn = AlloppNode.Abstract.simpletree2allopptree(apsp, tetnodes, tetnextn, root, true, leg);
		tetrootn = tetnextn - 1;
		tetnodes[tetrootn].fillinUnionsInSubtree(apsp.numberOfSpSeqs());
		
		nodes[nextn].setUnion(tetnodes[tetrootn].getUnion());
		nodes[nextn].setHeight(hybheight);
		nodes[nextn].setTaxon(name);
		nodes[nextn].tettree = tt;
		nodes[nextn].leg = leg;
		nextn++;
	}	

	
	
	
	// for testing 
	private SimpleTree makesimpletree() {
    	SimpleNode[] snodes = new SimpleNode[nodes.length];
    	for (int n = 0; n < nodes.length; n++) {
    		snodes[n] = new SimpleNode();
    		snodes[n].setTaxon(null); // I use taxon==null to identify joined leg node when removing hybtips
    	}
		makesimplesubtree(snodes, 0, nodes[rootn]);
		return new SimpleTree(snodes[nodes.length-1]);
	}
 	
	
	// for makesimpletree()
	private int makesimplesubtree(SimpleNode[] snodes, int nextsn, DipHistNode dhnode) {
		if (dhnode.child.length == 0) {
			Taxon tx = new Taxon(dhnode.taxon.getId());
			snodes[nextsn].setTaxon(tx);
			FixedBitSet union;
			union = dhnode.union;
			snodes[nextsn].setAttribute("union", union);
			if (dhnode.tettree >= 0) {
				snodes[nextsn].setAttribute("tettree", dhnode.tettree);
				snodes[nextsn].setAttribute("leg", dhnode.leg);
			}
		} else {
			nextsn = makesimplesubtree(snodes, nextsn, dhnode.child[0]);
			int subtree0 = nextsn-1;
			nextsn = makesimplesubtree(snodes, nextsn, dhnode.child[1]);
			int subtree1 = nextsn-1;
			snodes[nextsn].addChild(snodes[subtree0]);
			snodes[nextsn].addChild(snodes[subtree1]);
			FixedBitSet union0 = (FixedBitSet)snodes[subtree0].getAttribute("union");
			FixedBitSet union1 = (FixedBitSet)snodes[subtree1].getAttribute("union");
			FixedBitSet union = apsp.speciesseqEmptyUnion();
			union.union(union0);
			union.union(union1);
			snodes[nextsn].setAttribute("union", union);
		}
		snodes[nextsn].setHeight(dhnode.height);
		return nextsn+1;
	}

	
	
	/* ****************** For after MCMC move*****************************************/

	
	// for extractHybHistory()
	private HybTipPair getHybTipPair(int tt) {
		int n0 = -1;
		int n1 = -1;
		for (int n = 0; n < nodes.length; ++n) {
			if (nodes[n].tettree == tt) {
				if (n0 < 0) { n0 = n; } else { n1 = n; }
			}
		}
		assert(n0 >= 0  &&  n1 >= 0);
		assert(nodes[n0].parent != null);
		assert(nodes[n1] != null);
		return new HybTipPair(n0, n1);
	}
	
	
	
	// for ditreeSubtree()
	private DipHistNode findKeeperNode(DipHistNode dhnode) {
		if (!dhnode.del) {
			return dhnode;
		} else {
			if (dhnode.child.length > 0) {
				assert(dhnode.child.length == 2);
				DipHistNode ch0 = findKeeperNode(dhnode.child[0]);
				DipHistNode ch1 = findKeeperNode(dhnode.child[1]);
				if (ch0 != null) {
					return ch0;
				} else {
					return ch1;
				}
			} else {
				return null;
			}
		}
	}

	// for ditreeFromDipHist()
	private int ditreeSubtree(SimpleNode[] dnodes, int nextdn, DipHistNode dhnode) {
		int ch0 = -1;
		int ch1 = -1;
		if (dhnode.child.length == 0) {
		} else {
			nextdn = ditreeSubtree(dnodes, nextdn, findKeeperNode(dhnode.child[0]));
			ch0 = nextdn - 1;
			nextdn = ditreeSubtree(dnodes, nextdn, findKeeperNode(dhnode.child[1]));
			ch1 = nextdn - 1;
		}
		if (nextdn >= dnodes.length) {
			System.out.println("******** " + nextdn);
		}
		dnodes[nextdn] = new SimpleNode();
		dnodes[nextdn].setHeight(dhnode.height);
		dnodes[nextdn].setTaxon(new Taxon(dhnode.taxon.getId()));
		if (ch0 >= 0) {
			dnodes[nextdn].addChild(dnodes[ch0]);
			dnodes[nextdn].addChild(dnodes[ch1]);
		}
		nextdn++;
		return nextdn;
	}



	@Override
	public NodeRef getSlidableRoot() {
		assert nodes[rootn].parent == null;
		return nodes[rootn];
	}



	@Override
	public void replaceSlidableRoot(NodeRef root) {
		rootn = root.getNumber();
		nodes[rootn].parent = null;
	}



	@Override
	public int getSlidableNodeCount() {
		return nodes.length;
	}



	@Override
	public double getSlidableNodeHeight(NodeRef node) {
		return nodes[node.getNumber()].getHeight();
	}

    @Override
    public Taxon getSlidableNodeTaxon(NodeRef node) {
        return nodes[node.getNumber()].getTaxon();
    }

	@Override
	public void setSlidableNodeHeight(NodeRef node, double height) {
		nodes[node.getNumber()].height = height;
		
	}



	@Override
	public boolean isExternalSlidable(NodeRef node) {
		return (nodes[node.getNumber()].child.length == 0);
	}



	@Override
	public NodeRef getSlidableChild(NodeRef node, int j) {
		return nodes[node.getNumber()].child[j];
	}



	@Override
	public void replaceSlidableChildren(NodeRef node, NodeRef lft, NodeRef rgt) {
		int nn = node.getNumber();
		int lftn = lft.getNumber();
		int rgtn = rgt.getNumber();
		assert nodes[nn].child.length == 2;
		nodes[nn].child[0] = nodes[lftn];
		nodes[nn].child[1] = nodes[rgtn];	
		nodes[lftn].parent = nodes[nn];
		nodes[rgtn].parent = nodes[nn];
	}











   
}
