

package dr.evomodel.speciation;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.Formatter;
import java.util.List;
import java.util.Locale;
import java.util.Stack;

import jebl.util.FixedBitSet;
import dr.evolution.tree.SimpleNode;
import dr.evolution.tree.SimpleTree;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxon;
import dr.inference.model.Parameter;
import dr.util.AlloppMisc;


/**
 * An AlloppMulLabTree represents the species network as single 
 * binary tree with tips that can be multiply labelled with species. 
 * 
 * @author Graham Jones
 *         Date: 13/09/2011
 */


/*
 * An AlloppMulLabTree represents the species network as single 
 * binary tree with tips that can be multiply labelled with species. 
 * 
 * classes LegLink and FootLinks are for gathering and organising
 * the links between trees of different ploidy, so that the 
 * rootward-pointing legs can become tipward-pointing branches.
 * 
 * SpSqUnion is used for sorting the nodes in an AlloppMulLabTree. It is 
 * used by Comparator SPUNION_ORDER, and hence indirectly by
 * fillinpopvals().
 * 
 * class BranchPopulationAndLineages records the information needed
 * to calculate the probability of coalescences in a single branch of the
 * AlloppMulLabTree.
 * 
 */


/* nodes[], rootn implement the tree; nextn is for building it
 * 
 * apsp references the (species, indivs, sequences) structure.
 * 
 * popvals references the population parameters. fillinpopvals() assigns
 * them to branches.
 * 
 * simptree is so that AlloppSpeciesNetworkModel, which contains a AlloppMulLabTree,
 * can implement the Tree interface.
 */
public class AlloppMulLabTree  {
	private MulLabNode[] nodes;
	private int rootn;
	private int nextn;
	private AlloppSpeciesBindings apsp;
	private Parameter popvals;
	public SimpleTree simptree;
	
	/*
	 * parent, child[] join the nodes into a binary tree.
	 * 
	 * height is time into past
	 * 
	 * popsize is population at start of branch, and for tips, also at end.
	 * 
	 * union is a set of species for a single choice of sequence copy from
	 * each individual of the species. There is one bit for each Taxon like 
	 * "c0" or "c1"
	 */	
	private class MulLabNode extends AlloppNode.Abstract implements AlloppNode {
		private MulLabNode parent;
		private MulLabNode child[];
		private double height;
		private boolean tetraroot;
		private double hybridheight;
		private FixedBitSet union;
		private ArrayList<Double> coalheights;
		private int nlineages;
		private Taxon taxon;
		private int tippopindex;
		private int hybpopindex;
		private int rootpopindex;

		// dud constuctor
		MulLabNode() {
			parent = null;
			child = new MulLabNode[0];
			height = -1.0;    
			tetraroot = false;
			hybridheight = -1.0;   
			coalheights = new ArrayList<Double>();
			taxon = new Taxon("");
			tippopindex = -1;
			hybpopindex	= -1;
			rootpopindex = -1;
		}


		public double tippop() {
			return popvals.getParameterValue(tippopindex);
		}

		public double hybpop() {
			return popvals.getParameterValue(hybpopindex);
		}

		public double rootpop() {
			return popvals.getParameterValue(rootpopindex);
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
			double tippop = (tippopindex >= 0) ? tippop() : -1.0;
			formatter.format("%s %s ", AlloppMisc.nonnegIntIn2Chars(tippopindex), AlloppMisc.nonnegIn8Chars(tippop));
			double hybpop = (hybpopindex >= 0) ? hybpop() : -1.0;
			formatter.format("%s %s ", AlloppMisc.nonnegIntIn2Chars(hybpopindex), AlloppMisc.nonnegIn8Chars(hybpop));
			double rootpop = (rootpopindex >= 0) ? rootpop() : -1.0;
			formatter.format("%s %s ", AlloppMisc.nonnegIntIn2Chars(rootpopindex), AlloppMisc.nonnegIn8Chars(rootpop));
			formatter.format("%s ", tetraroot ? "tetroot" : "       ");
			formatter.format("%s ", AlloppMisc.nonnegIn8Chars(hybridheight));		
			formatter.format("%3d  ", nlineages);
			for (int c = 0; c < coalheights.size(); c++) {
				formatter.format(AlloppMisc.nonnegIn8Chars(coalheights.get(c)) + ",");
			}
			return s.toString();
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
		public double getHeight() {
			return height;
		}


		@Override
		public FixedBitSet getUnion() {
			return union;
		}


		@Override
		public void setChild(int ch, AlloppNode newchild) {
			child[ch] = (MulLabNode) newchild;
            newchild.setAnc(this);
		}


		@Override
		public void setAnc(AlloppNode anc) {
			parent = (MulLabNode) anc;
		}

		@Override
		public Taxon getTaxon() {
			return taxon;
		}


		@Override
		public void setTaxon(String name) {
			taxon = new Taxon(name);
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
			child = new MulLabNode[2];
			child[0] = (MulLabNode)c0;
			child[0].setAnc(this);
			child[1] = (MulLabNode)c1;
			child[1].setAnc(this);
			}
		
		
	
	}


	private class SpSqUnion {
		public FixedBitSet spsqunion;
		public FixedBitSet spunion;
		public SpSqUnion(FixedBitSet spsqunion) {
			this.spsqunion = spsqunion;
			spunion = apsp.spsqunion2spunion(spsqunion);
		}
	}


	private class PopulationAndLineages {
		public double t[];
		public double tippop;
		public double rootpop;
		public int tipnlin;

		public PopulationAndLineages(double[] t2, double tippop, double rootpop, 
				int tipnlin) {
			this.t = t2;
			this.tippop = tippop;
			this.rootpop = rootpop;
			this.tipnlin = tipnlin;
		}

		public double populationAt(double x) {
			final double begt = t[0];
			final double endt = t[t.length - 1];
			return ((endt-x)*tippop + (x-begt)*rootpop) / (endt-begt);
		}
	}
		
	    
	   
	    

	/*
	 * This constructor makes a single multiply labelled tree from the set
	 * of homoploid SimpleTrees which is passed to it. It is called directly
	 * by testing code.
	 */
	AlloppMulLabTree(AlloppLeggedTree[][] trees, AlloppSpeciesBindings apsp, Parameter popvals) {
		this.apsp = apsp;
		this.popvals = popvals;
		// Count tips for the tree to be made
		int ntips = 0;
		AlloppLeggedTree[] ditrees = trees[AlloppSpeciesNetworkModel.DITREES];
		AlloppLeggedTree[]  tetratrees = trees[AlloppSpeciesNetworkModel.TETRATREES];
		if (ditrees.length > 0) {
			ntips += ditrees[0].getTaxonCount(); // never more than one diploid tree
		}
		for (int i = 0; i < tetratrees.length; i++) {
			ntips += 2 * tetratrees[i].getTaxonCount();
		}
		// Make array of dud nodes
		nodes = new MulLabNode[2 * ntips - 1];
		for (int i = 0; i < nodes.length; i++) {
			nodes[i] = new MulLabNode();
		}
		// Copy homoploid trees into array
		nextn = 0;
		// Ditree if it exists
		MulLabNode diroot = null;
		if (ditrees.length > 0) {
			SimpleNode root = (SimpleNode) ditrees[0].getRoot();
			nextn = AlloppNode.Abstract.simpletree2allopptree(apsp, nodes, nextn, root, false, 0);
			rootn = nextn - 1;
			diroot = nodes[rootn];
		}
		if (diroot != null &&  diroot.height <= 0.0) {
			System.err.println("AlloppMulLabTree constructor: bug");
		}
		assert (diroot != null  ||  tetratrees.length == 1);
		if (diroot != null)
		  { diroot.fillinUnionsInSubtree(apsp.numberOfSpSeqs()); }   // added 2012-03-22 when doing DipHist

		List<AlloppLegLink> leglinks = new ArrayList<AlloppLegLink>();
		for (int i = 0; i < tetratrees.length; i++) {
			// Tetratrees. Two copies each.
			SimpleNode root = (SimpleNode)tetratrees[i].getRoot();
			nextn = AlloppNode.Abstract.simpletree2allopptree(apsp, nodes, nextn, root, true, 0);
			int r0 = nextn - 1;
			nodes[r0].hybridheight = tetratrees[i].getHybridHeight();
			nodes[r0].tetraroot = true;
			nextn = AlloppNode.Abstract.simpletree2allopptree(apsp, nodes, nextn, root, true, 1);
			int r1 = nextn - 1;
			nodes[r1].hybridheight = tetratrees[i].getHybridHeight();
			nodes[r1].tetraroot = true;
			nextn = AlloppNode.Abstract.collectLegLinksAndMakeLegJoin(leglinks, tetratrees[i], diroot, nodes, nextn, r0, r1);
		}
		if (diroot == null) {
			rootn = nextn - 1;
		}
		AlloppNode.Abstract.convertLegLinks(nodes, nextn, leglinks);
		nodes[rootn].fillinUnionsInSubtree(apsp.numberOfSpSeqs());
		makesimpletree();
	}

	
	
	
	public String mullabTreeAsNewick() {
		String s = Tree.Utils.uniqueNewick(simptree, simptree.getRoot());
		return s;
	}

	public String asText() {
		String header = "topology             height         union         []  tippop  []  hybpop  [] rootpop  tetroot   hybhgt nlin coalheights" + System.getProperty("line.separator");

		String s = "";
		Stack<Integer> x = new Stack<Integer>();
		return header + subtreeAsText(nodes[rootn], s, x, 0, "");
	}


	public void clearCoalescences() {
		clearSubtreeCoalescences(nodes[rootn]);
	}	

	public void recordLineageCounts() {
		recordSubtreeLineageCounts(nodes[rootn]);
	}

	public boolean coalescenceIsCompatible(double height, FixedBitSet union) {
		MulLabNode node = (MulLabNode) nodes[rootn].nodeOfUnionInSubtree(union);
		return (node.height <= height);
	}


	public void recordCoalescence(double height, FixedBitSet union) {
        MulLabNode node = (MulLabNode) nodes[rootn].nodeOfUnionInSubtree(union);
		assert (node.height <= height);
		while (node.parent != null  &&  node.parent.height <= height) {
			node = node.parent;
		}
		node.coalheights.add(height);
	}

	public void sortCoalescences() {
		for (MulLabNode node : nodes) {
			Collections.sort(node.coalheights);
		}
	}


	public double geneTreeInMULTreeLogLikelihood(boolean noDiploids) {
		fillinpopvals(noDiploids);
		//System.out.println(asText());
		return geneTreeInMULSubtreeLogLikelihood(nodes[rootn]);
	}


	
	
	
/*
 * 
 * ***************************************************
 * 	
 */
	
	
	
	
	
	
	private void makesimpletree() {
    	SimpleNode[] snodes = new SimpleNode[nodes.length];
    	for (int n = 0; n < nodes.length; n++) {
    		snodes[n] = new SimpleNode();
    	}
		makesimplesubtree(snodes, 0, nodes[rootn]);
		simptree = new SimpleTree(snodes[nodes.length-1]);
	}
 	
	
	private int makesimplesubtree(SimpleNode[] snodes, int nextsn, MulLabNode mnode) {
		if (mnode.child.length == 0) {
			snodes[nextsn].setTaxon(new Taxon(mnode.taxon.getId()));
		} else {
			nextsn = makesimplesubtree(snodes, nextsn, mnode.child[0]);
			int subtree0 = nextsn-1;
			nextsn = makesimplesubtree(snodes, nextsn, mnode.child[1]);
			int subtree1 = nextsn-1;
			snodes[nextsn].addChild(snodes[subtree0]);
			snodes[nextsn].addChild(snodes[subtree1]);
		}
		snodes[nextsn].setHeight(mnode.height);
		return nextsn+1;
	}
	


	private String subtreeAsText(MulLabNode node, String s, Stack<Integer> x, int depth, String b) {
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



	/* Recursively copies the topology from subtree rooted at node into
	 * mullabtree implemented as array nodes[].
	 * 
	 * Fill in the unions at the tips: convert Taxon Id, like "c" into 
	 * an index for a species, then combined with seq to make an index 
	 * to be set in union.
	 */
    /*  ood
	private void simpletree2mullabtree(SimpleNode snode, int seq) {
		if (snode.isExternal()) {
			nodes[nextn].child = new MulLabNode[0];
			nodes[nextn].union = apsp.speciesseqToTipUnion(snode.getTaxon(), seq);
			nodes[nextn].taxon = new Taxon(snode.getTaxon().getId() + seq);
		} else {
			simpletree2mullabtree(snode.getChild(0), seq);
			int c0 = nextn - 1;
			simpletree2mullabtree(snode.getChild(1), seq);
			int c1 = nextn - 1;
			nodes[nextn].child = new MulLabNode[2];
			nodes[nextn].child[0] = nodes[c0];
			nodes[c0].parent = nodes[nextn];
			nodes[nextn].child[1] = nodes[c1];
			nodes[c1].parent = nodes[nextn];
		}
		nodes[nextn].height = snode.getHeight();
		nextn++;
	} */




	private void clearSubtreeCoalescences(MulLabNode node) {
		for (int i = 0; i < node.child.length; i++) {
			clearSubtreeCoalescences(node.child[i]);
		}
		node.coalheights.clear();
	}


	private void recordSubtreeLineageCounts(MulLabNode node) {
		if (node.child.length == 0) {
			node.nlineages = apsp.nLineages(apsp.spseqindex2sp(union2spseqindex(node.union)));
		} else {
			node.nlineages = 0;
			for (int i = 0; i < node.child.length; i++) {
				recordSubtreeLineageCounts(node.child[i]);				
				node.nlineages += node.child[i].nlineages - node.child[i].coalheights.size();
			}
		}			
	}		



	/*
	 * This copies population values in the Parameter popvalues 
	 * to nodes in the AlloppMulLabTree. The population values are
	 * per-species-clade (per-branch in network), but of course more than 
	 * one node in AlloppMulLabTree may correspond to the same species.
	 * 
	 * The other complications are that tips are different from internal 
	 * nodes, and that nodes which roots of tetratrees or just below,
	 * as well as the root are special cases.
	 * 
	 * It collects unions (which represent sets whose elements 
	 * identify a species and a sequence) from the nodes and then
	 * sorts them primarily using identities of the species, so
	 * that sets of node with same species clade are grouped together. The sort 
	 * also puts the node sets corresponding to tips first in the array and sorts
	 * nodes within node sets in a well-defined way.
	 * This mainly does what is required, since nodes with the same
	 * species clade are treated the same.
	 * 
	 * fillinpopvalsforspunionNoDiploids() and
	 * fillinpopvalsforspunion() deal with a set of nodes
	 * with same species clade.
	 */
	private void fillinpopvals(boolean noDiploids) {
		ArrayList<SpSqUnion> unionarraylist = new ArrayList<SpSqUnion>();
		for (int n = 0; n < nodes.length; n++) {
			unionarraylist.add(new SpSqUnion(nodes[n].union));
		}
		Collections.sort(unionarraylist, SPUNION_ORDER);
		SpSqUnion[] unionarray = new SpSqUnion[unionarraylist.size()];
		unionarray = unionarraylist.toArray(unionarray);
		int p = 0;
		for (int n0 = 0; n0 < unionarray.length; ) {
			int n1 = n0+1; 
			while (n1 < unionarray.length && unionarray[n1].spunion.equals(unionarray[n0].spunion)) {
				n1++;
			}
			if (noDiploids) {
				p = fillinpopvalsforspunionNoDiploids(unionarray, n0, n1, p);
			} else {
				p = fillinpopvalsforspunion(unionarray, n0, n1, p);			
			}
			n0 = n1;
		}
		if ( p != popvals.getDimension()) {
			System.out.println("BUG in fillinpopvals()");
		}
		assert p == popvals.getDimension();
	}

		


	private int fillinpopvalsforspunion(SpSqUnion[] unionarray, int n0, int n1, int p) {
		MulLabNode nodeset[] = new MulLabNode[3];
		int n = n1-n0;
		assert n == 1 || n==2 || n==3;
		// Get one of:
		// one diploid node (tip or root)
		// one foot node, where feet meet different diploid branches
		// two nodes from different tettrees
		// two nodes which are two feet of tettree meeting diploid branch
		// two nodes where one is a foot of a tetratree other contains some diploids and the partner tetratree
		// three nodes which are two tetroots and a leg-join.

		// get set of nodes with same species clade
		for (int i = n0; i < n1; i++) {
			nodeset[i-n0] = (MulLabNode) nodes[rootn].nodeOfUnionInSubtree(unionarray[i].spsqunion);
		}
		// set all pop values to dud values
		for (int i = 0; i < n; i++) {
			nodeset[i].tippopindex = -1;
			nodeset[i].hybpopindex = -1;
			nodeset[i].rootpopindex = -1;

		}

		if (n == 3) {
			// must be two tetraroots and one leg-join
			assert nodeset[0].tetraroot  &&  nodeset[0].parent != null;
			assert nodeset[1].tetraroot  &&  nodeset[1].parent != null;
			assert nodeset[2].child.length == 2;
			assert nodeset[2].child[0] == nodeset[0]  ||  nodeset[2].child[0] == nodeset[1];
			assert nodeset[2].child[1] == nodeset[0]  ||  nodeset[2].child[1] == nodeset[1];
			if (nodeset[0].child.length == 0) {
				for (int i = 0; i < 2; i++) {
					nodeset[i].tippopindex = p;
					nodeset[i].hybpopindex = p+1;
					nodeset[i].rootpopindex = p+2;	
				}
				p += 3;
			} else {
				for (int i = 0; i < 2; i++) {
					nodeset[i].hybpopindex = p;
					nodeset[i].rootpopindex = p+1;	
				}
				p += 2;
			}
			// copy to leg-join
			nodeset[2].rootpopindex = nodeset[0].rootpopindex;
		} else if (n == 2) {
			// either two nodes at the 'same place' in two matching tetratrees, or 
			// the two feet of same tetratree joining same diploid branch
			// one foot of a tetratree joining some diploids and the partner tetratree
			assert nodeset[0].parent != null;
			assert nodeset[1].parent != null;
			boolean [] isfoot = new boolean [2];
			for (int i=0; i<2; i++) {
				isfoot[i] = nodeset[i].child.length == 2  &&  (nodeset[i].child[0].tetraroot || nodeset[i].child[1].tetraroot);
			}
			if (isfoot[0]  || isfoot[1]) {
				// one or two feet.
				// assert that neither node has two tetraroot children (which is joined legs case with n== 3)
				assert !nodeset[0].child[0].tetraroot || !nodeset[0].child[1].tetraroot;
				assert !nodeset[1].child[0].tetraroot || !nodeset[1].child[1].tetraroot;
				nodeset[0].rootpopindex = p;
				nodeset[1].rootpopindex = p+1;
				p += 2;
				// For each foot, copy rootpop from non-tetraroot child to tetraroot child
				for (int i=0; i<2; i++) {
					if (isfoot[i]) {
						if (nodeset[i].child[0].tetraroot) {
							nodeset[i].child[0].rootpopindex = nodeset[i].child[1].rootpopindex;
						} else {
							nodeset[i].child[1].rootpopindex = nodeset[i].child[0].rootpopindex;
						}
					}					
				}
			} else {
				// two nodes, neither a foot, so they are matching nodes in two tetratrees.
				// If it is a tip, add tippop.
				// Whether or not it is a tip, add one more param, as hybpop if it is
				// a tetraroot, or rootpop if it isn't
				if (nodeset[0].child.length == 0) {
					for (int i = 0; i < 2; i++) {
						nodeset[i].tippopindex = p;
					}
					p += 1;
				} 
				if (nodeset[0].tetraroot) {
					assert nodeset[1].tetraroot;
					for (int i = 0; i < 2; i++) {
						nodeset[i].hybpopindex = p;	
					}
					p += 1;
				} else {					
					for (int i = 0; i < 2; i++) {
						nodeset[i].rootpopindex = p;	
					}
					p += 1;
				}
			}
		} else {
			assert n==1;
			// a diploid node, may be the root, maybe a tip; or a foot.
			// Nothing happens at root. otherwise: tips get tippop, all get rootpop
			if (nodeset[0].parent != null) {
				if (nodeset[0].child.length == 0) {
					nodeset[0].tippopindex = p;
					p += 1;
				}
				nodeset[0].rootpopindex = p;	
				p += 1;
				// If it is a  foot, copy rootpop from non-tetraroot child to tetraroot child
				if (nodeset[0].child.length == 2 &&
						(nodeset[0].child[0].tetraroot || nodeset[0].child[1].tetraroot)) {
					if (nodeset[0].child[0].tetraroot) {
						nodeset[0].child[0].rootpopindex = nodeset[0].child[1].rootpopindex;
					} else {
						nodeset[0].child[1].rootpopindex = nodeset[0].child[0].rootpopindex;
					}							
				}
			}					
		}
		return p;
	}

		
		

	/*
	 * In NODIPLOIDS case, there are 5 kinds of nodes. 
	 * 
	 * 1. tip AND tetroot (only in case of a single tetraploid)
	 * 2. tip
	 * 3. NOT tip NOT tetroot NOT diproot (internal branch when >=3 tetraploids)
	 * 4. tetroot
	 * 5. diproot
	 * 
	 */
	private int fillinpopvalsforspunionNoDiploids(SpSqUnion[] unionarray, int n0, int n1, int p) {
		MulLabNode nodeset[] = new MulLabNode[3];
		int n = n1-n0;
		assert n == 3  || n == 2;
		// In no diploids case, either get two nodes from different tettrees
		// or two tetroots and the root.

		// get set of nodes with same species clade
		for (int i = n0; i < n1; i++) {
			nodeset[i-n0] = (MulLabNode) nodes[rootn].nodeOfUnionInSubtree(unionarray[i].spsqunion);

		}
		// set all pop values to dud values
		for (int i = 0; i < n; i++) {
			nodeset[i].tippopindex = -1; 
			nodeset[i].hybpopindex = -1; 
			nodeset[i].rootpopindex = -1; 
		}

		if (n == 3) {
			assert nodeset[0].tetraroot  &&  nodeset[0].parent != null;
			assert nodeset[1].tetraroot  &&  nodeset[1].parent != null;
			assert nodeset[2].parent == null;
			if (nodeset[0].child.length == 0) {
				for (int i = 0; i < 2; i++) {
					nodeset[i].tippopindex = p;
					nodeset[i].hybpopindex = p+1;
					nodeset[i].rootpopindex = p+2;	
				}
				p += 3;
			} else {
				for (int i = 0; i < 2; i++) {
					nodeset[i].hybpopindex = p;
					nodeset[i].rootpopindex = p+1;	
				}
				p += 2;
			}
		} else {
			assert n == 2;
			assert nodeset[0].parent != null;
			assert nodeset[1].parent != null;
			if (nodeset[0].child.length == 0) {
				for (int i = 0; i < 2; i++) {
					nodeset[i].tippopindex = p;
					nodeset[i].rootpopindex = p+1;	
				}
				p += 2;

			} else {
				for (int i = 0; i < 2; i++) {
					nodeset[i].rootpopindex = p;	
				}
				p += 1;
			}
		}
		return p;
	}



	/*
	 * Visits each node in MULtree and accumulates LogLikelihood
	 */
	private double geneTreeInMULSubtreeLogLikelihood(MulLabNode node) {
		double loglike = 0.0;
		for (int i = 0; i < node.child.length; i++) {
			loglike += geneTreeInMULSubtreeLogLikelihood(node.child[i]);
		}
		if (apsp.SpeciesWithinPloidyLevel(2).length == 0) {
			loglike += branchLLInMULtreeNoDiploids(node);
		} else {
			loglike += branchLLInMULtreeTwoDiploids(node);
		}

		return loglike;
	}


	/*
	 * Does likelihood calculation for a single node in the case
	 * of two diploids.
	 */
	private double branchLLInMULtreeTwoDiploids(MulLabNode node) {
		double loglike = 0.0;

		double tippop = 0.0;
		if (node.child.length == 0) {
			tippop = node.tippop();
		} else {
			tippop = node.child[0].rootpop() + node.child[1].rootpop();
		}

		PopulationAndLineages pal;
		double t[];
		if (node.tetraroot) {
			// since hybridization
			int nsince = 0;
			for ( ; nsince < node.coalheights.size() && 
			node.coalheights.get(nsince) < node.hybridheight; nsince++) {}
			t = new double[nsince + 2];
			t[0] = node.height;
			for (int i = 0; i < nsince; i++) {
				t[i+1] = node.coalheights.get(i);
			}
			t[t.length-1] = node.hybridheight;
			pal = new PopulationAndLineages(t, tippop, node.hybpop(), node.nlineages);
			loglike += limbLogLike(pal);
			// before hybridization
			int nbefore = node.coalheights.size() - nsince;
			t = new double[nbefore + 2];
			t[0] = node.hybridheight;
			for (int i = 0; i < nbefore; i++) {
				t[i+1] = node.coalheights.get(nsince+i);
			}				
			t[t.length-1] = node.parent.height;

			pal = new PopulationAndLineages(t, node.rootpop(), node.rootpop(), node.nlineages - nsince);
			loglike += limbLogLike(pal);				
		} else if (node.parent == null) {
			t = new double[node.coalheights.size() + 2];
			t[0] = node.height;
			t[t.length-1] = apsp.maxGeneTreeHeight();
			for (int i = 0; i < node.coalheights.size(); i++) {
				t[i+1] = node.coalheights.get(i);
			}
			pal = new PopulationAndLineages(t, tippop, tippop, node.nlineages);
			loglike += limbLogLike(pal);					
		} else {
			t = new double[node.coalheights.size() + 2];
			t[0] = node.height;
			t[t.length-1] = node.parent.height;
			for (int i = 0; i < node.coalheights.size(); i++) {
				t[i+1] = node.coalheights.get(i);
			}
			pal = new PopulationAndLineages(t, tippop, node.rootpop(), node.nlineages);
			loglike += limbLogLike(pal);				
		}

		return loglike;
	}


	/*
	 * Does likelihood calculation for a single node in the case
	 * of no diploids
	 */
	private double branchLLInMULtreeNoDiploids(MulLabNode node) {
		double loglike = 0.0;

		double tippop = 0.0;
		if (node.child.length == 0) {
			tippop = node.tippop();
		} else {
			tippop = node.child[0].rootpop() + node.child[1].rootpop();
		}

		PopulationAndLineages pal;
		double t[];
		if (node.tetraroot) {
			// since hybridization
			int nsince = 0;
			for ( ; nsince < node.coalheights.size() && 
			node.coalheights.get(nsince) < node.hybridheight; nsince++) {}
			t = new double[nsince + 2];
			t[0] = node.height;
			for (int i = 0; i < nsince; i++) {
				t[i+1] = node.coalheights.get(i);
			}
			t[t.length-1] = node.hybridheight;
			pal = new PopulationAndLineages(t, tippop, node.hybpop(), node.nlineages);
			loglike += limbLogLike(pal);
			// before hybridization
			int nbefore = node.coalheights.size() - nsince;
			t = new double[nbefore + 2];
			t[0] = node.hybridheight;
			for (int i = 0; i < nbefore; i++) {
				t[i+1] = node.coalheights.get(nsince+i);
			}				
			t[t.length-1] = node.parent.height;
			double dippop = node.rootpop();
			pal = new PopulationAndLineages(t, dippop, dippop, node.nlineages - nsince);
			loglike += limbLogLike(pal);				
		} else if (node.parent == null) {
			t = new double[node.coalheights.size() + 2];
			t[0] = node.height;
			t[t.length-1] = apsp.maxGeneTreeHeight();
			for (int i = 0; i < node.coalheights.size(); i++) {
				t[i+1] = node.coalheights.get(i);
			}
			pal = new PopulationAndLineages(t, tippop, tippop, node.nlineages);
			loglike += limbLogLike(pal);					
		} else {
			t = new double[node.coalheights.size() + 2];
			t[0] = node.height;
			t[t.length-1] = node.parent.height;
			for (int i = 0; i < node.coalheights.size(); i++) {
				t[i+1] = node.coalheights.get(i);
			}
			pal = new PopulationAndLineages(t, tippop, node.rootpop(), node.nlineages);
			loglike += limbLogLike(pal);				
		}
		return loglike;
	}





	/*
	 * limbLogLike calculates the log-likelihood for 
	 * the coalescences at t[1],t[2],...t[k] within a limb
	 * from t[0] to t[k+1]. ('limb' means a branch or part of one.)
	 */
	private double limbLogLike(PopulationAndLineages pal) {
		double loglike = 0.0;
		int k = pal.t.length - 2;
		for (int i = 1; i <= k; i++) {
			loglike -= Math.log(pal.populationAt(pal.t[i]));
		}
		for (int i = 0; i <= k; i++) {
			final double y = (pal.tipnlin-i) * (pal.tipnlin-i-1) / 2 ;
			final double z = limbLinPopIntegral(pal, pal.t[i], pal.t[i+1]);
			loglike -= y * z;
		}
		return loglike;
	}



	// integral from t0 to t1 of (endt-begt)/((endt-x)begPop + (x-begt)endPop) 
	// with respect to x
	private double limbLinPopIntegral(PopulationAndLineages b, double t0, double t1) {
		final double begt = b.t[0];
		final double endt = b.t[b.t.length-1];
		final double d = b.rootpop - b.tippop;
		final double c = endt * b.tippop - begt * b.rootpop;
		final double x = Math.abs(d / b.tippop);
		if (x > 0.001) {
			return ((endt - begt) / d) * Math.log((c + d * t1) / (c + d * t0));
		} else {
			double y = d * (t1 - t0) / (c + d * t0);
			double ys = (1.0 - y/2.0 + y*y * (1.0/3.0 - y/4.0 + y*y * (1.0/5.0 - y/6.0)));
			return ((endt - begt) * (t1 - t0) / (c + d * t0)) * ys; 
		}
	}		




	private static int union2spseqindex(FixedBitSet union) {
		assert union.cardinality() == 1;
		return union.nextOnBit(0);
	}




	/*
	 * This is for ordering the unions in the nodes of the AlloppMulLabTree.
	 * Those unions are of (species, sequence) pairs.
	 * 
	 * The comparator sorts the unions of (species, sequence) pairs (SpSqUnions)
	 * so that all unions containing the same set of species (ignoring sequence)
	 * are grouped together. Call the sets of SpSqUnions for the same species a
	 * `group'. There can be 1,2 or 3 SpSqUnions in a group. 
	 * 
	 * The groups are sorted in order of increasing number of species (clade size).
	 * All groups for a single species (a tip in the network) come first, then 
	 * those groups for two species, and so on to the root for all species.
	 * For groups that have equal numbers of species, a lexicographical
	 * ordering using species indices is used.
	 * 
	 * Within each group, species and sequence information is used to sort the 1 to 3
	 * SpSqUnions. The size of the `clade' of (species, sequence) pairs is used
	 * first in the comparison, which ensures that the three nodes with the same species
	 * - corresponding to two roots of tetratrees in the AlloppMulLabTree plus a leg-join - 
	 * are ordered so that the two roots come first.
	 * 
	 */

	static final Comparator<SpSqUnion> SPUNION_ORDER = new Comparator<SpSqUnion>() {
		public int compare(SpSqUnion a, SpSqUnion b) {
			int ac = a.spunion.cardinality();
			int bc = b.spunion.cardinality();
			if (ac != bc) {
				return ac-bc;
			} else {
				int an = a.spunion.nextOnBit(0);
				int bn = b.spunion.nextOnBit(0);
				while (an >= 0 || bn >= 0) {
					if (an != bn) {
						return an-bn;
					}
					an = a.spunion.nextOnBit(an+1);
					bn = b.spunion.nextOnBit(bn+1);
				}
				// spunions compare equal; do spsqunions
				ac = a.spsqunion.cardinality();
				bc = b.spsqunion.cardinality();
				if (ac != bc) {
					return ac-bc;
				} else {
					an = a.spsqunion.nextOnBit(0);
					bn = b.spsqunion.nextOnBit(0);
					while (an >= 0 || bn >= 0) {
						if (an != bn) {
							return an-bn;
						}
						an = a.spsqunion.nextOnBit(an+1);
						bn = b.spsqunion.nextOnBit(bn+1);
					}				
					return 0;
				}
			}
		}
	};



	
}
