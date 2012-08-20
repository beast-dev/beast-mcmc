

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
import dr.evolution.util.Taxon;
import dr.inference.model.Parameter;
import dr.util.AlloppMisc;


/**
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
	private class MulLabNode {
		private MulLabNode parent;
		private MulLabNode child[];
		private double height;
		private boolean tetraroot;
		private double hybridheight;
		private FixedBitSet union;
		private ArrayList<Double> coalheights;
		private int nlineages;
		private String name;
		private int tippopindex;
		private int hybpopindex;
		private int rootpopindex;

		// dud constuctor
		MulLabNode() {
			height = -1.0;    
			tetraroot = false;
			hybridheight = -1.0;   
			coalheights = new ArrayList<Double>();
			name = "";
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
				formatter.format("%s ", name);
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
	}

	private class LegLink {
		public MulLabNode foot;
		public MulLabNode hip;
		public double footheight;
		private boolean done;

		public LegLink(MulLabNode hip, MulLabNode foot, double footheight) {
			this.hip = hip;
			this.foot = foot;
			this.footheight =  footheight;
			done = false;
		}
	}

	private class FootLinks {
		public List<LegLink> hips;
		public MulLabNode foot;

		public FootLinks(List<LegLink> hips) {
			this.hips = hips;
			if (hips.size() > 0) {
				foot = hips.get(0).foot;
				for (LegLink x : hips) {
					assert foot == x.foot;
				}
			}
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
			simpletree2mullabtree(root, 0);
			rootn = nextn - 1;
			diroot = nodes[rootn];
		}
		if (diroot != null &&  diroot.height <= 0.0) {
			System.err.println("AlloppMulLabTree constructor: bug");
		}

		assert (diroot != null  ||  tetratrees.length == 1);

		// Tetratrees. Two copies each.
		// When there is one leg (JOINED case) make an extra node to join
		// the two tetratrees into one.
		// While doing this, make a list of LegLinks = roots and footUnion's
		// which specify the join to diploid branch.
		List<LegLink> leglinks = new ArrayList<LegLink>();
		for (int i = 0; i < tetratrees.length; i++) {
			SimpleNode root = (SimpleNode) tetratrees[i].getRoot();
			simpletree2mullabtree(root, 0);
			int r0 = nextn - 1;
			nodes[r0].hybridheight = tetratrees[i].getHybridHeight();
			nodes[r0].tetraroot = true;
			simpletree2mullabtree(root, 1);
			int r1 = nextn - 1;
			nodes[r1].hybridheight = tetratrees[i].getHybridHeight();
			nodes[r1].tetraroot = true;
			if (tetratrees[i].getNumberOfLegs() == 1) {
				// in this case, the two tettrees just made join at
				// another node (the root in no-diploid case, and a node
				// which goes on join to diploid tree otherwise).
				// grjtodo code shared here and in simpletree2mullabtree()
				nodes[nextn].child = new MulLabNode[2];
				nodes[nextn].child[0] = nodes[r0];
				nodes[r0].parent = nodes[nextn];
				nodes[nextn].child[1] = nodes[r1];
				nodes[r1].parent = nodes[nextn];
				nodes[nextn].height = tetratrees[i].getSplitHeight();
				if (diroot != null) {
					FixedBitSet u = tetratrees[i].getFootUnion(0);
					double footheight = tetratrees[i].getFootHeight(0);
					leglinks.add(new LegLink(nodes[nextn], nodeOfUnionInSubtree(diroot, u), footheight));
				}
				nextn++;
			} else {
				assert diroot != null;
				FixedBitSet u0 = tetratrees[i].getFootUnion(0);
				FixedBitSet u1 = tetratrees[i].getFootUnion(1);
				double footheight0 = tetratrees[i].getFootHeight(0);
				double footheight1 = tetratrees[i].getFootHeight(1);
				leglinks.add(new LegLink(nodes[r0], nodeOfUnionInSubtree(diroot, u0), footheight0));
				leglinks.add(new LegLink(nodes[r1], nodeOfUnionInSubtree(diroot, u1), footheight1));
			}
		}
		if (diroot == null) {
			rootn = nextn - 1;
		}
		// Re-organise the LegLinks into per-diploid-branch lists (FootLinks)
		List<FootLinks> footlinkslist = new ArrayList<FootLinks>();
		for (LegLink x : leglinks) {
			if (!x.done) {
				List<LegLink> hips = new ArrayList<LegLink>();
				for (LegLink y : leglinks) {
					if (y.foot == x.foot) {
						hips.add(y);
						y.done = true;
					}
				}
				footlinkslist.add(new FootLinks(hips));
			}
		}
		// sort the feet in time order.
		for (FootLinks x : footlinkslist) {
			Collections.sort(x.hips, FOOTHEIGHT_ORDER);
		}

		// For all x, add nodes down from x.foot for all x.hips
		for (FootLinks x : footlinkslist) {
			MulLabNode f = x.foot;
			for (LegLink hip : x.hips) {
				MulLabNode h = hip.hip;
				nodes[nextn].parent = f.parent;
				nodes[nextn].child = new MulLabNode[2];
				nodes[nextn].child[0] = f;
				nodes[nextn].child[0].parent = nodes[nextn]; 
				nodes[nextn].child[1] = h;
				nodes[nextn].child[1].parent = nodes[nextn]; 
				if (nodes[nextn].parent.child[0] == f) {
					nodes[nextn].parent.child[0] = nodes[nextn];
				} else {
					assert nodes[nextn].parent.child[1] == f;
					nodes[nextn].parent.child[1] = nodes[nextn];
				}
				nodes[nextn].height = hip.footheight;
				f = nodes[nextn]; 
				nextn++;
			}
		}

		fillinUnionsInSubtree(nodes[rootn]);
		makesimpletree();
	}

	
	public String mullabTreeAsNewick() {
		return mullabSubtreeAsNewick(nodes[rootn], new String(""));
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
		MulLabNode node = nodeOfUnion(union);
		return (node.height <= height);
	}


	public void recordCoalescence(double height, FixedBitSet union) {
		MulLabNode node = nodeOfUnion(union);
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


	public double geneTreeInMULTreeLogLikelihood(boolean noDiploids, boolean twoDiploids) {
		fillinpopvals(noDiploids, twoDiploids);
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
			snodes[nextsn].setTaxon(new Taxon(mnode.name));
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
	private void simpletree2mullabtree(SimpleNode snode, int seq) {
		if (snode.isExternal()) {
			nodes[nextn].child = new MulLabNode[0];
			int sp = apsp.apspeciesId2index(snode.getTaxon().getId());
			int spseq = apsp.spandseq2spseqindex(sp, seq);
			nodes[nextn].union = new FixedBitSet(apsp.numberOfSpSeqs());
			nodes[nextn].union.set(spseq);
			nodes[nextn].name = snode.getTaxon().getId() + seq;
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
	}



	private void fillinUnionsInSubtree(MulLabNode node) {
		if (node.child.length > 0) {
			fillinUnionsInSubtree(node.child[0]);
			fillinUnionsInSubtree(node.child[1]);
			node.union = new FixedBitSet(apsp.numberOfSpSeqs());
			node.union.union(node.child[0].union);
			node.union.union(node.child[1].union);
		}
	}

	private MulLabNode nodeOfUnion(FixedBitSet x) {
		return nodeOfUnionInSubtree(nodes[rootn], x);
	}

	/* Searches subtree rooted at node for most tipward node 
	 * whose union contains x. If x is known to be a union of one of the nodes,
	 * it finds that node, so acts as a map union -> node
	 */
	private MulLabNode nodeOfUnionInSubtree(MulLabNode node, FixedBitSet x) {
		if (node.child.length == 0) {
			return node;
		}
		if (x.setInclusion(node.child[0].union)) {
			return nodeOfUnionInSubtree(node.child[0], x);
		} else if (x.setInclusion(node.child[1].union)) {
			return nodeOfUnionInSubtree(node.child[1], x);
		} else {
			return node;
		}
	}



	private String mullabSubtreeAsNewick(MulLabNode node, String newick) {
		if (node.child.length == 0) {
			int spseq = union2spseqindex(node.union);
			int sp = apsp.spseqindex2sp(spseq);
			int seq = apsp.spseqindex2seq(spseq);
			newick += apsp.apspeciesName(sp);
			assert seq < 10;
			newick += "0123456789".charAt(seq);
		} else {
			newick = "(" + mullabSubtreeAsNewick(node.child[0], newick) + "," + 
			mullabSubtreeAsNewick(node.child[1], newick) + ")";
		}
		return newick;
	}



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
	 * fillinpopvalsforspunionTwoDiploids() deal with a set of nodes
	 * with same species clade.
	 */
	private void fillinpopvals(boolean noDiploids, boolean twoDiploids) {
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
			} else if (twoDiploids) {
				p = fillinpopvalsforspunionTwoDiploids(unionarray, n0, n1, p);
			} else {
				// grjtodo tetraonly Not sure if TwoDiploids code generalises (not sure of model either)
				assert false;					
			}
			n0 = n1;
		}
		assert p == popvals.getDimension();
	}

		


	private int fillinpopvalsforspunionTwoDiploids(SpSqUnion[] unionarray, int n0, int n1, int p) {
		MulLabNode nodeset[] = new MulLabNode[3];
		int n = n1-n0;
		assert n == 1 || n==2 || n==3;
		// In two diploids case, get one of:
		// one diploid node (tip or root)
		// one foot node, where feet meet different diploid branches
		// two nodes from different tettrees
		// two nodes which are two feet of tettree meeting diploid branch
		// three nodes which are two tetroots and a leg-join.

		// get set of nodes with same species clade
		for (int i = n0; i < n1; i++) {
			nodeset[i-n0] = nodeOfUnion(unionarray[i].spsqunion);
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
			// either two nodes at the 'same place' in two tetratrees, or two
			// feet of a tetratree meeting same diploid branch
			assert nodeset[0].parent != null;
			assert nodeset[1].parent != null;
			if (nodeset[0].child.length == 2 &&
					(nodeset[0].child[0].tetraroot || nodeset[0].child[1].tetraroot)) {
				// if one node has a tetraroot child, assert that it has only one
				// such child (ie it is a foot) and that the other node is also a foot
				assert !nodeset[0].child[0].tetraroot || !nodeset[0].child[1].tetraroot;
				assert nodeset[1].child[0].tetraroot || nodeset[1].child[1].tetraroot;
				assert !nodeset[1].child[0].tetraroot || !nodeset[1].child[1].tetraroot;
				nodeset[0].rootpopindex = p;
				nodeset[1].rootpopindex = p+1;
				p += 2;
				// For each foot, copy rootpop from non-tetraroot child to tetraroot child
				for (int i = 0; i < 2; i++) {
					if (nodeset[i].child[0].tetraroot) {
						nodeset[i].child[0].rootpopindex = nodeset[i].child[1].rootpopindex;
					} else {
						nodeset[i].child[1].rootpopindex = nodeset[i].child[0].rootpopindex;
					}						
				}
			} else {
				// two nodes, not two feet. 'normal' case of matching nodes in two tetratrees.
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
			nodeset[i-n0] = nodeOfUnion(unionarray[i].spsqunion);
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
	 * this is for sorting the links from higher ploidy tree  to a single
	 * branch within a lower ploidy tree.
	 */
	static final Comparator<LegLink> FOOTHEIGHT_ORDER = new Comparator<LegLink>() {
		public int compare(LegLink a, LegLink b) {
			if (a.footheight == b.footheight) {
				return 0;
			} else {
				return (a.footheight > b.footheight) ? 1 : -1;
			}
		}
	};    



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
