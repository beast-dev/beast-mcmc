package dr.evomodel.speciation;



import dr.evolution.tree.*;
import dr.evolution.util.Taxon;
import dr.evolution.util.Units;
import dr.evomodelxml.speciation.AlloppSpeciesNetworkModelParser;
import dr.inference.loggers.LogColumn;
import dr.inference.loggers.Loggable;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.inference.operators.OperatorFailedException;
import dr.inference.operators.Scalable;
import dr.math.MathUtils;
import dr.util.AlloppMisc;
import dr.util.Author;
import dr.util.Citable;
import dr.util.Citation;
import jebl.util.FixedBitSet;
import java.util.*;
import java.util.logging.Logger;

import test.dr.evomodel.speciation.AlloppSpeciesNetworkModelTEST;


/**
 * 
 * Implements an allopolyploid species network as a collection of `trees with legs'.
 * 
 * @author Graham Jones
 *         Date: 19/04/2011
 */


/*
 * class AlloppSpeciesNetworkModel
 * 
 * Implements the species network as a collection of `trees with legs'.
 * and converts this representation into a multiply labelled 
 * binary tree.
 * 
 * General idea is that the network is easiest to change (eg detach
 * and re-attach tetraploid subtrees) while likelihood calculations
 * are easiest to do in the multiply labelled tree.
 * 
 * The individual `trees with legs' are implemented by AlloppLeggedTree's.
 * 
 * *********************
 * 
 * apsp is a reference to the AlloppSpeciesBindings which knows how 
 * species are made of individuals and individuals are made of taxa,
 * and which contains the list of gene trees.
 * 
 * trees[][] represents the network as a set of homoploid trees
 * 
 * mullabtree represents the network as single tree with tips that
 * can be multiply labelled with species. 
 * 
 */
//  grjtodo JH's SpeciesTreeModel implements 
// MutableTree, TreeTraitProvider, TreeLogger.LogUpon, Scalable
// not clear how much of those sensible here.
// AlloppLeggedTree implements MutableTree, TreeLogger.LogUpon.
// Nothing so far does TreeTraitProvider.
public class AlloppSpeciesNetworkModel extends AbstractModel implements
		Scalable, Units, Citable, Loggable {

	private final AlloppSpeciesBindings apsp;
	private AlloppLeggedTree[][] trees;
	private AlloppLeggedTree[][] oldtrees;
	private MulLabTree mullabtree;
		
	// grjtodo this is public, copying JH - is that necessary? TreeNodeSlide accesses it.
	// 2011-06-30 parser accesses it too. 
	// Parameter or Parameter.Default ?? (a Java thing I don't get)
    public final Parameter popvalues;

	public final static boolean DBUGTUNE = false;


	public enum LegType {
		NONE, TWOBRANCH, ONEBRANCH, JOINED, NODIPLOIDS
	};

	public final static int DITREES = 0;
	public final static int TETRATREES = 1;
	public final static int NUMBEROFPLOIDYLEVELS = 2;

	

/* ***************** subclasses *********************************** */    
	
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
			return popvalues.getParameterValue(tippopindex);
		}
		
		public double hybpop() {
			return popvalues.getParameterValue(hybpopindex);
		}
		
		public double rootpop() {
			return popvalues.getParameterValue(rootpopindex);
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

	
	

	/*
	 * A MulLabTree represents the species network as single 
	 * binary tree with tips that can be multiply labelled with species. 
	 * 
	 * classes LegLink and FootLinks are for gathering and organising
	 * the links between trees of different ploidy, so that the 
	 * rootward-pointing legs can become tipward-pointing branches.
	 * 
	 * SpSqUnion is used for sorting the nodes in a MulLabTree. It is 
	 * used by Comparator SPUNION_ORDER, and hence indirectly by
	 * fillinpopvals().
	 * 
	 * class BranchPopulationAndLineages records the information needed
	 * to calculate the probability of coalescences in a single branch of the
	 * MulLabTree.
	 * 
	 */
	private class MulLabTree {
		private MulLabNode[] nodes;
		private int rootn;
		private int nextn;

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
		 * of homoploid SimpleTrees.
		 */
	    MulLabTree() {
			 this(trees);
		 }
	    
		
		/*
		 * This constructor makes a single multiply labelled tree from the set
		 * of homoploid SimpleTrees which is passed to it. It is called directly
		 * by testing code.
		 */
		 MulLabTree(AlloppLeggedTree[][] trees) {
			// Count tips for the tree to be made
			int ntips = 0;
		    if (trees[DITREES].length > 0) {
		    	ntips += trees[DITREES][0].getTaxonCount(); // never more than one diploid tree
		    }
			for (int i = 0; i < trees[TETRATREES].length; i++) {
				ntips += 2 * trees[TETRATREES][i].getTaxonCount();
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
			if (trees[DITREES].length > 0) {
				SimpleNode root = (SimpleNode) trees[DITREES][0].getRoot();
				simpletree2mullabtree(root, 0);
				rootn = nextn - 1;
				diroot = nodes[rootn];
			}
			if (diroot != null &&  diroot.height <= 0.0) {
				System.err.println("MulLabTree constructor: bug");
			}
				
			assert (diroot != null  ||  trees[TETRATREES].length == 1);
 
			// Tetratrees. Two copies each.
			// When there is one leg (JOINED case) make an extra node to join
			// the two tetratrees into one.
			// While doing this, make a list of LegLinks = roots and footUnion's
			// which specify the join to diploid branch.
			List<LegLink> leglinks = new ArrayList<LegLink>();
			for (int i = 0; i < trees[TETRATREES].length; i++) {
				SimpleNode root = (SimpleNode) trees[TETRATREES][i].getRoot();
				simpletree2mullabtree(root, 0);
				int r0 = nextn - 1;
				nodes[r0].hybridheight = trees[TETRATREES][i].getHybridHeight();
				nodes[r0].tetraroot = true;
				simpletree2mullabtree(root, 1);
				int r1 = nextn - 1;
				nodes[r1].hybridheight = trees[TETRATREES][i].getHybridHeight();
				nodes[r1].tetraroot = true;
				if (trees[TETRATREES][i].getNumberOfLegs() == 1) {
					// in this case, the two tettrees just made join at
					// another node (the root in no-diploid case, and a node
					// which goes on join to diploid tree otherwise).
					// grjtodo code shared here and in simpletree2mullabtree()
					nodes[nextn].child = new MulLabNode[2];
					nodes[nextn].child[0] = nodes[r0];
					nodes[r0].parent = nodes[nextn];
					nodes[nextn].child[1] = nodes[r1];
					nodes[r1].parent = nodes[nextn];
					nodes[nextn].height = trees[TETRATREES][i].getSplitHeight();
					if (diroot != null) {
						FixedBitSet u = trees[TETRATREES][i].getFootUnion(0);
						double footheight = trees[TETRATREES][i].getFootHeight(0);
						leglinks.add(new LegLink(nodes[nextn], nodeOfUnionInSubtree(diroot, u), footheight));
					}
					nextn++;
				} else {
					assert diroot != null;
					FixedBitSet u0 = trees[TETRATREES][i].getFootUnion(0);
					FixedBitSet u1 = trees[TETRATREES][i].getFootUnion(1);
					double footheight0 = trees[TETRATREES][i].getFootHeight(0);
					double footheight1 = trees[TETRATREES][i].getFootHeight(1);
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
		


		public double geneTreeInMULTreeLogLikelihood() {
			fillinpopvals();
			//System.out.println(asText());
			return geneTreeInMULSubtreeLogLikelihood(nodes[rootn]);
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
				nodes[nextn].name = snode.getTaxon().getId();
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
		 * to nodes in the MulLabTree. The population values are
		 * per-species-clade (per-branch in network), but of course more than 
		 * one node in MulLabTree may correspond to the same species.
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
		private void fillinpopvals() {
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
				
				if (trees[DITREES].length == 0) {
					p = fillinpopvalsforspunionNoDiploids(unionarray, n0, n1, p);
				} else if (trees[DITREES][0].getExternalNodeCount() == 2) {
					p = fillinpopvalsforspunionTwoDiploids(unionarray, n0, n1, p);
				} else {
					// grjtodo tetraonly Not sure if TwoDiploids code generalises (not sure of model either)
					assert false;					
				}
				n0 = n1;
			}
			assert p == popvalues.getDimension();
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
				pal = new PopulationAndLineages(t, tippop, node.rootpop(), node.nlineages);
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
		
		
		
	} // end MulLabTree
	
	
	

	/*
	 * this is for sorting the links from higher ploidy tree  to a single
	 * branch within a lower ploidy tree.
	 */
	static final Comparator<MulLabTree.LegLink> FOOTHEIGHT_ORDER = new Comparator<MulLabTree.LegLink>() {
		public int compare(MulLabTree.LegLink a, MulLabTree.LegLink b) {
			if (a.footheight == b.footheight) {
				return 0;
			} else {
				return (a.footheight > b.footheight) ? 1 : -1;
			}
		}
	};    
    
	
	
	
	/*
	 * This is for ordering the unions in the nodes of the MulLabTree.
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
	 * - corresponding to two roots of tetratrees in the MulLabTree plus a leg-join - 
	 * are ordered so that the two roots come first.
	 * 
	*/

	static final Comparator<MulLabTree.SpSqUnion> SPUNION_ORDER = new Comparator<MulLabTree.SpSqUnion>() {
		public int compare(MulLabTree.SpSqUnion a, MulLabTree.SpSqUnion b) {
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
    
	
	
		
	
/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 * 
 * 
 * 
 *                      AlloppSpeciesNetworkModel
 * 
 * +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 */
	
	
	
	/*
	 * Constructors. 
	 * 
	 */
	// This one only works for one tetra tree case. 
	// 
	public AlloppSpeciesNetworkModel(AlloppSpeciesBindings apspecies, double popvalue, boolean onehyb) {
		super(AlloppSpeciesNetworkModelParser.ALLOPPSPECIESNETWORK);
		apsp = apspecies;
		addModel(apsp);
		
		if (onehyb) {
			if (apsp.SpeciesWithinPloidyLevel(2).length == 0) {
				makeInitialOneTetraTreeNetwork();	
			} else if (apsp.SpeciesWithinPloidyLevel(2).length == 2) {
				makeInitialOneTetraTreeTwoDiploidsNetwork(LegType.JOINED);
			} else {
				// the restriction to two diploids is to avoid the need for topology changes in diptree
				assert false; // grjtodo tetraonly need to deal with other cases
			}				
		} else {
			assert false; // grjtodo morethanonetree need to deal with other cases
		}
		

		double maxrootheight = 0.0;
		for (int i = 0; i < trees[DITREES].length; i++) {
			double height = trees[DITREES][i].getMaxHeight();
			if (height > maxrootheight) { maxrootheight = height; }
		}
		for (int i = 0; i < trees[TETRATREES].length; i++) {
			double height = trees[TETRATREES][i].getMaxHeight();
			if (height > maxrootheight) { maxrootheight = height; }
		}
		double scale = 0.99 * apsp.initialMinGeneNodeHeight() / maxrootheight;
		scaleAllHeights(scale);
		
		int npopparams = numberOfPopParameters();
		this.popvalues = new Parameter.Default(npopparams, popvalue);
		addVariable(popvalues);
		
		mullabtree = new MulLabTree();
		
        Logger.getLogger("dr.evomodel.speciation.allopolyploid").info("\tConstructing an allopolyploid network,  please cite:\n"
                + Citable.Utils.getCitationString(this));
	}	
	
		
	
	/*
	 * This constructor is for testing.
	 */
	public AlloppSpeciesNetworkModel(AlloppSpeciesBindings apsp, 
			AlloppSpeciesNetworkModelTEST.NetworkToMulLabTreeTEST nmltTEST) {
		super(AlloppSpeciesNetworkModelParser.ALLOPPSPECIESNETWORK);
		this.apsp = apsp;
		popvalues = null;
	}
	
	
	/*
	 * This constructor is for testing.
	 */
	public AlloppSpeciesNetworkModel(AlloppSpeciesBindings testASB,
			AlloppSpeciesNetworkModelTEST.LogLhoodGTreeInNetworkTEST llgtnTEST) {
		super(AlloppSpeciesNetworkModelParser.ALLOPPSPECIESNETWORK);
		apsp = testASB;
		
		makeInitialOneTetraTreeNetwork(llgtnTEST);			
		
		this.popvalues = new Parameter.Default(llgtnTEST.popvalues);
		
		mullabtree = new MulLabTree();
	}

	
	public List<Citation> getCitations() {
		List<Citation> citations = new ArrayList<Citation>();
		citations.add(new Citation(
				new Author[]{
						new Author("GR", "Jones")
				},
				Citation.Status.IN_PREPARATION
		));
		return citations;
	}
	
	
	public String toString() {
		int ngt = apsp.numberOfGeneTrees();
		String nl = System.getProperty("line.separator");
		String s = nl + mullabtree.asText() + nl;
		for (int g = 0; g < ngt; g++) {
			s += "Gene tree " + g + nl;
			s += apsp.genetreeAsText(g) + nl;
			s += apsp.seqassignsAsText(g) + nl;
		}
		s += nl;
		return s;
	}
	
	
	
	public LogColumn[] getColumns() {
		LogColumn[] columns = new LogColumn[1];
		columns[0] = new LogColumn.Default("    MUL-tree and gene trees", this);
		return columns;
	}

	
	/*
	 * Stretches or squashes the whole network. Used by constructors and
	 * MCMC operators.
	 */
	public int scaleAllHeights(double scale) {
		int count = 0;
		for (int i = 0; i < trees[DITREES].length; i++) {
			count += trees[DITREES][i].scaleAllHeights(scale);
		}
		for (int i = 0; i < trees[TETRATREES].length; i++) {
			count += trees[TETRATREES][i].scaleAllHeights(scale);
		}
		return count;
	}


	/*
	 * Called from AlloppSpeciesBindings to check if a node in a gene tree
	 * is compatible with the network. 
	 */
	public boolean coalescenceIsCompatible(double height, FixedBitSet union) {
		MulLabNode node = mullabtree.nodeOfUnion(union);
		return (node.height <= height);
	}
	
	
    /*
     * Called from AlloppSpeciesBindings to remove coalescent information
     * from branches of mullabtree. Required before call to recordCoalescence
     */
	public void clearCoalescences() {
		mullabtree.clearCoalescences();
	}	

	
	/*
	 * Called from AlloppSpeciesBindings to add a node from a gene tree
	 * to its branch in mullabtree.
	 */
	public void recordCoalescence(double height, FixedBitSet union) {
		MulLabNode node = mullabtree.nodeOfUnion(union);
		assert (node.height <= height);
		while (node.parent != null  &&  node.parent.height <= height) {
			node = node.parent;
		}
		node.coalheights.add(height);
	}
	
	
	public void sortCoalescences() {
		for (MulLabNode node : mullabtree.nodes) {
			Collections.sort(node.coalheights);
		}
	}

	
	/*
	 * Records the number of gene lineages at nodes of mullabtree.
	 */
	public void recordLineageCounts() {
		mullabtree.recordLineageCounts();
	}	
	
	
	/*
	 * Calculates the log-likelihood for a single gene tree in the network
	 * 
	 * Requires that clearCoalescences(), recordCoalescence(), recordLineageCounts()
	 * called to fill mullabtree with information about gene tree coalescences first.
	 */
	public double geneTreeInNetworkLogLikelihood() {
		return mullabtree.geneTreeInMULTreeLogLikelihood();
	}	

	
	public int getTipCount() {
		int n = 0;
		for (int i = 0; i < trees[DITREES].length; i++) {
			n += trees[DITREES][i].getExternalNodeCount();
		}
		for (int i = 0; i < trees[TETRATREES].length; i++) {
			n += trees[TETRATREES][i].getExternalNodeCount();
		}
		return n;
	}
		
	/*
	 * Returns an array of heights for the events (speciations
	 * including those of extinct diploids in some cases). Used
	 * in calculation of proir for network.
	 */
	public double[] getEventHeights() {
		double heights[];
		if (trees[DITREES].length == 0  &&  trees[TETRATREES].length == 1) {
			AlloppLeggedTree astree = trees[TETRATREES][0];
			int n = astree.getExternalNodeCount();
			heights = new double[n];
			double[] inthgts = astree.getInternalHeights(); 
			assert inthgts.length == n-1;
			for (int i = 0; i < n-1; i++) {
				heights[i] = inthgts[i];
			}
			heights[n-1] = astree.getSplitHeight();
		} else if (trees[DITREES].length == 1  &&  trees[TETRATREES].length == 1) {
			double[] dipinthgts = trees[DITREES][0].getInternalHeights(); 
			double[] tetinthgts = trees[TETRATREES][0].getInternalHeights(); 
			int nlegs = trees[TETRATREES][0].getNumberOfLegs();
			assert nlegs > 0;
			double leghgt0 = trees[TETRATREES][0].getFootHeight(0);
			double leghgt1 = 0.0;
			if (nlegs == 2) {
				leghgt1 = trees[TETRATREES][0].getFootHeight(1);
			} else {
				leghgt1 = trees[TETRATREES][0].getSplitHeight();
			}
			assert leghgt0 > 0.0;
			assert leghgt1 > 0.0;
			heights = new double[dipinthgts.length + tetinthgts.length + 2];
			for (int i = 0; i < dipinthgts.length; i++) {
				heights[i] = dipinthgts[i];
			}
			for (int i = 0; i < tetinthgts.length; i++) {
				heights[dipinthgts.length + i] = tetinthgts[i];
			}			
			heights[dipinthgts.length + tetinthgts.length]	= leghgt0;		
			heights[dipinthgts.length + tetinthgts.length + 1]	= leghgt1;		
		} else {
			// grjtodo morethanonetree
			heights = new double[0];
		}
		return heights;
	}


	public int getNumberOfDiTrees() {
		return trees[DITREES].length;
	}
	
	public int getNumberOfTetraTrees() {
		return trees[TETRATREES].length;
	}

	public int getNumberOfNodeHeightsInTree(int pl, int n) {
		return trees[pl][n].getInternalNodeCount();
	}	
	
	
	public AlloppLeggedTree getHomoploidTree(int pl, int t) {
		return trees[pl][t];
	}
	
	
	
	public String mullabTreeAsText() {
		return mullabtree.asText();
	}

	

	// grjtodo beginNetworkEdit(), endNetworkEdit():
	// I am doing an analogous thing to what speciesTreeModel methods
	// beginTreeEdit(), endTreeEdit() do.
	// But I don't understand the purpose.
	// 2011-07-06, OK I am starting to: remake mullabtree() after edits.
	// Could do oldmullabtree = mullabtree, etc, instead.
	// Could use dirty flags instead, and remake on demand.
    public boolean beginNetworkEdit() {
    	boolean beingEdited = false;
    	for (int i=0; i<trees.length; i++) {
    		for (int j=0; j<trees[i].length; j++) {
    			beingEdited = beingEdited || trees[i][j].beginTreeEdit();
    		}
    	} 
    	return beingEdited;
    }


    public void endNetworkEdit() {
     	for (int i=0; i<trees.length; i++) {
    		for (int j=0; j<trees[i].length; j++) {
    			trees[i][j].endTreeEdit();
    		}
    	}
     	mullabtree = new MulLabTree();
     	fireModelChanged();
    }

	

	
	public String getName() {
		return getModelName();
	}

	
    //  based on SpeciesTreeModel
    //  grjtodo internalTreeOP remaining to do: scaling without enforcing consitency
 	public int scale(double scaleFactor, int nDims) throws OperatorFailedException {
    	assert scaleFactor > 0;
    	if (nDims <= 0) {
    		beginNetworkEdit();
    		int count = scaleAllHeights(scaleFactor);
    		endNetworkEdit();
    		fireModelChanged(this, 1);
    		return count;
    	} else {
    		if (nDims != 1) {
    			throw new OperatorFailedException("not implemented for count != 1");
    		}
    		/*
            if (internalTreeOP == null) {
                internalTreeOP = new TreeNodeSlide(this, species, 1);
            }

            internalTreeOP.operateOneNode(scaleFactor);*/
    		fireModelChanged(this, 1);
    		return nDims;
    	}	
    }


    // MCMC operator which moves the legs of a tetraploid subtree, 
    // ie, changes the way it joins the diploid tree
	public void moveLegs() {
		// grjtodo tetraonly morethanonetree
		// 2011-08-31 For now, the new legs are chosen from the same distribution regardless of current state.
		// I am not clear what this distribution should be. There are always two times,
		// both between the hybridizaton time and the diploid root. For the two-diploids case,
		// they are straightforward. Then there is the topology. There are five distinguishable
		// topologies in the two-diploids case: both legs left, both legs right, joined then left, 
		// joined then right, and one leg to each. I think the latter should be regarded as two cases
		// (most recent to left vs most recent to right) and it is here.
		assert trees[DITREES].length == 1;
		assert trees[DITREES][0].getExternalNodeCount() == 2;
		assert trees[TETRATREES].length == 1;
		if (MathUtils.nextInt(2) == 0) {
			// change times but not topology
			if (MathUtils.nextInt(2) == 0) {
				trees[TETRATREES][0].moveMostRecentLegHeight();
			} else {
				trees[TETRATREES][0].moveMostAncientLegHeight(trees[DITREES][0].getRootHeight());
			}
		} else {
			// change topology but not times 
			trees[TETRATREES][0].moveLegTopology(diploidtipbitset(0), diploidtipbitset(1));
		}
		
	}    

    
	protected void handleModelChangedEvent(Model model, Object object, int index) {
		if (DBUGTUNE)
			System.err.println("AlloppSpeciesNetworkModel.handleModelChangedEvent " + model.getId());
		fireModelChanged();
	}

	protected final void handleVariableChangedEvent(Variable variable,
							int index, Parameter.ChangeType type) {
		if (DBUGTUNE)
			System.err.println("AlloppSpeciesNetworkModel.handleVariableChangedEvent" + variable.getId());

	}

	
	
	
	protected void storeState() {
		oldtrees = new AlloppLeggedTree[NUMBEROFPLOIDYLEVELS][];
		for (int i=0; i<trees.length; i++) {
			oldtrees[i] = new AlloppLeggedTree[trees[i].length];
		}
    	for (int i=0; i<trees.length; i++) {
    		for (int j=0; j<trees[i].length; j++) {
    			oldtrees[i][j] = new AlloppLeggedTree(trees[i][j]);
    		}
    	}
    	// addVariable(popvalues) deals with popvalues
    	
		if (DBUGTUNE)
			System.err.println("AlloppSpeciesNetworkModel.storeState()");
	}

	
	
	protected void restoreState() {
		trees = new AlloppLeggedTree[NUMBEROFPLOIDYLEVELS][];
		for (int i=0; i<oldtrees.length; i++) {
			trees[i] = new AlloppLeggedTree[oldtrees[i].length];
		}
    	for (int i=0; i<oldtrees.length; i++) {
    		for (int j=0; j<oldtrees[i].length; j++) {
    			trees[i][j] = new AlloppLeggedTree(oldtrees[i][j]);
    		}
    	}		
    	mullabtree = new MulLabTree();
    	// addVariable(popvalues) deals with popvalues
    	
		if (DBUGTUNE)
			System.err.println("AlloppSpeciesNetworkModel.restoreState()");
	}

	
	protected void acceptState() {
	}
	
	
	// grjtodo not sure this is the best way to implement Units
    public Type getUnits() {
        return trees[DITREES][0].getUnits();
    }

    public void setUnits(Type units) {
        trees[DITREES][0].setUnits(units);
    }
	
	
	
	
	
	
	
    /*
     * For simple case of one tetraploid tree and no diploids.
     * Assume a history before root of a diploid speciating
     * at time s, the two diploids (or two descendants) forming
     * a hybrid at time h, which speciates at time r, the root
     * of the tetraploid tree.
     * 
     * heights is used for testing in cases of 2 or 3 tetraploids
     * (2011-06-13). Here heights[] supply r,h,s or d,r,h,s where
     * r,h,s are as above and d in the 3 tetraploid case is the
     * most recent tetraploid divergence.
     * 
     * The diploids go extinct or are not sampled.
     * 
     * There will be one population parameter at each tip,
     * one at the bottom of each branch within tetraploid tree
     * including the root node which goes down to h, and one
     * more population parameter for both diploids. The root
     * of the diploids will have population derived from this.
     */
	private void makeInitialOneTetraTreeNetwork() {
		
		trees = new AlloppLeggedTree[NUMBEROFPLOIDYLEVELS][];
		trees[DITREES] = new AlloppLeggedTree[0];
		
		Taxon[] spp = apsp.SpeciesWithinPloidyLevel(4);
		AlloppLeggedTree tettree = new AlloppLeggedTree(
				                       spp, LegType.NODIPLOIDS);
		trees[TETRATREES] = new AlloppLeggedTree[1];
		trees[TETRATREES][0] = tettree;
		
	}
	
	
	
	
	
    /*
     * For case of one tetraploid tree and two diploids.
     * Assume that a single diploid splits at the root,
     * and that a single hybridization event takes place between
     * the two initial diploids or descendants of them. Both initial
     * diploids leave exactly one descendant in the sample.
     * 
     * In other words, the data consist of two diploids and one
     * or more allotetraploids, and it is assumed that all
     * species arose from the MRCA of the two diploids, and that all
     * allotetraploids arose from a single hybridization event.
     * 
     * For a single tetraploid, these (plus reflections) are the 
     * possible evolutionary histories.
     * 
     * \     ||     /    \     ||     /   \     ||     /
     *  \....||..../      \....||..| /     \  |.||..| /
     *   \        /        \       |/       \ |     |/
     *    \      /          \      /         \|     /
     *     \    /            \    /           \    /
     *      \  /              \  /             \  /
     *       \/                \/               \/            
     *      
     * \  ||        /   \    ||      /
     *  \.||.|     /     \ |.||.|   /
     *   \   |    /       \|    |  /
     *    \  |   /         \   /  /
     *     \ |  /           \ |  /
     *      \| /             \| /
     *       \/               \/      
     *      
     * \    ||      /
     *  \ |.||.|   /
     *   \ \  /   /
     *    \ \/   /
     *     \/   /
     *      \  /
     *       \/           
     */
	private void makeInitialOneTetraTreeTwoDiploidsNetwork(LegType legtype) {
		
		trees = new AlloppLeggedTree[NUMBEROFPLOIDYLEVELS][];
		
		Taxon[] dipspp = apsp.SpeciesWithinPloidyLevel(2);		
		AlloppLeggedTree diptree = new AlloppLeggedTree(dipspp, LegType.NONE);
		trees[DITREES] = new AlloppLeggedTree[1];
		trees[DITREES][0] = diptree;
		
		Taxon[] tetspp = apsp.SpeciesWithinPloidyLevel(4);
		AlloppLeggedTree tettree = new AlloppLeggedTree(tetspp, legtype);
		// make the diploid root older than the earliest foot
		diptree.setNodeHeight(diptree.getRoot(), 1.1*tettree.getMaxHeight());
		FixedBitSet dip;
		switch (legtype) {
		case TWOBRANCH:
			// attach feet, one to each diploid
			for (int i = 0; i < 2; i++) {
				dip = diploidtipbitset(i);
				tettree.setFootUnion(i, dip);
			}
			break;
		case ONEBRANCH:
			// attach feet, both to one diploid
			dip = diploidtipbitset(0);
			for (int i = 0; i < 2; i++) {
				tettree.setFootUnion(i, dip);
			}
			break;
		case JOINED:
			// attach foot to a diploid
			dip = diploidtipbitset(0);
			tettree.setFootUnion(0, dip);
		}
		trees[TETRATREES] = new AlloppLeggedTree[1];
		trees[TETRATREES][0] = tettree;
	}
	
	
	private FixedBitSet diploidtipbitset(int i) {
		FixedBitSet dip = new FixedBitSet(apsp.numberOfSpSeqs());
		// note that ditree has been  constructed randomly, so getExternalNode(i) 
		// chooses an arbitrary node.
		String dipname = trees[DITREES][0].getNodeTaxon(trees[DITREES][0].getExternalNode(i)).getId();
		int sp = apsp.apspeciesId2index(dipname);
		int spseq = apsp.spandseq2spseqindex(sp, 0);
		// grjtodo tetraonly Don't like this way of finding the spseq index. 
		// I do similar code when constructing mullab tree, but that's OK
		// because seq is passed in. Here I am relying on diploids only 
		// having seq==0. It won't work for hexaploids joining tetras.
		// Maybe I want a apsp.firstspseqindexofspecies(species name)
		dip.set(spseq);
		return dip;
	}
	

	
	/*
	 * 2011-08-04. For no diploids, I am using one population parameter for
	 * of diploid history (1), one at tips of tetras (n), one at all branch bottoms
	 * in tettree (2n-2), and one for post-hybridization (1). Total 3n
	 * 
	 * For d>1 diploids, 3d-2 for ditree, and 3t+1 for a tetratree with t tips.
	 * 
	 * ditree is:  d (tips) + 2d-2 ( branch bottoms) = 3d-2
	 * 
	 * A tetratree is: t (tips) + 2t-2 ( branch bottoms) + 
	 *             2 (feet or foot+split) + 1 (hybridization) = 3t+1
	 * 
	 * Not sure what I want for other scenarios. 
	 */
	private int numberOfPopParameters() {
		int dim = 0;
		if (trees[DITREES].length == 0) {
			assert trees[TETRATREES].length == 1;
			int ntetratips = trees[TETRATREES][0].getExternalNodeCount();
			dim = 3*ntetratips; 
		} else {
			int nditips = trees[DITREES][0].getExternalNodeCount();
			if (nditips > 1) {
				dim += 3*nditips - 2;
				for (int i = 0; i < trees[TETRATREES].length; i++) {
					int ntetratips = trees[TETRATREES][i].getExternalNodeCount();
					dim += 3*ntetratips + 1;
				}
			} else {
				assert false; // grjtodo useful case??
			}
		}
		return dim;
	}


	
	private static int union2spseqindex(FixedBitSet union) {
		assert union.cardinality() == 1;
		return union.nextOnBit(0);
	}

	
	
/* *********************** TEST CODE **********************************/	
	
	
	
	/*
	 * Test of conversion from network to mullab tree
	 * 
	 * 2011-05-07 It is called from testAlloppSpeciesNetworkModel.java.
	 * I don't know how to put the code in there without
	 * making lots public here.
	 */
	// grjtodo. should be possible to pass stuff in nmltTEST. Currently
	// it just signals that this is indeed a test.
	public String testExampleNetworkToMulLabTree(
			AlloppSpeciesNetworkModelTEST.NetworkToMulLabTreeTEST nmltTEST) {
		AlloppLeggedTree[][] testtrees = new AlloppLeggedTree[NUMBEROFPLOIDYLEVELS][];
		
		Taxon[] spp = new Taxon[5];
		spp[0] = new Taxon("a");
		spp[1] = new Taxon("b");
		spp[2] = new Taxon("c");
		spp[3] = new Taxon("d");
		spp[4] = new Taxon("e");
		// 1,2,3, or b,c,d to be tets, 0 and 4 dips

		double tetheight = 0.0;
		// case 1. one tettree with one foot in each diploid branch
		// case 2. one tettree with both feet in one diploid branch
		// case 3. one tettree with one joined
		// case 4. two tettrees, 2+1, first with one foot in each diploid
		// branch, second joined
		// case 5. three tettrees, 1+1+1, one of each type of feet, as in cases 1-3

		int ntettrees = 0;
		switch (nmltTEST.testcase) {
		case 1:
		case 2:
		case 3:
			ntettrees = 1;
			break;
		case 4:
			ntettrees = 2;
			break;
		case 5:
			ntettrees = 3;
			break;
		}
		AlloppLeggedTree tettrees[] = new AlloppLeggedTree[ntettrees];

		Taxon[] tets123 = { spp[1], spp[2], spp[3] };
		Taxon[] tets12 = { spp[1], spp[2] };
		Taxon[] tets1 = { spp[1] };
		Taxon[] tets2 = { spp[2] };
		Taxon[] tets3 = { spp[3] };
		switch (nmltTEST.testcase) {
		case 1:
			tettrees[0] = new AlloppLeggedTree(tets123, nmltTEST, LegType.TWOBRANCH, 0.0);
			tetheight = tettrees[0].getMaxFootHeight();
			break;
		case 2:
			tettrees[0] = new AlloppLeggedTree(tets123, nmltTEST, LegType.ONEBRANCH, 0.0);
			tetheight = tettrees[0].getMaxFootHeight();
			break;
		case 3:
			tettrees[0] = new AlloppLeggedTree(tets123, nmltTEST, LegType.JOINED, 0.0);
			tetheight = tettrees[0].getMaxFootHeight();
			break;
		case 4:
			tettrees[0] = new AlloppLeggedTree(tets12, nmltTEST, LegType.TWOBRANCH, 0.0);
			tettrees[1] = new AlloppLeggedTree(tets3, nmltTEST, LegType.JOINED, 0.0);
			tetheight = tettrees[0].getMaxFootHeight();
			tetheight = Math.max(tetheight, tettrees[1].getMaxFootHeight());
			break;
		case 5:
			tettrees[0] = new AlloppLeggedTree(tets1, nmltTEST, LegType.TWOBRANCH, 0.0);
			tettrees[1] = new AlloppLeggedTree(tets2, nmltTEST, LegType.ONEBRANCH, 0.0);
			tettrees[2] = new AlloppLeggedTree(tets3, nmltTEST, LegType.JOINED, 0.0);
			tetheight = tettrees[0].getMaxFootHeight();
			tetheight = Math.max(tetheight, tettrees[1].getMaxFootHeight());
			tetheight = Math.max(tetheight, tettrees[2].getMaxFootHeight());
			break;
		}

		AlloppLeggedTree ditrees[] = new AlloppLeggedTree[1];
		Taxon[] dips = { spp[0], spp[4] };
		ditrees[0] = new AlloppLeggedTree(dips, nmltTEST, LegType.NONE, tetheight + 4.0);

		testtrees[0] = ditrees;
		testtrees[1] = tettrees;


		FixedBitSet a = new FixedBitSet(8);
		a.set(0);
		FixedBitSet e = new FixedBitSet(8);
		e.set(7);
		switch (nmltTEST.testcase) {
		case 1:
			// leg 0 to node a, leg to node e
			testtrees[TETRATREES][0].setFootUnion(0, a);
			testtrees[TETRATREES][0].setFootUnion(1, e);
			break;
		case 2:
			// both legs to node a,
			testtrees[TETRATREES][0].setFootUnion(0, a);
			testtrees[TETRATREES][0].setFootUnion(1, a);
			break;
		case 3:
			// only leg to node a
			testtrees[TETRATREES][0].setFootUnion(0, a);
			break;
		case 4:
			// first tet tree (with two tips): leg 0 to node a, leg to node e
			testtrees[TETRATREES][0].setFootUnion(0, a);
			testtrees[TETRATREES][0].setFootUnion(1, e);
			// second tet tree, only leg to node a
			testtrees[TETRATREES][1].setFootUnion(0, a);
			break;
		case 5:
			// first tet tree. leg 0 to node a, leg to node e
			testtrees[TETRATREES][0].setFootUnion(0, a);
			testtrees[TETRATREES][0].setFootUnion(1, e);
			// second tet tree. both legs to node a,
			testtrees[TETRATREES][1].setFootUnion(0, a);
			testtrees[TETRATREES][1].setFootUnion(1, a);
			// third tet tree. only leg to node a
			testtrees[TETRATREES][2].setFootUnion(0, a);
			break;
		}
		MulLabTree testmullabtree = new MulLabTree(testtrees);
		String newick = testmullabtree.mullabTreeAsNewick();
		return newick;
	}

	
	
	/*
	 * for testing
	 */
	private void makeInitialOneTetraTreeNetwork(AlloppSpeciesNetworkModelTEST.LogLhoodGTreeInNetworkTEST llgtnTEST) {
		
		trees = new AlloppLeggedTree[NUMBEROFPLOIDYLEVELS][];
		trees[DITREES] = new AlloppLeggedTree[0];

		Taxon[] spp = apsp.SpeciesWithinPloidyLevel(4);
		AlloppLeggedTree tettree = new AlloppLeggedTree(
				                       spp, LegType.NODIPLOIDS, llgtnTEST);
		trees[TETRATREES] = new AlloppLeggedTree[1];
		trees[TETRATREES][0] = tettree;
		
	}

	
	

}
