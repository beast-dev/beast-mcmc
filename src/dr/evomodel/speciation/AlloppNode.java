package dr.evomodel.speciation;


import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.SimpleNode;
import dr.evolution.util.Taxon;
import dr.math.MathUtils;

import jebl.util.FixedBitSet;


/**
 * An AlloppNode is an interface implemented by DipHistNode in AlloppDiploidHistory, 
 * and MulLabNode in AlloppMulLabTree.
 * The ABC AlloppNode.Abstract contains some common functionality.
 * 
 * @author Graham Jones
 *         Date: 20/03/2012
 */


public interface AlloppNode {
	int nofChildren();
	AlloppNode getChild(int ch);
	AlloppNode getAnc();
	Taxon getTaxon();
	double getHeight();
	FixedBitSet getUnion();
	
	void setChild(int ch, AlloppNode newchild);
	void setAnc(AlloppNode anc);
	void setTaxon(String name);
	void setHeight(double height);
	void setUnion(FixedBitSet union);
	void addChildren(AlloppNode c0, AlloppNode c1);

	
	void fillinUnionsInSubtree(AlloppSpeciesBindings apsp);
	AlloppNode nodeOfUnionInSubtree(FixedBitSet x);
	
	
	public abstract class Abstract implements AlloppNode {
		 

		public void fillinUnionsInSubtree(AlloppSpeciesBindings apsp) {
			if (nofChildren() > 0) {
				getChild(0).fillinUnionsInSubtree(apsp);
				getChild(1).fillinUnionsInSubtree(apsp);
				FixedBitSet union = apsp.speciesseqEmptyUnion();
				union.union(getChild(0).getUnion());
				union.union(getChild(1).getUnion());
				setUnion(union);
			}
		}	
		
		
		/* 
		 * For constructor of AlloppDiploidHistory, AlloppMulLabTree.
		 * Searches subtree rooted at node for most tipward node 
		 * whose union contains x. If x is known to be a union of one of the nodes,
		 * it finds that node, so acts as a map union -> node
		 */
		public AlloppNode nodeOfUnionInSubtree(FixedBitSet x) {
			if (nofChildren() == 0) {
				return this;
			}
			if (x.setInclusion(getChild(0).getUnion())) {
				return getChild(0).nodeOfUnionInSubtree(x);
			} else if (x.setInclusion(getChild(1).getUnion())) {
				return getChild(1).nodeOfUnionInSubtree(x);
			} else {
				return this;
			}
		}



		/*
		 * Specialized method which is used by AlloppMulLabTree, AlloppDiploidHistory constructors.
		 * 
		 * leglinks is a list of AlloppLegLink's which are (hip, foot, footheight)'s. These will
		 * be used later to insert the tetraploid subtrees into the diploid tree.
		 * 
		 * tetratree is the tetraploid subtree being dealt with
		 * 
		 * diroot is root node of diploid tree
		 * 
		 * nodes[] is the list of nodes representing the partially built tree, either a
		 * AlloppMulLabTree or a AlloppDiploidHistory.
		 *  
		 * nextn is the next free element in nodes[]
		 * 
		 * r0,r1 are the nodes in nodes[] which are the roots of the tetratree.
		 *  
		 * If tetratree has one leg (JOINED case) this makes an extra node to join
		 * the two tetratrees into one, and stores one AlloppLegLink. In the other
		 * cases, no new node made, and two AlloppLegLink are stored.
		 */
		static int collectLegLinksAndMakeLegJoin(List<AlloppLegLink> leglinks, 
				                   AlloppLeggedTree tetratree, AlloppNode diroot,
				                   AlloppNode nodes[], int nextn, int r0, int r1) {
			if (tetratree.getNumberOfLegs() == 1) {
				// in this case, the two hyb-tips just made join at
				// another node (the root in no-diploid case, and a node
				// which goes on join to diploid tree otherwise).
				nodes[nextn].addChildren(nodes[r0], nodes[r1]);
				nodes[nextn].setHeight(tetratree.getSplitHeight());
				if (diroot != null) {
					FixedBitSet u = tetratree.getFootUnion(0);
					double footheight = tetratree.getFootHeight(0);
					leglinks.add(new AlloppLegLink(nodes[nextn], diroot.nodeOfUnionInSubtree(u), footheight));
				}
				nextn++;
			} else {
				assert diroot != null;
				FixedBitSet u0 = tetratree.getFootUnion(0);
				FixedBitSet u1 = tetratree.getFootUnion(1);
				double footheight0 = tetratree.getFootHeight(0);
				double footheight1 = tetratree.getFootHeight(1);
				leglinks.add(new AlloppLegLink(nodes[r0], diroot.nodeOfUnionInSubtree(u0), footheight0));
				leglinks.add(new AlloppLegLink(nodes[r1], diroot.nodeOfUnionInSubtree(u1), footheight1));
			}
			return nextn;
		}
		
		

		/* 
		 * For constructor of AlloppDiploidHistory, AlloppMulLabTree.
		 */
		static void convertLegLinks(AlloppNode[] nodes, int nextn, List<AlloppLegLink> leglinks) {
			// Re-organise the LegLinks into per-diploid-branch lists (FootLinks)
			List<AlloppFootLinks> footlinkslist = new ArrayList<AlloppFootLinks>();
			for (AlloppLegLink x : leglinks) {
				if (!x.isDone()) {
					List<AlloppLegLink> hips = new ArrayList<AlloppLegLink>();
					for (AlloppLegLink y : leglinks) {
						if (y.foot == x.foot) {
							hips.add(y);
							y.setIsDone();
						}
					}
					footlinkslist.add(new AlloppFootLinks(hips));
				}
			}
			// sort the feet in time order.
			for (AlloppFootLinks x : footlinkslist) {
				Collections.sort(x.hips, AlloppLegLink.FOOTHEIGHT_ORDER);
			}

			// For all x, add nodes down from x.foot for all x.hips
			for (AlloppFootLinks x : footlinkslist) {
				AlloppNode f = x.foot;
				for (AlloppLegLink hip : x.hips) {
					AlloppNode h = hip.hip;
					nodes[nextn].setAnc(f.getAnc());
					nodes[nextn].addChildren(f, h);

					if (nodes[nextn].getAnc().getChild(0) == f) {
						nodes[nextn].getAnc().setChild(0, nodes[nextn]);
					} else {
						assert nodes[nextn].getAnc().getChild(1) == f;
						nodes[nextn].getAnc().setChild(1, nodes[nextn]);
					}
					nodes[nextn].setHeight(hip.footheight);
					f = nodes[nextn]; 
					nextn++;
				}
			}
		}
		
		
		
		/* For constructor of AlloppMulLabTree and AlloppDiploidHistory.
		 * Recursively copies the topology from subtree rooted at node into
		 * tree implemented as array nodes[]. Fills in the unions at the tips: 
		 * using apsp which converts species name and sequence index into union.
		 */
		static int simpletree2allopptree(AlloppSpeciesBindings apsp, AlloppNode[] nodes, int nextn, 
				SimpleNode snode, int seq) {
			if (snode.isExternal()) {
				nodes[nextn].setTaxon(snode.getTaxon().getId() + seq);
				nodes[nextn].setUnion(apsp.speciesseqToTipUnion(snode.getTaxon(), seq));
			} else {
				nextn = simpletree2allopptree(apsp, nodes, nextn, snode.getChild(0), seq);
				int c0 = nextn - 1;
				nextn = simpletree2allopptree(apsp, nodes, nextn, snode.getChild(1), seq);
				int c1 = nextn - 1;
				nodes[nextn].addChildren(nodes[c0], nodes[c1]);
			}
			nodes[nextn].setHeight(snode.getHeight());
			nextn++;
			return nextn;
		}
		
		
		/* For constructor of AlloppDiploidHistory. Like simpletree2allopptree but
		 * assumes diploids, so seq==0. Does NOT add seq to taxon names.
		 */
		static int simpletree2ditree(AlloppSpeciesBindings apsp, AlloppNode[] nodes, int nextn, 
				SimpleNode snode) {
			if (snode.isExternal()) {
				nodes[nextn].setTaxon(snode.getTaxon().getId());
				nodes[nextn].setUnion(apsp.speciesseqToTipUnion(snode.getTaxon(), 0));
			} else {
				nextn = simpletree2ditree(apsp, nodes, nextn, snode.getChild(0));
				int c0 = nextn - 1;
				nextn = simpletree2ditree(apsp, nodes, nextn, snode.getChild(1));
				int c1 = nextn - 1;
				nodes[nextn].addChildren(nodes[c0], nodes[c1]);
			}
			nodes[nextn].setHeight(snode.getHeight());
			nextn++;
			return nextn;
		}
		
	}


}
