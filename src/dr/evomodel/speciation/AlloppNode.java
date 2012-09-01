package dr.evomodel.speciation;



import java.util.Stack;

import dr.evolution.tree.SimpleNode;
import dr.evolution.util.Taxon;

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

    String asText(int indentlen);
	
	void fillinUnionsInSubtree(int unionsize);
	AlloppNode nodeOfUnionInSubtree(FixedBitSet x);


    public abstract class Abstract implements AlloppNode {
		 

		public void fillinUnionsInSubtree(int unionsize) {
			if (nofChildren() > 0) {
				getChild(0).fillinUnionsInSubtree(unionsize);
				getChild(1).fillinUnionsInSubtree(unionsize);
				FixedBitSet union = new FixedBitSet(unionsize);
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






        public static String subtreeAsText(AlloppNode node, String s, Stack<Integer> x, int depth, String b) {
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
            if (node.nofChildren() > 0) {
                x.push(depth);
                subs += subtreeAsText(node.getChild(0), "", x, depth+1, "-");
                x.pop();
                subs += subtreeAsText(node.getChild(1), "", x, depth+1, "`-");
            }
            return s + subs;
        }



        /* For PopsIOSpeciesTreeModel, to restore state after MCMC move
           * Recursively copies the topology from subtree rooted at node into
           * tree implemented as array nodes[]. Fills in the unions at the tips:
           * using piosb which converts species name into union.
           */
        static int simpletree2piotree(PopsIOSpeciesBindings piosb, AlloppNode[] nodes, int nextn,
                                         SimpleNode snode) {
            if (snode.isExternal()) {
                Taxon tx = snode.getTaxon();
                nodes[nextn].setTaxon(tx.getId());
                nodes[nextn].setUnion(piosb.tipUnionFromTaxon(tx));
            } else {
                nextn = simpletree2piotree(piosb, nodes, nextn, snode.getChild(0));
                int c0 = nextn - 1;
                nextn = simpletree2piotree(piosb, nodes, nextn, snode.getChild(1));
                int c1 = nextn - 1;
                nodes[nextn].addChildren(nodes[c0], nodes[c1]);
            }
            nodes[nextn].setAnc(null);
            nodes[nextn].setHeight(snode.getHeight());
            nextn++;
            return nextn;
        }
		
	}


}
