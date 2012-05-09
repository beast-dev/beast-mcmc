
package dr.evolution.tree;

import dr.evolution.util.Taxon;
import dr.math.MathUtils;

/**
 * Interface for a binary tree which can support a Mau-et-al MCMC move 
 * in which node heights `slide' and can change topology.
 * 
 * @author Graham Jones
 * Date 2012-03-30
 *
 */

public interface SlidableTree {

	NodeRef getRoot();
	
	void setRoot(NodeRef root);
	
	int getNodeCount();
	
	Taxon getNodeTaxon(NodeRef node);
	
	double getNodeHeight(NodeRef node);
	
	void setNodeHeight(NodeRef node, double height);
	
	boolean isExternal(NodeRef node);
	
	NodeRef getChild(NodeRef node, int j);
	
	void replaceChildren(NodeRef node, NodeRef lft, NodeRef rgt);
	
	
	
	
	public class Utils {
	
		/*
		 * This is an implementation of the ideas of
         * 
         * Bayesian Phylogenetic Inference via Markov Chain Monte Carlo Methods
         * Bob Mau, Michael A. Newton, Bret Larget
         * Biometrics, Vol. 55, No. 1 (Mar., 1999), pp. 1-12 
         * 
		 * similar to implementation by Joseph Heled in TreeNodeSlide.
		 * This version works on a SlidableTree rather than a SimpleTree and is 
		 * simpler than JH's version. It does not keep track of what was swapped.
		 * Instead I assign popvalues to nodes via an ordering of species (see 
		 * fillinpopvals() in network). 
		 */

		
		/*
		 * For Mau-etal MCMC moves on tree.
		 * Do this, then change a node height, then call mauReconstruct()
		 * 
		 */
	    static public NodeRef[] mauCanonical(SlidableTree tree) {
			final int count = tree.getNodeCount();  
			NodeRef[] order = new NodeRef[count];
	        mauCanonicalSub(tree, tree.getRoot(), 0, order);
	        return order;
	    }

	    
		/*
		 * For Mau-etal MCMC moves on tree
		 */
	    static public void mauReconstruct(SlidableTree tree, NodeRef[] order) {
	        final NodeRef root = mauReconstructSub(tree, 0, (order.length - 1) / 2, order);
	        tree.setRoot(root);
	    }
	    
	    
	/*
	 * ****************** private methods **********************    
	 */
	    
	    
	    static private int mauCanonicalSub(SlidableTree tree, NodeRef node,
	    		                     int nextloc, NodeRef[] order) {
	        if( tree.isExternal(node) ) {
	            order[nextloc] = node;     assert (nextloc & 0x1) == 0;
	            nextloc++;
	            return nextloc;
	        }
	        final boolean swap = MathUtils.nextBoolean();
	        nextloc = mauCanonicalSub(tree, tree.getChild(node, swap ? 1 : 0), nextloc, order);
	        order[nextloc] = node;   assert (nextloc & 0x1) == 1;
	        nextloc++;
	        nextloc = mauCanonicalSub(tree, tree.getChild(node, swap ? 0 : 1), nextloc, order);
	        return nextloc;
	    }

	    

	    static private NodeRef mauReconstructSub(SlidableTree tree, int from, int to, 
	    		                     NodeRef[] order) {
	        if (from == to) {
	            return order[2*from];
	        }
	        int rootIndex = highestNode(tree, order, from, to);
	        NodeRef root = order[2 * rootIndex + 1];
	        NodeRef lft = mauReconstructSub(tree, from, rootIndex, order);
	        NodeRef rgt = mauReconstructSub(tree, rootIndex+1, to, order);
	        tree.replaceChildren(root, lft, rgt);
	        return root;
	    }
		
		
		
		
		

	    static private int highestNode(SlidableTree tree, NodeRef[] order, int from, int to) {
	    	int rootIndex = -1;
	    	double maxh = -1.0;

	    	for (int i = from; i < to; ++i) {
	    		final double h = tree.getNodeHeight(order[2 * i + 1]);
	    		if (h > maxh) {
	    			maxh = h;
	    			rootIndex = i;
	    		}
	    	}
	    	return rootIndex;
	    }
	    
	    

	}
}
