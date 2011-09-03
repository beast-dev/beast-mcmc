package dr.evomodel.operators;

import jebl.util.FixedBitSet;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.speciation.AlloppLeggedTree;
import dr.evomodel.speciation.AlloppSpeciesBindings;
import dr.evomodel.speciation.AlloppSpeciesNetworkModel;
import dr.evomodel.speciation.SpeciesTreeModel;
import dr.evolution.tree.SimpleTree;
import dr.evomodelxml.operators.AlloppNetworkNodeSlideParser;
import dr.inference.model.Parameter;
import dr.inference.operators.OperatorFailedException;
import dr.inference.operators.SimpleMCMCOperator;
import dr.math.MathUtils;


/**
 * 
 * @author Graham Jones
 *         Date: 01/07/2011
 */


/*
 * Extension of TreeNodeSlide to deal with allopolyploid networks.
 * TreeNodeSlide is an operator written by Joseph Heled for a species tree based on the ideas of Mau et all,  
 * <a href="http://citeseer.ist.psu.edu/rd/27056960%2C5592%2C1%2C0.25%2CDownload/http://citeseer.ist.psu.edu/cache/papers/cs/5768/ftp:zSzzSzftp.stat.wisc.eduzSzpubzSznewtonzSzlastfinal.pdf/mau98bayesian.pdf">
 * Bayesian Phylogenetic Inference via Markov Chain Monte Carlo Methods</a>
 */


public class AlloppNetworkNodeSlide extends SimpleMCMCOperator {

    private final AlloppSpeciesNetworkModel apspnet;
    private final AlloppSpeciesBindings apsp;
    
    public AlloppNetworkNodeSlide(AlloppSpeciesNetworkModel apspnet, AlloppSpeciesBindings apsp, double weight) {
		this.apspnet = apspnet;
		this.apsp = apsp;
		setWeight(weight);
	}

	public String getPerformanceSuggestion() {
		return "None";
	}

	@Override
	public String getOperatorName() {
		return AlloppNetworkNodeSlideParser.NETWORK_NODE_REHEIGHT + "(" + apspnet.getId() +
                      "," + apsp.getId() + ")";
		}

	@Override
	public double doOperation() throws OperatorFailedException {
		operateOneNodeInNet(0.0);
		return 0;
	}

	
	
	
	private class NodeHeightInNetIndex {
		public int pl;
		public int t;
		public int n;
		public int r;
		
		public NodeHeightInNetIndex(int pl, int t, int n, int r) {
			this.pl = pl;
			this.t = t;
			this.n = n;
			this.r = r;
		}
	}
	
	
	private NodeHeightInNetIndex randomnode() {
		int nofditrees = apspnet.getNumberOfDiTrees();
		int noftettrees = apspnet.getNumberOfTetraTrees();
		int count;
		if (nofditrees <= 1  &&  noftettrees == 1) {
			//grjtodo 2011-08-11. tetraonly morethanonetree. 
			// For now, don't change diploid tree with networkslide
			int nheights = apspnet.getNumberOfNodeHeightsInTree(AlloppSpeciesNetworkModel.TETRATREES, 0);
			count = nheights + 2; // internal/root heights, plus hybrid and split
			int which = MathUtils.nextInt(count);
			if (which < nheights) {
				return new NodeHeightInNetIndex(AlloppSpeciesNetworkModel.TETRATREES, 0, which, -1);
			} else {
				which -= nheights;
				return new NodeHeightInNetIndex(AlloppSpeciesNetworkModel.TETRATREES, 0, -1, which);
			}
		} else {
			assert nofditrees == 1;
			//grjtodo morethanonetree. This is untested, and will not be tested for some time...
			count = 0;
			int nheights = apspnet.getNumberOfNodeHeightsInTree(AlloppSpeciesNetworkModel.DITREES, 0);
			count += nheights-1;
			for (int i = 0; i < noftettrees; i++) {
				nheights = apspnet.getNumberOfNodeHeightsInTree(AlloppSpeciesNetworkModel.TETRATREES, i);
				count += nheights+3; // internal/root heights, plus hybrid plus two feet or split and foot
			}
			int which = MathUtils.nextInt(count);
			nheights = apspnet.getNumberOfNodeHeightsInTree(AlloppSpeciesNetworkModel.DITREES, 0);
			if (which < nheights) {
				return new NodeHeightInNetIndex(AlloppSpeciesNetworkModel.DITREES, 0, which, -1);
			} else {
				which -= nheights;
			}
			for (int i = 0; i < noftettrees; i++) {
				nheights = apspnet.getNumberOfNodeHeightsInTree(AlloppSpeciesNetworkModel.TETRATREES, i);
				if (which <nheights) {
					return new NodeHeightInNetIndex(AlloppSpeciesNetworkModel.TETRATREES, i, which, -1);
				} else {
					which -= nheights;
				}
				if (which < 3) {
					return new NodeHeightInNetIndex(AlloppSpeciesNetworkModel.TETRATREES, i, -1, which);
				} else {
					which -= 3;
				}
			}
		}
		assert false;
		return new NodeHeightInNetIndex(-1, -1, -1, -1);
	}
	

	private void operateOneNodeInNet(double factor)
			throws OperatorFailedException {
		
		NodeHeightInNetIndex nhi = randomnode();
		AlloppLeggedTree altree = apspnet.getHomoploidTree(nhi.pl, nhi.t);
		if (nhi.n >= 0) {
			// change node height nhi.n within tree trees[nhi.pl][nhi.t]
			operateOneNodeInLeggedTree(altree, nhi.n, factor);
		} else {
			// change hybrid (nhi.r==0) or split/foot heights (nhi.r=1,2) within tree trees[nhi.pl][nhi.t]
			if (nhi.r == 0) {
				operateHybridHeightInLeggedTree(altree);
			} else {
				operateLegsInLeggedTree(altree, nhi.r, factor);
			}
		}
	}
	
	
	// grjtodo morethanonetree
	private void operateHybridHeightInLeggedTree(AlloppLeggedTree tree) {
		double h = tree.getHybridHeight();	
		double r = 	tree.getRootHeight();
		double maxh = Double.MAX_VALUE;
		for (int i=0; i< tree.getNumberOfLegs(); i++) {
			double fi = tree.getFootHeight(i);
			if (fi >= 0.0) {
				maxh = Math.min(maxh, fi);
			}
		}
		double s = 	tree.getSplitHeight();
		if (s > 0.0) {
			maxh = Math.min(maxh, s);
		}
		assert maxh < Double.MAX_VALUE;
		double newh = h;
		newh = r + (maxh-r) * MathUtils.nextDouble();
		apspnet.beginNetworkEdit();
		tree.setHybridHeight(newh);	
		apspnet.endNetworkEdit();
	}
	
	
	// grjtodo morethanonetree
	private void operateLegsInLeggedTree(AlloppLeggedTree tree, int which, double factor) {
		if (tree.getNumberOfLegs() == 1) {
			double s = tree.getSplitHeight();
			double h = tree.getHybridHeight();
			double f = tree.getFootHeight(0);
			double news = s;
			if (f < 0.0) {
				if (factor > 0.0) {
					news = h + (s-h) * factor;
				} else {
					double limit = apsp.diploidSplitUpperBound();
					news = h + (limit-h) * MathUtils.nextDouble();
				}
			}
			apspnet.beginNetworkEdit();
			tree.setSplitHeight(news);
			apspnet.endNetworkEdit();
		}
	}
	
	
	private void operateOneNodeInLeggedTree(AlloppLeggedTree tree, int which, double factor) {
		
		// As TreeNodeSlide(). Randomly flip children at each node,
		// keeping track of node order (in-order order, left to right)
		// and what got swapped.
		final int count = tree.getExternalNodeCount();  
		NodeRef[] order = new NodeRef[2 * count - 1];
		boolean[] swapped = new boolean[count - 1];
		mauCanonical(tree, order, swapped);

		// Find the time of the most recent gene coalescence which
		// has species to left and right of this node. Note these are 
		// species, not species plus sequence, since inside tetra tree.
		FixedBitSet left = new FixedBitSet(count);
		FixedBitSet right = new FixedBitSet(count);
		for (int k = 0; k < 2 * which + 1; k += 2) {
			left.set(apsp.apspeciesId2index(tree.getNodeTaxon(order[k]).getId()));
		}
		for (int k = 2 * (which + 1); k < 2 * count; k += 2) {
			right.set(apsp.apspeciesId2index(tree.getNodeTaxon(order[k]).getId()));
		}
		double genelimit = apsp.speciationUpperBound(left, right);
		// also keep this node more recent than the hybridization event that led to this tree.
		double hybridheight = tree.getHybridHeight();
	    final double limit = Math.min(genelimit, hybridheight);
	    
	    // On direct call, factor==0.0 and use limit. Else use passed in scaling factor
	    double newHeight = -1.0;
        if( factor > 0 ) {
            newHeight = tree.getNodeHeight(order[2*which+1]) * factor;
          } else {
            newHeight = MathUtils.nextDouble() * limit;
          }
	    
	    // As TreeNodeSlide(). Change the height and reconstruct the tree.
	    // Note I don't swap popvalues. I assign popvalues to nodes
	    // via an ordering of species (fillinpopvals() in network). 
	    // I don't know how well that works.
	    apspnet.beginNetworkEdit();
		final NodeRef node = order[2 * which + 1];
		tree.setNodeHeight(node, newHeight);
		mauReconstruct(tree, order, swapped);
		apspnet.endNetworkEdit();
	}
	

	
	
	// rest is as TreeNodeSlide
	
	
    /**
     * Obtain an ordering of tree tips from randomly swaping the children order in internal nodes.
     *
     * @param tree    tree to create order from
     * @param order   Nodes in their random order (only odd indices are filled)
     * @param wasSwapped  true if internal node was swapped
     */
    static private void mauCanonical(AlloppLeggedTree tree, NodeRef[] order, boolean[] wasSwapped) {
        mauCanonicalSub(tree, tree.getRoot(), 0, order, wasSwapped);
    }

    
    static private int mauCanonicalSub(AlloppLeggedTree tree, NodeRef node,
    		                     int loc, NodeRef[] order, boolean[] wasSwapped) {
        if( tree.isExternal(node) ) {
            order[loc] = node;     assert (loc & 0x1) == 0;
            return loc + 1;
        }

        final boolean swap = MathUtils.nextBoolean();
        //wasSwaped[(loc-1)/2] = swap;
       
        int l = mauCanonicalSub(tree, tree.getChild(node, swap ? 1 : 0), loc, order, wasSwapped);

        order[l] = node;   assert (l & 0x1) == 1;
        wasSwapped[(l-1)/2] = swap;

        l = mauCanonicalSub(tree, tree.getChild(node, swap ? 0 : 1), l+1, order, wasSwapped);
        return l;
    }

    
    static private void mauReconstruct(AlloppLeggedTree tree, NodeRef[] order, boolean[] swapped) {
        final NodeRef root = mauReconstructSub(tree, 0, swapped.length, order, swapped);
        if( tree.getRoot() != root ) {
            tree.setRoot(root);
        }
    }
    

    static private NodeRef mauReconstructSub(AlloppLeggedTree tree, int from, int to, 
    		                     NodeRef[] order, boolean[] wasSwapped) {
        if( from == to ) {
            return order[2*from];
        }

        int rootIndex = -1;
        {
            double h = -1;

            for(int i = from; i < to; ++i) {
                final double v = tree.getNodeHeight(order[2 * i + 1]);
                if( h < v ) {
                    h = v;
                    rootIndex = i;
                }
            }
        }

        final NodeRef root = order[2 * rootIndex + 1];
        
        final NodeRef lchild = tree.getChild(root, 0);
        final NodeRef rchild = tree.getChild(root, 1);

        NodeRef lTargetChild = mauReconstructSub(tree, from, rootIndex, order, wasSwapped);
        NodeRef rTargetChild = mauReconstructSub(tree, rootIndex+1, to, order, wasSwapped);

        if( wasSwapped[rootIndex] ) {
            NodeRef z = lTargetChild;
            lTargetChild = rTargetChild;
            rTargetChild = z;
        }

        if( lchild != lTargetChild ) {
            tree.replaceChild(root, lchild, lTargetChild);
        }

        if( rchild != rTargetChild ) {
            tree.replaceChild(root, rchild, rTargetChild);
        }

        return root;
    }
	
}
