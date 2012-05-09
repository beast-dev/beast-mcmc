package dr.evomodel.operators;

import jebl.util.FixedBitSet;
import dr.evolution.tree.MutableTree;
import dr.evolution.tree.MutableTreeUtils;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.SlidableTree;
import dr.evolution.tree.Tree;
import dr.evomodel.speciation.AlloppDiploidHistory;
import dr.evomodel.speciation.AlloppLeggedTree;
import dr.evomodel.speciation.AlloppNode;
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
 * An operator for with allopolyploid networks. Uses mauCanonical() and mauReconstruct()
 * in SlidableTree.Utils to do some of the work. These are called to change one node in
 * diploid history or one node in a tetratree. There is also a move to change hyb height.
 * 
 * 
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
		public int ploidy;  // AlloppSpeciesNetworkModel.DITREES or TETRATREES
		public int tree;
		public int node;
		public boolean doHybheight;
		
		public NodeHeightInNetIndex(int ploidy, int tree, int node, boolean doHybheight) {
			this.ploidy = ploidy;
			this.tree = tree;
			this.node = node;
			this.doHybheight = doHybheight;
		}
	}
	
	
	private NodeHeightInNetIndex randomnode() {
		assert(1 == apspnet.getNumberOfDiTrees());
		int noftettrees = apspnet.getNumberOfTetraTrees();
		int dhcount;
		int count;
		
		dhcount = count = 0;
		// Diploid history has 2 extra tips for each tettree
		dhcount += apspnet.getNumberOfNodeHeightsInTree(AlloppSpeciesNetworkModel.DITREES, 0);
		dhcount += 2 * noftettrees;
		count = dhcount;
		// For each tetratree, the internal/root heights, plus hybrid time
		for (int i = 0; i < noftettrees; i++) {
			int n = apspnet.getNumberOfNodeHeightsInTree(AlloppSpeciesNetworkModel.TETRATREES, i);
			count += n+1; 
		}		
		int which = MathUtils.nextInt(count);
		if (which < dhcount) {
			return new NodeHeightInNetIndex(AlloppSpeciesNetworkModel.DITREES, 0, which, false);
		} else {
			which -= dhcount;
			for (int i = 0; i < noftettrees; i++) {
				int n = apspnet.getNumberOfNodeHeightsInTree(AlloppSpeciesNetworkModel.TETRATREES, i);
				if (which < n+1) {
					return new NodeHeightInNetIndex(AlloppSpeciesNetworkModel.TETRATREES, i, which, which==n);
				} else {
					which -= n+1;
				}
			}
		}
		return new NodeHeightInNetIndex(-1, -1, -1, false);
	}
	
	
	

	private void operateOneNodeInNet(double factor)
			throws OperatorFailedException {
		
		NodeHeightInNetIndex nhi = randomnode();
		if (nhi.ploidy == AlloppSpeciesNetworkModel.DITREES) {
			operateOneNodeInDiploidHistory(nhi.node, factor);
		} else {
			AlloppLeggedTree altree = apspnet.getHomoploidTree(nhi.ploidy, nhi.tree);
			if (nhi.doHybheight) {
				operateHybridHeightInLeggedTree(altree);
			} else {
			operateOneNodeInTetraTree(altree, nhi.node, factor);
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
	
	
	
	private void operateOneNodeInTetraTree(AlloppLeggedTree tree, int which, double factor) {
		
		// As TreeNodeSlide(). Randomly flip children at each node,
		// keeping track of node order (in-order order, left to right).

		NodeRef[] order = SlidableTree.Utils.mauCanonical(tree);
		
		// Find the time of the most recent gene coalescence which
		// has (species,sequence)'s to left and right of this node. 
		FixedBitSet left = apsp.speciesseqEmptyUnion();
		FixedBitSet right = apsp.speciesseqEmptyUnion();
		for (int k = 0; k < 2 * which + 1; k += 2) {
			FixedBitSet left0 = apsp.speciesseqToTipUnion(tree.getNodeTaxon(order[k]), 0);
			FixedBitSet left1 = apsp.speciesseqToTipUnion(tree.getNodeTaxon(order[k]), 1);
			left.union(left0);
			left.union(left1);
		}
		for (int k = 2 * (which + 1); k < order.length; k += 2) {
			FixedBitSet right0 = apsp.speciesseqToTipUnion(tree.getNodeTaxon(order[k]), 0);
			FixedBitSet right1 = apsp.speciesseqToTipUnion(tree.getNodeTaxon(order[k]), 1);
			right.union(right0);
			right.union(right1);
		}
		double genelimit = apsp.spseqUpperBound(left, right);
	
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
	    
	    apspnet.beginNetworkEdit();
		final NodeRef node = order[2 * which + 1];
		tree.setNodeHeight(node, newHeight);
		SlidableTree.Utils.mauReconstruct(tree, order);
		apspnet.endNetworkEdit();
	}
	
	
	
	
	private void operateOneNodeInDiploidHistory(int which, double factor) {
		int slidingn =  2 * which + 1;
		AlloppDiploidHistory diphist = new AlloppDiploidHistory(apspnet);

		NodeRef[] order = SlidableTree.Utils.mauCanonical(diphist);
		
		// Find the time of the most recent gene coalescence which
		// has (species,sequence)'s to left and right of this node. 
		FixedBitSet left = apsp.speciesseqEmptyUnion();
		FixedBitSet right = apsp.speciesseqEmptyUnion();
		for (int k = 0;  k < slidingn;  k += 2) {
			FixedBitSet u = diphist.getNodeUnion(order[k]);
			left.union(u);
		}
		for (int k = slidingn + 1;  k < order.length;  k += 2) {
			FixedBitSet u = diphist.getNodeUnion(order[k]);
			right.union(u);
		}
		double genelimit = apsp.spseqUpperBound(left, right);
		
		// find limit due to hyb-tips - adjacent heights must be bigger 
		double hybtiplimit = 0.0;
		if (slidingn-1 >= 0) {
			hybtiplimit = Math.max(hybtiplimit, diphist.getNodeHeight(order[slidingn-1]));
		}
		if (slidingn+1 < order.length) {
			hybtiplimit = Math.max(hybtiplimit, diphist.getNodeHeight(order[slidingn+1]));
		}	
		
		// find limit to keep root a diploid
		// 1. If node to slide is the root, and the second highest node is to left or 
		// right of all diploids, then the root must stay the root: lowerlimit =  second highest.
		// 2. If node to slide is not the root, and is to left or right of all diploids,
		// then it must not become the root: upperlimit = root height.
		double drootlowerlimit = 0.0;
		double drootupperlimit = Double.MAX_VALUE;
		int rootn = -1;
		double maxhgt = 0.0;
		for (int k = 1;  k < order.length;  k += 2) {
			double hgt = diphist.getNodeHeight(order[k]);
			if (hgt > maxhgt) {
				maxhgt = hgt;
				rootn = k;
			}
		}
		int secondn = -1;
		double secondhgt = 0.0;
		for (int k = 1;  k < order.length;  k += 2) {
			if (k != rootn) {
				double hgt = diphist.getNodeHeight(order[k]);
				if (hgt > secondhgt) {
					secondhgt = hgt;
					secondn = k;
				}				
			}
		}		
		int leftmostdip = -1;
		int rightmostdip = -1;
		for (int k = 0;  k < order.length;  k += 2) {
			if (diphist.nodeIsDiploidTip(order[k])) {
				if (leftmostdip < 0) {
					leftmostdip = k;
				}
				rightmostdip = k;
			}
		}		
		if (slidingn == rootn  &&  (secondn < leftmostdip  ||  secondn > rightmostdip)) {
			drootlowerlimit = diphist.getNodeHeight(order[secondn]);
		}
		if (slidingn < leftmostdip  ||  slidingn > rightmostdip) {
			drootupperlimit = diphist.getNodeHeight(order[rootn]);
		}

	    final double upperlimit = Math.min(genelimit, drootupperlimit);
	    final double lowerlimit = Math.max(hybtiplimit, drootlowerlimit);
	    
	    // On direct call, factor==0.0 and use limit. Else use passed in scaling factor
	    double newHeight = -1.0;
        if( factor > 0 ) {
            newHeight = diphist.getNodeHeight(order[slidingn]) * factor;
          } else {
            newHeight = MathUtils.uniform(lowerlimit, upperlimit);
          }
	    
        assert diphist.diphistOK();
		final NodeRef node = order[slidingn];
		diphist.setNodeHeight(node, newHeight);
		SlidableTree.Utils.mauReconstruct(diphist, order);
		if (!diphist.diphistOK()) {
			System.out.println("BUG in operateOneNodeInDiploidHistory()");
		}
        assert diphist.diphistOK();
		
	    apspnet.beginNetworkEdit();
		apspnet.replaceDiploidHistory(diphist); 
		apspnet.endNetworkEdit();		
	}
	
}
