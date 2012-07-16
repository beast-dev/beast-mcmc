package dr.evomodel.speciation;

import jebl.util.FixedBitSet;

/**
 * A leg for a AlloppLeggedTree.
 * 
 * @author Graham Jones
 *         Date: 20/03/2012
 */

public class AlloppTreeLeg {
	/*
	 *  footUnion is the node in a lower ploidy tree Y whose branch 
	 *  contains the foot leading to this tree. footUnion specifies
	 *  the clade (of species) in Y at the node.
	 *  
	 *  grjtodo tetraonly. With hexaploids, etc, may need to identify Y?
	 *  Or clade is enough?
	 *  
	 *  height is the time of the foot (which is during a branch in Y).
	 *
	 */
	public FixedBitSet footUnion;
	public double height;
	
	/**
	 * clone constructor
	 */
	public AlloppTreeLeg(AlloppTreeLeg leg) {
		this.height = leg.height;
		this.footUnion = new FixedBitSet(leg.footUnion);
	}
	
	// Partial constructor. The leg dangles, unattached
	public AlloppTreeLeg(double height) {
		this.footUnion = new FixedBitSet(0);
		this.height = height;
	}
}

