package dr.evomodel.speciation;

import java.util.Comparator;

/**
 * Classes AlloppLegLink and AlloppFootLinks are for gathering and organising
 * the links between trees of different ploidy, so that the 
 * rootward-pointing legs can become tipward-pointing branches.
 * 
 * @author Graham Jones
 *         Date: 20/03/2012
 */


public class AlloppLegLink {
	
	public AlloppNode foot;
	public AlloppNode hip;
	public double footheight;
	private boolean done;

	public AlloppLegLink(AlloppNode hip, AlloppNode foot, double footheight) {
		this.hip = hip;
		this.foot = foot;
		this.footheight =  footheight;
		done = false;
	}
	
	boolean isDone() {
		return done;
	}
	
	void setIsDone() {
		done = true;
	}
	
	/*
	 * this is for sorting the links from higher ploidy tree  to a single
	 * branch within diploid tree.
	 */
	static final Comparator<AlloppLegLink> FOOTHEIGHT_ORDER = new Comparator<AlloppLegLink>() {
		public int compare(AlloppLegLink a, AlloppLegLink b) {
			if (a.footheight == b.footheight) {
				return 0;
			} else {
				return (a.footheight > b.footheight) ? 1 : -1;
			}
		}
	};

}
