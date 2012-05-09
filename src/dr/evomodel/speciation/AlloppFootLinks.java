package dr.evomodel.speciation;

import java.util.List;

/**
 * Classes AlloppLegLink and AlloppFootLinks are for gathering and organising
 * the links between trees of different ploidy, so that the 
 * rootward-pointing legs can become tipward-pointing branches.
 * 
 * @author Graham Jones
 *         Date: 20/03/2012
 */


public class AlloppFootLinks {
	public List<AlloppLegLink> hips;
	public AlloppNode foot;

	public AlloppFootLinks(List<AlloppLegLink> hips) {
		this.hips = hips;
		if (hips.size() > 0) {
			foot = hips.get(0).foot;
			for (AlloppLegLink x : hips) {
				assert foot == x.foot;
			}
		}
	}
}
