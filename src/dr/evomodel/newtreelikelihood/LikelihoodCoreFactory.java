/**
 * 
 */
package dr.evomodel.newtreelikelihood;

import dr.evomodel.newtreelikelihood.GeneralLikelihoodCore;

/**
 * @author Marc Suchard
 *
 */
public class LikelihoodCoreFactory {
	
	public static LikelihoodCore getLikelihoodCore(int stateCount) {
				
		return new GeneralLikelihoodCore(stateCount);
	}

}
