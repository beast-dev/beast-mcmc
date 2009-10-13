package test.dr.distibutions;

import dr.math.distributions.GammaDistribution;
import dr.math.functionEval.GammaFunction;
import junit.framework.TestCase;

public class GammaDistributionTest extends TestCase{

	GammaDistribution gamma;
	
	public void setUp(){
		gamma = new GammaDistribution(1.0,1.0);
	}
	
	/**
	 * This test stochastically draws 100 gamma
	 * variates and compares the coded pdf 
	 * with the actual pdf.  
	 * The tolerance is required to be at most 1e-10.
	 */
	
	public void testPdf(){
		
		int numberOfTests = 100;
		
			for(int i = 0; i < numberOfTests; i++){
			double j = i/500.0 + 1.0;
			
			double shape = 1.0/j;
			double scale = j;
				
			gamma.setShape(shape);
			gamma.setScale(scale);
			
			double value = gamma.nextGamma();
			
			double mypdf = Math.pow(value, shape-1)/GammaFunction.gamma(shape)
					*Math.exp(-value/scale)/Math.pow(scale, shape);
			
			
			assertEquals(mypdf,gamma.pdf(value),1e-10);
		}
		
		gamma.setScale(1.0);
		gamma.setShape(1.0);

	}
	
	
	
	
}
