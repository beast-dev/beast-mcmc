package dr.math.distributions;

/**
 * @author Guy Baele
 */
public class MultivariateKDEDistribution implements MultivariateDistribution {
	
	public static final String TYPE = "multivariateKDE";
    public static final boolean DEBUG = true;
	
	private Distribution[] multivariateKDE;
	private int dimension;
	//private boolean[] flags;

	public MultivariateKDEDistribution (Distribution[] multivariateKDE) {
		
		if (multivariateKDE.length <= 0) {
			throw new RuntimeException("Creation error in MultivariateKDEDistribution(Distribution[] multivariateKDE)");
		}
		
		this.multivariateKDE = multivariateKDE;
		this.dimension = multivariateKDE.length;
		/*for (int i = 0; i < dimension; i++) {
			flags[i] = true;
		}*/

	}
	
	public MultivariateKDEDistribution (Distribution[] multivariateKDE, boolean[] flags) {
		
		if (multivariateKDE.length <= 0) {
			throw new RuntimeException("Creation error in MultivariateKDEDistribution(Distribution[] multivariateKDE, boolean[] flags)");
		}
		
		this.multivariateKDE = multivariateKDE;
		this.dimension = multivariateKDE.length;
		//this.flags = flags;

	}

	public double logPdf(double[] x) {
		
		double logPdf = 0;
		
		if (x.length != dimension) {
            throw new IllegalArgumentException("data array is of the wrong dimension");
        }
		
		for (int i = 0; i < dimension; i++) {
			//if (flags[i]) {
			logPdf += multivariateKDE[i].logPdf(x[i]);
			//}
		}

        if (DEBUG){
            for (int i = 0; i < dimension; i++) {
                System.err.println(i + ", " + "x[i] = " + x[i] + ", logPdf = " + multivariateKDE[i].logPdf(x[i]));
                System.err.println("    mean = " + multivariateKDE[i].mean() + ", variance = " + multivariateKDE[i].variance());
            }
        }
		
		return logPdf;
	}

	public double[][] getScaleMatrix() {
		throw new RuntimeException("Not yet implemented");
	}

	public double[] getMean() {
		throw new RuntimeException("Not yet implemented");
	}

	public String getType() {
		return TYPE;
	}
	
}
