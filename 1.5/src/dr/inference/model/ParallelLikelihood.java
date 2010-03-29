package dr.inference.model;

/**
 * @author Marc A. Suchard
 */

public interface ParallelLikelihood extends Likelihood {

	public boolean getLikelihoodKnown();

	public void setLikelihood(double likelihood);

}
