package dr.inference.model;

/**
 * Created by IntelliJ IDEA.
 * User: msuchard
 * Date: Feb 3, 2007
 * Time: 2:09:26 PM
 * To change this template use File | Settings | File Templates.
 */
public interface ParallelLikelihood extends Likelihood {

	public boolean getLikelihoodKnown();

	public void setLikelihood(double likelihood);

}
