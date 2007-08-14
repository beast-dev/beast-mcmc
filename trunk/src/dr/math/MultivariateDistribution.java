package dr.math;

import dr.inference.model.Parameter;

/**
 * Created by IntelliJ IDEA.
 * User: msuchard
 * Date: Jun 13, 2007
 * Time: 1:35:29 PM
 * To change this template use File | Settings | File Templates.
 */
public interface MultivariateDistribution {

    public double logPdf(Parameter x);

    public double[][] getScaleMatrix();

    public double[] getMean();

    public String getType();

}
