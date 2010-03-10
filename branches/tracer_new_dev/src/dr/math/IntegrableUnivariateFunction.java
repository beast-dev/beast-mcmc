package dr.math;

/**
 * @author Marc A. Suchard
 */
public interface IntegrableUnivariateFunction extends UnivariateFunction {

    /* computer \int_{a}^{b} f(x) dx
     *
     * @param a lower bound of integration
     * @param b upper bound of integration
	 *
	 * @return function value
     */
    double evaluateIntegral(double a, double b);
}
