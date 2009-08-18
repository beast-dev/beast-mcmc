package dr.inference.model;

/**
 * @author Marc Suchard
 */

public interface BayesianStochasticSearchVariableSelection {

    public Parameter getIndicators();

    public boolean validState();

    public class Utils {

        public static boolean connectedAndWellConditioned(double[] probability) {
            for(int i=0; i<probability.length; i++) {
                if(probability[i] == 0 || probability[i] > 1)
                    return false;
            }
            return true;
        }
    }
}
