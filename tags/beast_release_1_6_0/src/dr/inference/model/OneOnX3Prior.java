package dr.inference.model;

import java.util.ArrayList;

/**
 * @author Chieh-Hsi Wu
 *
 * This an improper prior 1/(x_i^3) which is the Jeffrey's prior for parameters of Normal distribution
 * (and Log-normal distribution?)
 *
 */
public class OneOnX3Prior extends Likelihood.Abstract{

    public OneOnX3Prior() {

        super(null);
    }

    /**
     * Adds a statistic, this is the data for which the Prod_i (1/x_i^3) prior is calculated.
     *
     * @param data the statistic to compute density of
     */
    public void addData(Statistic data) {
        dataList.add(data);
    }


    protected ArrayList<Statistic> dataList = new ArrayList<Statistic>();

    /**
     * Overridden to always return false.
     */
    protected boolean getLikelihoodKnown() {
        return false;
    }

    /**
     * Calculate the log likelihood of the current state.
     *
     * @return the log likelihood.
     */
    public double calculateLogLikelihood() {

        double logL = 0.0;

        for (Statistic statistic : dataList) {
            for (int j = 0; j < statistic.getDimension(); j++) {
                logL -= 3*Math.log(statistic.getStatisticValue(j));
            }
        }
        return logL;
    }


    public String prettyName() {
        String s = "OneOnX3" + "(";
        for (Statistic statistic : dataList) {
            s = s + statistic.getStatisticName() + ",";
        }
        return s.substring(0, s.length() - 1) + ")";
    }
}
