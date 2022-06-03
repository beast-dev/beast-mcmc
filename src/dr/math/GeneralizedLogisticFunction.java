package dr.math;

/**
 * @author Marc A. Suchard
 */
public class GeneralizedLogisticFunction {

    public static double evaluate(double x,
                                  double leftTime,
                                  double rightTime,
                                  double leftAsymptote,
                                  double rightAsymptote,
                                  double growthRate,
                                  double inflectionTime,
                                  double inflectionHeight) {
        return evaluate(x, leftTime, rightTime, leftAsymptote, rightAsymptote,
                growthRate, inflectionTime, inflectionHeight, 1.0, 1.0);
    }

    /**
     *  Based on https://en.wikipedia.org/wiki/Generalised_logistic_function
     */
    public static double evaluate(double x,
                                  double leftTime,
                                  double rightTime,
                                  double leftAsymptote, // A
                                  double rightAsymptote, // K
                                  double growthRate, // B
                                  double inflectionTime, // M
                                  double inflectionHeight, // Q
                                  double nu,
                                  double C) {

        final double value;
        if (x >= leftTime) {
            value =  leftAsymptote;
        } else if (x <= rightTime) {
            value = rightTime;
        } else {
            double timeOdds = (leftTime - x) / (x - rightTime);
            double inflectionOdds = (leftTime - inflectionTime) / (inflectionTime - rightTime);
            double inflation = Math.exp(-growthRate * (timeOdds - inflectionOdds));
            double weight = C + inflectionHeight * inflation;
            if (nu != 1.0) {
                weight = Math.pow(weight, 1.0 / nu);
            }
            value = leftAsymptote + (rightAsymptote - leftAsymptote) / weight;
        }
        return value;
    }
}
