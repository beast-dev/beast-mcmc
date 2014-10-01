package dr.inference.model;

import dr.stats.DiscreteStatistics;

/**
 * @author Simon Greenhill
 * @author Alexei Drummond
 */
public class PearsonCorrelation extends Statistic.Abstract {

    Parameter X, Y;
    boolean log = false;


    public PearsonCorrelation(Parameter X, Parameter Y, boolean log) {

        if (X.getDimension() != Y.getDimension()) throw new IllegalArgumentException();
        this.X = X;
        this.Y = Y;
        this.log = log;
    }

    public int getDimension() {
        return 1;
    }

    public double getStatisticValue(int dim) {

        double[] xvalues = X.getParameterValues();
        double[] yvalues = Y.getParameterValues();

        if (log) {
            for (int i = 0; i < xvalues.length; i++) {
                xvalues[i] = Math.log(xvalues[i]);
                yvalues[i] = Math.log(yvalues[i]);
            }
        }

        double meanX = DiscreteStatistics.mean(xvalues);
        double meanY = DiscreteStatistics.mean(yvalues);
        double stdevX = DiscreteStatistics.stdev(xvalues);
        double stdevY = DiscreteStatistics.stdev(yvalues);

        double corr = 0;
        for (int i = 0; i < xvalues.length; i++) {
            double deviateX = xvalues[i] - meanX;
            double deviateY = yvalues[i] - meanY;
            corr += deviateX*deviateY;
        }
        corr /= X.getDimension();

        corr /= (stdevX*stdevY);
        return corr;


    }
}
