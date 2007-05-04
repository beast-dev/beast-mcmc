package dr.app.tracer.traces;

/**
 * A class that stores the correlation statistics for a trace
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: TraceCorrelation.java,v 1.2 2006/11/29 14:53:53 rambaut Exp $
 */
public class TraceCorrelation extends TraceDistribution {

    public TraceCorrelation(double[] values, int stepSize) {
        super(values, stepSize);

        if (isValid) {
            analyseCorrelation(values, stepSize);
        }
    }

    public double getStdErrorOfMean() { return stdErrorOfMean; }
    public double getACT() { return ACT; }

    private final static int MAX_OFFSET = 1000;

    /**
     * Analyze trace
     */
    private void analyseCorrelation(double[] values, int stepSize) {

        int maxLag = MAX_OFFSET;
        int samples = values.length;
        if ((samples/3) < maxLag) {
            maxLag = (samples/3);
        }

        double[] gammaStat = new double[maxLag];
        double[] varGammaStat = new double[maxLag];
        double varStat = 0.0;
        //double varVarStat = 0.0;
        //double assVarCor = 1.0;
        double del1, del2;

        for (int lag=0; lag < maxLag; lag++) {
            for (int j = 0; j < samples-lag; j++) {
                del1=values[j] - mean;
                del2=values[j + lag] - mean;
                gammaStat[lag] += ( del1*del2 );
                varGammaStat[lag] += (del1*del1*del2*del2);
            }

            gammaStat[lag] /= ((double)(samples-lag));
            varGammaStat[lag] /= ((double) samples-lag);
            varGammaStat[lag] -= (gammaStat[0] * gammaStat[0]);

            if (lag==0) {
                varStat = gammaStat[0];
                //varVarStat = varGammaStat[0];
                //assVarCor = 1.0;
            }
            else if (lag%2==0)
            {
                // fancy stopping criterion :)
                if (gammaStat[lag-1] + gammaStat[lag] > 0) {
                    varStat    += 2.0*(gammaStat[lag-1] + gammaStat[lag]);
                   // varVarStat += 2.0*(varGammaStat[lag-1] + varGammaStat[lag]);
                   // assVarCor  += 2.0*((gammaStat[lag-1] * gammaStat[lag-1]) + (gammaStat[lag] * gammaStat[lag])) / (gammaStat[0] * gammaStat[0]);
                }
                // stop
                else
                    maxLag=lag;
            }
        }

        // standard error of mean
        stdErrorOfMean = Math.sqrt(varStat/samples);

        // auto correlation time
        ACT = stepSize * varStat / gammaStat[0];

        // effective sample size
        ESS = (stepSize * samples) / ACT;

        // standard deviation of autocorrelation time
        stdErrOfACT = (2.0* Math.sqrt(2.0*(2.0*(double) (maxLag+1))/samples)*(varStat/gammaStat[0])*stepSize);

        isValid = true;
    }

    //************************************************************************
    // private methods
    //************************************************************************

    protected double stdErrorOfMean;
    protected double stdErrorOfVariance;
    protected double ACT;
    protected double stdErrOfACT;
}