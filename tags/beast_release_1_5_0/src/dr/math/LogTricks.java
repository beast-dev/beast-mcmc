package dr.math;

/**
 * @author Marc Suchard
 */
public class LogTricks {

    public static final double maxFloat = Double.MAX_VALUE; //3.40282347E+38;
    public static final double logLimit = -maxFloat / 100;
    public static final double logZero = -maxFloat;
    public static final double NATS =  400; //40;

    public static double logSumNoCheck(double x, double y) {
        double temp = y - x;
        if (Math.abs(temp) > NATS)
            return (x > y) ? x : y;
        else
            return x + StrictMath.log1p(StrictMath.exp(temp));
    }

    public static double logSum(double[] x) {
        double sum = x[0];
        final int len = x.length;
        for(int i=1; i<len; i++)
            sum = logSumNoCheck(sum,x[i]);
        return sum;
    }

    public static double logSum(double x, double y) {
        final double temp = y - x;
        if (temp > NATS || x < logLimit)
            return y;
        if (temp < -NATS || y < logLimit)
            return x;
        if (temp < 0)
            return x + StrictMath.log1p(StrictMath.exp(temp));
        return y + StrictMath.log1p(StrictMath.exp(-temp));
    }

    public static void logInc(Double x, double y) {
        double temp = y - x;
        if (temp > NATS || x < logLimit)
            x = y;
        else if (temp < -NATS || y < logLimit)
            ;
        else
            x += StrictMath.log1p(StrictMath.exp(temp));
    }

    public static double logDiff(double x, double y) {
        assert x > y;
        double temp = y - x;
        if (temp < -NATS || y < logLimit)
            return x;
        return x + StrictMath.log1p(-Math.exp(temp));
    }

}
