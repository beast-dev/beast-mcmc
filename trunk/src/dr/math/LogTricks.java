package dr.math;

/**
 * Created by IntelliJ IDEA.
 * User: msuchard
 * Date: Jan 18, 2007
 * Time: 9:00:50 AM
 * To change this template use File | Settings | File Templates.
 */
public class LogTricks {

    public static final double maxFloat = 3.40282347E+38;
    public static final double logLimit = -maxFloat / 100;
    public static final double logZero = -maxFloat;
    public static final double NATS = 40;

    public static double logSumNoCheck(double x, double y) {
        double temp = y - x;
        if (Math.abs(temp) > NATS)
            return (x > y) ? x : y;
        else
            return x + StrictMath.log1p(StrictMath.exp(temp));
    }

    public static double logSum(double x, double y) {
        double temp = y - x;
        if (temp > NATS || x < logLimit)
            return y;
        if (temp < -NATS || y < logLimit)
            return x;
        return x + StrictMath.log1p(StrictMath.exp(temp));
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
