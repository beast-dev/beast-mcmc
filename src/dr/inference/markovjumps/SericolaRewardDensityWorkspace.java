/*
 * SericolaRewardDensityWorkspace.java
 *
 * Reusable per-thread workspace for reward-density PDF and derivative calls.
 */

package dr.inference.markovjumps;

import java.util.Arrays;

final class SericolaRewardDensityWorkspace {

    private final int dim2;
    private final int phi;
    private final double epsilon;

    private final double[] scalarRewardProportion = new double[1];
    private final double[] scalarTime = new double[1];
    private final double[][] scalarOutput = new double[1][];

    private int[] H;
    private int[] NN;
    private double[] lt;
    private double[] premult;
    private double[] xh;
    private double[] oneMinus;
    private double[] ratio;
    private double[] w0;
    private boolean[] isZero;
    private boolean[] isOne;
    private double[] inc;
    private double[] contractionWeights;

    SericolaRewardDensityWorkspace(int dim, double epsilon) {
        this.dim2 = dim * dim;
        this.phi = dim - 1;
        this.epsilon = epsilon;
    }

    double[] scalarRewardProportions(double rewardProportion) {
        scalarRewardProportion[0] = rewardProportion;
        return scalarRewardProportion;
    }

    double[] scalarTimes(double time) {
        scalarTime[0] = time;
        return scalarTime;
    }

    double[][] scalarOutputs(double[] out) {
        scalarOutput[0] = out;
        return scalarOutput;
    }

    double preparePdf(
            double[] rewardProportions,
            double[] times,
            boolean parsimonious,
            double lambda,
            double[] sortedAlpha,
            double[] invAlphaDiff) {

        final int T = rewardProportions.length;
        final boolean singleTime = (times.length == 1);
        double maxT = 0.0;

        ensureCapacity(T);

        if (singleTime) {
            final double time = times[0];
            if (time <= 0.0) {
                throw new IllegalArgumentException("time must be > 0");
            }
            maxT = time;

            final double lambdaTime = lambda * time;
            final double premult0 = Math.exp(-lambdaTime);
            final int stepCount = (parsimonious ? determineNumberOfSteps(lambda, time) : Integer.MAX_VALUE);

            for (int t = 0; t < T; ++t) {
                lt[t] = lambdaTime;
                premult[t] = premult0;
                NN[t] = stepCount;

                H[t] = fillRewardProportionPrecomp(t, rewardProportions[t], sortedAlpha, invAlphaDiff);
            }
            return maxT;
        }

        for (int t = 0; t < T; ++t) {
            final double time = times[t];
            if (time <= 0.0) {
                throw new IllegalArgumentException("time must be > 0 at index " + t);
            }
            if (time > maxT) {
                maxT = time;
            }

            final double lambdaTime = lambda * time;
            lt[t] = lambdaTime;
            premult[t] = Math.exp(-lambdaTime);
            NN[t] = (parsimonious ? determineNumberOfSteps(lambda, time) : Integer.MAX_VALUE);

            H[t] = fillRewardProportionPrecomp(t, rewardProportions[t], sortedAlpha, invAlphaDiff);
        }
        return maxT;
    }

    int prepareDerivative(double rewardProportion, double[] sortedAlpha, double[] invAlphaDiff) {
        ensureCapacity(1);
        return fillRewardProportionPrecomp(0, rewardProportion, sortedAlpha, invAlphaDiff);
    }

    double prepareDerivative(
            double[] rewardProportions,
            double[] times,
            int count,
            double[] sortedAlpha,
            double[] invAlphaDiff) {

        if (count < 0 || count > rewardProportions.length) {
            throw new IllegalArgumentException("count must be in [0, rewardProportions.length]");
        }
        if (times.length != 1 && times.length < count) {
            throw new IllegalArgumentException("Either times.length==1 or times.length>=count");
        }

        ensureCapacity(count);

        final boolean singleTime = times.length == 1;
        double maxT = 0.0;

        if (singleTime) {
            final double time = times[0];
            if (time <= 0.0) {
                throw new IllegalArgumentException("time must be > 0");
            }
            maxT = time;

            for (int t = 0; t < count; ++t) {
                H[t] = fillRewardProportionPrecomp(t, rewardProportions[t], sortedAlpha, invAlphaDiff);
            }
            return maxT;
        }

        for (int t = 0; t < count; ++t) {
            final double time = times[t];
            if (time <= 0.0) {
                throw new IllegalArgumentException("time must be > 0 at index " + t);
            }
            if (time > maxT) {
                maxT = time;
            }

            H[t] = fillRewardProportionPrecomp(t, rewardProportions[t], sortedAlpha, invAlphaDiff);
        }
        return maxT;
    }

    int determineNumberOfSteps(double lambda, double time) {
        final double target = 1.0 - epsilon;
        final double lambdaTime = lambda * time;

        double term = Math.exp(-lambdaTime);
        double sum = term;

        int i = 0;
        final int hardCap = Math.max(5000, (int) (lambdaTime + 10.0 * Math.sqrt(lambdaTime + 1.0)));

        while (sum < target && i < hardCap) {
            i++;
            term *= lambdaTime / i;
            sum += term;
        }
        return i;
    }

    int[] intervals() {
        return H;
    }

    int[] stepCounts() {
        return NN;
    }

    double[] lambdaTimes() {
        return lt;
    }

    double[] premults() {
        return premult;
    }

    double[] ratios() {
        return ratio;
    }

    double[] bernsteinStartWeights() {
        return w0;
    }

    double[] oneMinus() {
        return oneMinus;
    }

    boolean[] isZero() {
        return isZero;
    }

    boolean[] isOne() {
        return isOne;
    }

    double[] increments() {
        return inc;
    }

    double[] contractionWeights() {
        return contractionWeights;
    }

    double xh(int index) {
        return xh[index];
    }

    boolean isZero(int index) {
        return isZero[index];
    }

    boolean isOne(int index) {
        return isOne[index];
    }

    private void ensureCapacity(int T) {
        if (H == null || H.length < T) H = new int[T];
        if (NN == null || NN.length < T) NN = new int[T];
        if (lt == null || lt.length < T) lt = new double[T];
        if (premult == null || premult.length < T) premult = new double[T];

        if (xh == null || xh.length < T) xh = new double[T];
        if (oneMinus == null || oneMinus.length < T) oneMinus = new double[T];
        if (ratio == null || ratio.length < T) ratio = new double[T];
        if (w0 == null || w0.length < T) w0 = new double[T];
        if (isZero == null || isZero.length < T) isZero = new boolean[T];
        if (isOne == null || isOne.length < T) isOne = new boolean[T];

        if (inc == null || inc.length != dim2) inc = new double[dim2];
        if (contractionWeights == null || contractionWeights.length != dim2) {
            contractionWeights = new double[dim2];
        }
    }

    private int fillRewardProportionPrecomp(
            int index,
            double rewardProportion,
            double[] sortedAlpha,
            double[] invAlphaDiff) {

        final int h = getHfromRewardProportion(rewardProportion, sortedAlpha);
        final double scaled = (rewardProportion - sortedAlpha[h - 1]) * invAlphaDiff[h];

        if (scaled <= 0.0) {
            xh[index] = 0.0;
            isZero[index] = true;
            isOne[index] = false;
            return h;
        }
        if (scaled >= 1.0) {
            xh[index] = 1.0;
            isZero[index] = false;
            isOne[index] = true;
            return h;
        }

        xh[index] = scaled;
        isZero[index] = false;
        isOne[index] = false;

        final double om = 1.0 - scaled;
        oneMinus[index] = om;
        ratio[index] = scaled / om;
        w0[index] = 1.0;

        return h;
    }

    private int getHfromRewardProportion(double rewardProportion, double[] sortedAlpha) {
        final double lo = sortedAlpha[0];
        final double hi = sortedAlpha[phi];

        if (rewardProportion < lo) {
            throw new IllegalArgumentException(
                    "rewardProportion=" + rewardProportion + " < min(alpha)=" + lo +
                            "; snapping to boundary is not supported; rewardProportion must be >= min(alpha)");
        }
        if (rewardProportion > hi) {
            throw new IllegalArgumentException(
                    "rewardProportion=" + rewardProportion + " > max(alpha)=" + hi +
                            "; snapping to boundary is not supported; rewardProportion must be <= max(alpha)");
        }

        final int index = Arrays.binarySearch(sortedAlpha, 0, phi + 1, rewardProportion);

        if (index >= 0) {
            if (index == 0) return 1;
            if (index == phi) return phi;
            return index + 1;
        }

        final int insertionPoint = -index - 1;
        return Math.max(1, Math.min(phi, insertionPoint));
    }
}
