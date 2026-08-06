/*
 * NormalKDEDistribution.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
 *
 * This file is part of BEAST.
 * See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * BEAST is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 *  BEAST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 *
 */

package dr.math.distributions;

import dr.inference.model.GradientProvider;
import dr.math.ComplexArray;
import dr.math.FastFourierTransform;
import dr.stats.DiscreteStatistics;
import dr.util.HeapSort;

import java.util.Random;

/**
 * @author Marc A. Suchard
 */
public class NormalKDEDistribution extends KernelDensityEstimatorDistribution implements GradientProvider {

    public static final int MINIMUM_GRID_SIZE = 512;
    public static final double DEFAULT_CUT = 3.0;
    public static final boolean DEBUG = false;

    public NormalKDEDistribution(double[] sample) {
        super(sample, null, null, null);
        setup(DEFAULT_CUT, MINIMUM_GRID_SIZE);
    }

    public NormalKDEDistribution(double[] sample, Double lowerBound, Double upperBound, Double bandWidth) {
        super(sample, lowerBound, upperBound, bandWidth);
        setup(DEFAULT_CUT, MINIMUM_GRID_SIZE);
    }

    public NormalKDEDistribution(Double[] sample) {
        this(sample, null, null, null, DEFAULT_CUT, MINIMUM_GRID_SIZE);
    }

    public NormalKDEDistribution(Double[] sample, Double lowerBound, Double upperBound, Double bandWidth) {
        this(sample, lowerBound, upperBound, bandWidth, DEFAULT_CUT, MINIMUM_GRID_SIZE);
    }

    public NormalKDEDistribution(Double[] sample, Double lowerBound, Double upperBound, Double bandWidth,
                                 int n) {
        this(sample, lowerBound, upperBound, bandWidth, DEFAULT_CUT, n);
    }

    public NormalKDEDistribution(Double[] sample, Double lowerBound, Double upperBound, Double bandWidth,
                                 double cut, int n) {
        super(sample, lowerBound, upperBound, bandWidth);
        setup(cut, n);
    }

    private void setup(double cut, int n) {
        this.gridSize = Math.max(n, MINIMUM_GRID_SIZE);
        if (this.gridSize > MINIMUM_GRID_SIZE) {
            this.gridSize = (int) Math.pow(2, Math.ceil(Math.log(this.gridSize) / Math.log(2.0)));
        }
        this.cut = cut;

        setBounds();

        if (DEBUG) {
            System.out.println("lo = " + lo);
            System.out.println("up = " + up);
        }

        densityKnown = false;
    }

    void setBounds() {
        from = DiscreteStatistics.min(this.sample) - this.cut * this.bandWidth;
        to = DiscreteStatistics.max(this.sample) + this.cut * this.bandWidth;

        if (DEBUG) {
            System.out.println("min: " + DiscreteStatistics.min(this.sample));
            System.out.println("max: " + DiscreteStatistics.max(this.sample));
            System.out.println("bandWidth = " + this.bandWidth);
            System.out.println("cut = " + this.cut);
            System.out.println("from = " + from);
            System.out.println("to = " + to);
        }

        lo = from - 4.0 * this.bandWidth;
        up = to + 4.0 * this.bandWidth;
    }

    public double getFromPoint() {
        return from;
    }

    public double getToPoint() {
        return to;
    }

    /**
     * Returns a linear approximation evaluated at pt
     *
     * @param x    data (assumed sorted increasingly
     * @param y    data
     * @param pt   evaluation point
     * @param low  return value if pt < x
     * @param high return value if pt > x
     * @return evaluated coordinate
     */
    double linearApproximate(double[] x, double[] y, double pt, double low, double high) {

        int i = 0;
        int j = x.length - 1;

        if (pt < x[i]) {
            return low;
        }
        if (pt > x[j]) {
            return high;
        }

        // Bisection search
        while (i < j - 1) {
            int ij = (i + j) / 2;
            if (pt < x[ij]) {
                j = ij;
            } else {
                i = ij;
            }
        }

        if (pt == x[j]) {
            return y[j];
        }
        if (pt == x[i]) {
            return y[i];
        }
        return y[i] + (y[j] - y[i]) * ((pt - x[i]) / (x[j] - x[i]));
    }

    /**
     * Returns the gradient of a log linear approximation evaluated at pt
     *
     * @param x    data (assumed sorted increasingly
     * @param y    data
     * @param pt   evaluation point
     * @return evaluated coordinate
     */
    private double gradientLogLinearApproximate(double[] x, double[] y, double pt) {

        int i = 0;
        int j = x.length - 1;

        if (pt < x[i] || pt > x[j]) {
            return 0.0;
        }

        // Bisection search
        while (i < j - 1) {
            int ij = (i + j) / 2;
            if (pt < x[ij]) {
                j = ij;
            } else {
                i = ij;
            }
        }

        double slope =  (y[j] - y[i]) / (x[j] - x[i]);
        double eval = y[i] + (y[j] - y[i]) * ((pt - x[i]) / (x[j] - x[i]));
//        return slope / eval;
        return 1 / ((pt - x[i]) + y[i] / slope);
    }

    private double[] rescaleAndTrim(double[] x) {
        final int length = x.length / 2;
        final double scale = 1.0 / x.length;
        double[] out = new double[length];
        for (int i = 0; i < length; ++i) {
            out[i] = x[i] * scale;
            if (out[i] < 0) {
                out[i] = 0;
            }
        }
        return out;
    }

    private double[] massdist(double[] x,
//                              double[] xmass,
                              double xlow, double xhigh, int ny) {

        int nx = x.length;
        double[] y = new double[ny * 2];

        final int ixmin = 0;
        final int ixmax = ny - 2;
        final double xdelta = (xhigh - xlow) / (ny - 1);

        for (int i = 0; i < ny; ++i) {
            y[i] = 0.0;
        }

        final double xmi = 1.0 / nx;
        for (int i = 0; i < nx; ++i) {
            final double xpos = (x[i] - xlow) / xdelta;
            final int ix = (int) Math.floor(xpos);
            final double fx = xpos - ix;
//            final double xmi = xmass[i];

            if (ixmin <= ix && ix <= ixmax) {
                y[ix] += (1 - fx) * xmi;
                y[ix + 1] += fx * xmi;
            } else if (ix == -1) {
                y[0] += fx * xmi;
            } else if (ix == ixmax + 1) {
                y[ix] += (1 - fx) * xmi;
            }
        }
        return y;
    }

    /**
     * Override for different kernels
     *
     * @param ordinates the points in complex space
     * @param bandWidth predetermined bandwidth
     */
    private void fillKernelOrdinates(ComplexArray ordinates, double bandWidth) {
                final int length = ordinates.length;
        final double a = 1.0 / (Math.sqrt(2.0 * Math.PI) * bandWidth);
        final double precision = -0.5 / (bandWidth * bandWidth);
        for (int i = 0; i < length; i++) {
            final double x = ordinates.real[i];
            ordinates.real[i] = a * Math.exp(x * x * precision);
        }
    }

    protected void computeDensity() {
        makeOrdinates();
        makeXGrid();
        transformData();
        densityKnown = true;
    }

    void transformData() {
        ComplexArray Y = new ComplexArray(massdist(this.sample, lo, up, this.gridSize));
        FastFourierTransform.fft(Y, false);

        ComplexArray product = Y.product(kOrdinates);
        FastFourierTransform.fft(product, true);

        densityPoints = rescaleAndTrim(product.real);
    }

    void makeOrdinates() {

        final int length = 2 * gridSize;
        if (kOrdinates == null) {
            kOrdinates = new ComplexArray(new double[length]);
        }

        // Fill with grid values
        final double max = 2.0 * (up - lo);
        double value = 0;
        final double inc = max / (length - 1);
        for (int i = 0; i <= gridSize; i++) {
            kOrdinates.real[i] = value;
            value += inc;
        }
        for (int i = gridSize + 1; i < length; i++) {
            kOrdinates.real[i] = -kOrdinates.real[length - i];
        }
        fillKernelOrdinates(kOrdinates, bandWidth);

        FastFourierTransform.fft(kOrdinates, false);
        kOrdinates.conjugate();
    }

    void makeXGrid() {
        // Make x grid
        xPoints = new double[gridSize];
        double x = lo;
        double delta = (up - lo) / (gridSize - 1);
        if (DEBUG) {
            System.out.println("X");
        }
        for (int i = 0; i < gridSize; i++) {
            xPoints[i] = x;
            x += delta;
            if (DEBUG) {
                System.out.println(xPoints[i]);
            }
        }
    }

    @Override
    protected double evaluateKernel(double x) {
        if (!densityKnown) {
            computeDensity();
        }
        return linearApproximate(xPoints, densityPoints, x, 0.0, 0.0);
    }

    @Override
    protected void processBounds(Double lowerBound, Double upperBound) {
        if ((lowerBound != null && lowerBound != Double.NEGATIVE_INFINITY) ||
                (upperBound != null && upperBound != Double.POSITIVE_INFINITY)) {
            throw new RuntimeException("NormalKDEDistribution must be unbounded");
        }
    }

    @Override
    protected void setBandWidth(Double bandWidth) {
        if (bandWidth == null) {
            // Default bandwidth
            this.bandWidth = bandwidthNRD(sample);
        } else
            this.bandWidth = bandWidth;

        densityKnown = false;
    }

//   bandwidth.nrd =
//   function (x)
//   {
//       r <- quantile(x, c(0.25, 0.75))
//       h <- (r[2] - r[1])/1.34
//       4 * 1.06 * min(sqrt(var(x)), h) * length(x)^(-1/5)
//   }

    private double bandwidthNRD(double[] x) {
        if (indices == null) {
            indices = new int[x.length];
            HeapSort.sort(x, indices);
        }

        final double h =
                (DiscreteStatistics.quantile(0.75, x, indices) - DiscreteStatistics.quantile(0.25, x, indices)) / 1.34;
        return 1.06 *
                Math.min(Math.sqrt(DiscreteStatistics.variance(x)), h) *
                Math.pow(x.length, -0.2);
    }

    void resetIndices(boolean transformIncreasing) {
        if (!transformIncreasing) { // if the transform in not increasing, reset the indices.
            indices = null;
        }
    }

    @Override
    public int getDimension() {
        return 1;
    }

    @Override
    public double[] getGradientLogDensity(Object obj) {
        double[] x = GradientProvider.toDoubleArray(obj);
        double[] result = new double[x.length];
        if (!densityKnown) {
            computeDensity();
        }
        for (int i = 0; i < x.length; ++i) {
            result[i] = gradientLogLinearApproximate(xPoints, densityPoints, x[i]);
        }
        return result;
    }

    double getGradLogDensity(double x) {
        if (!densityKnown) {
            computeDensity();
        }
        return gradientLogLinearApproximate(xPoints, densityPoints, x);
        // NOTE: This is the gradient of the log of the linear approximation.
        // NOTE: This is NOT the log linear approximation of the gradient of the kernels.
        // NOTE: Tis is NOT the KDE estimation of the true gradient.
    }

    private ComplexArray kOrdinates;
    double[] xPoints;
    double[] densityPoints;

    private int[] indices;

    private int gridSize;
    private double cut;
    double from;
    double to;
    double lo;
    double up;

    boolean densityKnown = false;

    public static void main(String[] args) {

        long start = System.currentTimeMillis();

        Random random = new Random(1234);

        Double[] samples = new Double[10000000];
        for (int i = 0; i < samples.length; i++) {
//            samples[i] = random.nextDouble();
            samples[i] = (random.nextGaussian() * 0.01) + 0.5;
//            samples[i] = 0.3 + (random.nextDouble() * 0.1);
        }

//        NormalKDEDistribution nKDE = new NormalKDEDistribution(samples);
        NormalKDEDistribution nKDE = new NormalKDEDistribution(TEST_DATA, null, null, 0.2);

        for (int i = 0; i < 100; i++) {
            nKDE.evaluateKernel(random.nextDouble());
        }

        double x = 0.00;
        double sum = 0.0;
        for (int i = 0; i <= 100; i++) {
            double y = nKDE.pdf(x);
            double z = nKDE.cdf(x);
            sum += y;
            System.out.println("kde[" + i + "]: " + x + ", " + y + ", " + z);
            x += 0.01;
        }
        System.out.println("mean: " + sum / 100);

        // play with quantiles...
        double interval = 0.95;
        double lower = nKDE.quantile((1.0 - interval) / 2);
        double upper = nKDE.quantile(interval + ((1.0 - interval) / 2));
        System.out.println("95% KDE Intervals: " + lower + ", " + upper);


        long end = System.currentTimeMillis();
        System.out.println("Time: " + (end - start));

    }

    private static final double[] TEST_DATA = new double[] { 0.3876371659982018, 0.3748482373285085, 0.321027946, 0.30250701238936056, 0.2945585759737494, 0.3359061504997863, 0.29743346527967907, 0.2801877283473099, 0.32263928820992305, 0.26913786114965305, 0.2723149726753001, 0.25077808227631887, 0.2767243155191876, 0.26451857984501803, 0.2610800301055078, 0.26106833009314606, 0.25904783166515183, 0.2786566659439036, 0.29165805499036124, 0.3403150730362148, 0.34806964823497794, 0.31626808116537813, 0.3249533691878216, 0.32351203861563926, 0.27552865667668297, 0.28311694896998063, 0.2757842371094443, 0.30123804177753505, 0.29635370805887123, 0.25812477422367786, 0.27548179049839566, 0.2645852366674732, 0.23345032222772633, 0.2565311364379997, 0.24704628620641236, 0.28164029699965265, 0.29915160474417185, 0.23658270834945855, 0.28403872243569867, 0.29801894794223605, 0.3038119886402988, 0.3021081848708458, 0.38907178508489154, 0.31542208328657695, 0.33175334499593584, 0.3499272299148587, 0.31941149225873305, 0.2896649153886082, 0.29557542551319116, 0.2473864535302983, 0.28331157509828436, 0.28506909491368193, 0.31267276460136717, 0.31903309122906715, 0.2700177321677554, 0.31293883669962275, 0.27716281579431606, 0.24558904640880344, 0.25630480515428394, 0.2566234793537018, 0.27652605484334475, 0.26225152623092013, 0.23569991663888828, 0.25591405217243857, 0.29823322006524805, 0.2529453932891342, 0.2670655617937811, 0.2696758627861866, 0.246232879, 0.2732587303557807, 0.3563017785644923, 0.3014708887489509, 0.3345168473322027, 0.3229090051249468, 0.28000415943216045, 0.2681517081778991, 0.2357912582292073, 0.25245016284565747, 0.25162938865957085, 0.2715805980742842, 0.23219340558608073, 0.26343841006869917, 0.2793827427728628, 0.32196293496761835, 0.27181740411565514, 0.33406923916575076, 0.36709466308390837, 0.31685554003264593, 0.3113475197840096, 0.36341043648376975, 0.3891567014108614, 0.4015847879019321, 0.46896516983919534, 0.4251664995951691, 0.49715310139773083, 0.4136660159410094, 0.4154024141726014, 0.423620689, 0.35756650305291116, 0.4390910785990308, 0.42634383538672393, 0.438570735};
}