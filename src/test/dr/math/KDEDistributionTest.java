package test.dr.math;

import dr.inference.distribution.DistributionLikelihood;
import dr.inference.model.Parameter;
import dr.inference.trace.LogFileTraces;
import dr.inference.trace.TraceException;
import dr.math.MathUtils;
import dr.math.distributions.LogTransformedNormalKDEDistribution;
import dr.math.distributions.LogitTransformedNormalKDEDistribution;
import dr.math.distributions.NormalKDEDistribution;
import dr.math.distributions.TransformedNormalKDEDistribution;
import dr.math.matrixAlgebra.Vector;

import java.io.File;
import java.io.IOException;

import static dr.math.distributions.TransformedNormalKDEDistribution.getLogTransformedNormalKDEDistribution;
import static dr.math.distributions.TransformedNormalKDEDistribution.getLogitTransformedNormalKDEDistribution;

/**
 * @author Marc A. Suchard
 * @author Paul Bastide
 */
public class KDEDistributionTest extends MathTestCase {

    public void testNormalCDF() {

        MathUtils.setSeed(666);
        int length = 100000;

        Double[] values = new Double[length];
        for (int i = 0; i < length; ++i) {
            values[i] = MathUtils.nextGaussian();
        }

        // Normal KDE
        NormalKDEDistribution distribution = new NormalKDEDistribution(values);

        double cdf = distribution.cdf(0);
        System.out.println("cdf(0.00) = " + cdf);
        assertEquals(cdf, 0.5, cdfTolerance);

        cdf = distribution.cdf(1.96);
        System.out.println("cdf(1.96) = " + cdf);
        assertEquals(cdf, 0.975, cdfTolerance);

        double quantile = distribution.quantile(0.68);
        System.out.println("quantile(0.68) = " + quantile);
        assertEquals(quantile, 0.468, cdfTolerance);

        quantile = distribution.quantile(0.25);
        System.out.println("quantile(0.25) = " + quantile);
        assertEquals(quantile, -0.674, cdfTolerance);
    }

    private final static double cdfTolerance = 1E-2;

    public void testLogNormalCDF() {
        // Log Normal KDE on truncated normal
        MathUtils.setSeed(666);
        int length = 100000;

        Double[] values = new Double[length];
        for (int i = 0; i < length; ++i) {
            values[i] = MathUtils.nextGaussian();
            if (values[i] < 0.0) values[i] = -values[i];
        }

        // Testing CDF
        TransformedNormalKDEDistribution distribution = getLogTransformedNormalKDEDistribution(values);
        double cdf = distribution.cdf(1E-3);
        System.out.println("cdfLog(0.00) = " + cdf);
        assertEquals(cdf, 0.0007978844, cdfTolerance);

        cdf = distribution.cdf(1.96);
        System.out.println("cdfLog(1.96) = " + cdf);
        assertEquals(cdf, 0.9500042, cdfTolerance);

        // Old is NOT working
        LogTransformedNormalKDEDistribution distributionOld = new LogTransformedNormalKDEDistribution(values);
        double cdfOld = distributionOld.cdf(1E-3);
        System.out.println("cdfLogOld(0.00) = " + cdfOld);
        cdfOld = distributionOld.cdf(1.96);
        System.out.println("cdfLogOld(1.96) = " + cdfOld);

    }

    public void testLogitNormalCDF() {
        // Logit Normal KDE on logit normal
        MathUtils.setSeed(666);
        int length = 100000;

        Double[] values = new Double[length];
        for (int i = 0; i < length; i++) {
            values[i] = MathUtils.nextGaussian();
            values[i] = Math.exp(values[i]) / (1 + Math.exp(values[i]));
        }

        // Testing CDF
        TransformedNormalKDEDistribution distribution = getLogitTransformedNormalKDEDistribution(values, 1.0);
        double cdf = distribution.cdf(0.5);
        System.out.println("cdfLogit(0.00) = " + cdf);
        assertEquals(cdf, 0.5, cdfTolerance);

        cdf = distribution.cdf(0.75);
        System.out.println("cdfLogit(0.75) = " + cdf);
        assertEquals(cdf, 0.8640314, cdfTolerance); // `R pnorm(logit(0.75))`

        cdf = distribution.cdf(0.25);
        System.out.println("cdfLogit(0.25) = " + cdf);
        assertEquals(cdf, 0.1359686, cdfTolerance); // `R pnorm(logit(0.25))`

        // Old is NOT working
        LogitTransformedNormalKDEDistribution distributionOld = new LogitTransformedNormalKDEDistribution(values, 1.0);
        double cdfOld = distributionOld.cdf(0.5);
        System.out.println("cdfLogitOld(0.00) = " + cdfOld);
        cdfOld = distributionOld.cdf(0.75);
        System.out.println("cdfLogitOld(0.75) = " + cdfOld);
        cdfOld = distributionOld.cdf(0.25);
        System.out.println("cdfLogitOld(0.25) = " + cdfOld);

    }

    public void testNormalKDE() {

        int length = 100;
        Double[] sample = new Double[length];
        for (int i = 0; i < length; i++) {
            sample[i] = (double) (i + 1);
        }

//        final int gridSize = 256;

        NormalKDEDistribution kde = new NormalKDEDistribution(sample, null, null, null);
        System.out.println("bw = " + kde.getBandWidth());
        double tolerance = 1E-4;
        assertEquals(rBandWidth[0], kde.getBandWidth(), tolerance);


        final int gridSize = NormalKDEDistribution.MINIMUM_GRID_SIZE;
        double[] testPoints = new double[gridSize];
        double x = kde.getFromPoint();
        double delta = (kde.getToPoint() - kde.getFromPoint()) / (gridSize - 1);
        for (int i = 0; i < gridSize; ++i) {
            testPoints[i] = x;
            x += delta;
        }
        System.out.println("Eval @ " + new Vector(testPoints));

        double[] testDensity = new double[gridSize];
        for (int i = 0; i < gridSize; ++i) {
            testDensity[i] = kde.pdf(testPoints[i]);
        }
        System.out.println("Den    " + new Vector(testDensity));

        System.out.println("den[0] = " + testDensity[0]);
        System.out.println("den[N] = " + testDensity[NormalKDEDistribution.MINIMUM_GRID_SIZE - 1]);

        //  System.exit(-1);
    }

    public void testGammaKDE() {

//        int length = 10000;
//        double shape = 4;
//        double scale = 5;

//        double[] values = new double[length];
//        for (int i = 0; i < length; i++) {
//            values[i] = MathUtils.nextGamma(shape, scale);
//        }


//        GammaKDEDistribution kde = new GammaKDEDistribution(values);
//             System.out.println("prediction: at 2.02: "+kde.pdf(2.02177)+" at 0.405: "+kde.pdf(0.4046729)+" at 0.15: "+kde.pdf(0.1502078));
//     System.out.println("sm: "+kde.sampleMean());
        // TODO Need test values

    }

    private static double[] rBandWidth = { 12.24266 };

    /**
     * Test KDE construction from a BEAST log file
     * @param args 0: BEAST log file; 1: trace to construct KDE; 2: burn-in
     */
    public static void main(String[] args) {

        File file = new File(args[0]);
        String parent = file.getParent();
        file = new File(parent, args[0]);
        String fileName = file.getAbsolutePath();

        try {
            LogFileTraces traces = new LogFileTraces(fileName, file);
            traces.loadTraces();
            long maxState = traces.getMaxState();

            long burnin = Long.parseLong(args[2]);
            traces.setBurnIn(burnin);

            int traceIndexParameter = -1;
            for (int i = 0; i < traces.getTraceCount(); i++) {
                String traceName = traces.getTraceName(i);
                if (traceName.trim().equals(args[1])) {
                    traceIndexParameter = i;
                }
            }

            Double[] parameterSamples = new Double[traces.getStateCount()];
            DistributionLikelihood likelihood = new DistributionLikelihood(new LogTransformedNormalKDEDistribution((Double[]) traces.getValues(traceIndexParameter).toArray(parameterSamples)));
            //DistributionLikelihood likelihood = new DistributionLikelihood(new NormalKDEDistribution((Double[]) traces.getValues(traceIndexParameter).toArray(parameterSamples)));
            //DistributionLikelihood likelihood = new DistributionLikelihood(new GammaKDEDistribution((Double[]) traces.getValues(traceIndexParameter).toArray(parameterSamples)));

            Parameter parameter = new Parameter.Default(args[1]);
            parameter.setDimension(1);

            likelihood.addData(parameter);

            for (Double d : parameterSamples) {
                parameter.setParameterValue(0, d);
                System.out.println(d + " > " + likelihood.calculateLogLikelihood());
            }

        } catch (IOException ioe) {
            System.out.println(ioe);
        } catch (TraceException te) {
            System.out.println(te);
        }

    }

}

/*

# R test

sample = seq(1:100)
d = density(sample, bw="nrd")


*/