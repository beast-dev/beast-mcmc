package dr.inference.trace;

import dr.math.LogTricks;
import dr.math.MathUtils;
import dr.util.Attribute;
import dr.xml.*;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.Reader;

/**
 * Created by IntelliJ IDEA.
 * User: msuchard
 * Date: Jan 17, 2007
 * Time: 6:45:12 PM
 * To change this template use File | Settings | File Templates.
 * <p/>
 * Source translated from model_P.c (a component of BAli-Phy by Benjamin Redelings and Marc Suchard
 */
public class    MarginalLikelihoodAnalysis {

    public static final String ML_ANALYSIS = "marginalLikelihoodAnalysis";
    public static final String FILE_NAME = "fileName";
    public static final String BURN_IN = "burnIn";
    public static final String COLUMN_NAME = "likelihoodColumn";
    public static final String DO_BOOTSTRAP = "bootstrap";
    public static final String ONLY_HARMONIC = "harmonicOnly";
    public static final String BOOTSTRAP_LENGTH = "bootstrapLength";

    private Trace trace;
    private int burnin;
    private boolean harmonicOnly;
    private boolean doBootstrap;
    private int bootstrapLength;


    private boolean marginalLikelihoodCalculated = false;
    private double logMarginalLikelihood;
    private double[] bootstrappedLogML;
    private double bootstrappedAverage;
    private double bootstrappedSE;


    public MarginalLikelihoodAnalysis(Trace trace, int burnin) {
        this.trace = trace;
        this.burnin = burnin;

        this.harmonicOnly = false;
        this.doBootstrap = true;
        this.bootstrapLength = 1000;
    }

    public MarginalLikelihoodAnalysis(Trace trace, int burnin, boolean harmonicOnly, boolean doBootstrap, int bootstrapLength) {
        this.trace = trace;
        this.burnin = burnin;
        this.harmonicOnly = harmonicOnly;
        this.doBootstrap = doBootstrap;
        this.bootstrapLength = bootstrapLength;
//        System.err.println("setting burnin to "+burnin);
    }

    public double calculateLogMarginalLikelihood(double[] sample) {
        if (harmonicOnly)
            return logMarginalLikelihoodHarmonic(sample);
        else
            return logMarginalLikelihoodSmoothed(sample);
    }

    /**
     * Calculates the log marginal likelihood of a model using Newton and Raftery's harmonic mean estimator
     *
     * @param v a posterior sample of logLikelihoods
     * @return the log marginal likelihood
     */

    public double logMarginalLikelihoodHarmonic(double[] v) {

        double sum = 0;
        final int size = v.length;
        for (int i = 0; i < size; i++)
            sum += v[i];

        double denominator = LogTricks.logZero;

        for (int i = 0; i < size; i++)
            denominator = LogTricks.logSum(denominator, sum - v[i]);

        return sum - denominator + StrictMath.log(size);
    }

    public void update() {
        double sample[] = trace.getValues(burnin);

//        System.err.println("using "+sample.length+" samples.") ;
//        if (!marginalLikelihoodCalculated)      {
        logMarginalLikelihood = calculateLogMarginalLikelihood(sample);
//        }
        if (doBootstrap) {
            final int bsLength = bootstrapLength;
            final int sampleLength = sample.length;
            double[] bsSample = new double[sampleLength];
            bootstrappedLogML = new double[bsLength];
            double sum = 0;
            for (int i = 0; i < bsLength; i++) {
                int[] indices = MathUtils.sampleIndicesWithReplacement(sampleLength);
                for (int k = 0; k < sampleLength; k++)
                    bsSample[k] = sample[indices[k]];
                bootstrappedLogML[i] = calculateLogMarginalLikelihood(bsSample);
                sum += bootstrappedLogML[i];
            }
            sum /= bsLength;
            bootstrappedAverage = sum;
            // Summarize bootstrappedLogML
            double var = 0;
            for (int i = 0; i < bsLength; i++) {
                var += (bootstrappedLogML[i] - sum) * (bootstrappedLogML[i] - sum);
            }
            var /= (bsLength - 1.0);
            bootstrappedSE = Math.sqrt(var);
        }
        marginalLikelihoodCalculated = true;
    }

    /**
     * Calculates the log marginal likelihood of a model using Newton and Raftery's smoothed estimator
     *
     * @param v     a posterior sample of logLilelihood
     * @param delta proportion of pseudo-samples from the prior
     * @param Pdata current estimate of the log marginal likelihood
     * @return the log marginal likelihood
     */
    public double logMarginalLikelihoodSmoothed(double[] v, double delta, double Pdata) {

        final double logDelta = StrictMath.log(delta);
        final double logInvDelta = StrictMath.log(1.0 - delta);
        final int n = v.length;
        final double logN = StrictMath.log(n);

        final double offset = logInvDelta - Pdata;

        double bottom = logN + logDelta - logInvDelta;
        double top = bottom + Pdata;

        for (int i = 0; i < n; i++) {
            double weight = -LogTricks.logSum(logDelta, offset + v[i]);
            top = LogTricks.logSum(top, weight + v[i]);
            bottom = LogTricks.logSum(bottom, weight);
        }

        return top - bottom;
    }

    public String toString() {
        if (!marginalLikelihoodCalculated)
            update();
        StringBuilder sb = new StringBuilder();
        sb.append("log P(").append(trace.getId()).append("|Data) = ").append(String.format("%5.4f", logMarginalLikelihood));
        if (doBootstrap) {
            sb.append(" +/- ").append(String.format("%5.4f", bootstrappedSE));
        } else {
            sb.append("           ");
        }
        if (harmonicOnly)
            sb.append(" (harmonic)");
        else
            sb.append(" (smoothed)");
        sb.append(" burnin=").append(burnin);
        if (doBootstrap)
            sb.append(" replicates=").append(bootstrapLength);
//        sb.append("\n");

        return sb.toString();

    }

    public double logMarginalLikelihoodSmoothed(double[] v) {

        final double delta = 0.01;  // todo make class adjustable by accessor/setter

        // Start with harmonic estimator as first guess
        double Pdata = logMarginalLikelihoodHarmonic(v);

        double deltaP = 1.0;

        int iterations = 0;

        double dx = 10.0;

        final double tolerance = 1E-3; // todo make class adjustable by accessor/setter

        while (Math.abs(deltaP) > tolerance) {
            double g1 = logMarginalLikelihoodSmoothed(v, delta, Pdata) - Pdata;
            double Pdata2 = Pdata + g1;
            dx = g1 * 10.0;
            double g2 = logMarginalLikelihoodSmoothed(v, delta, Pdata + dx) - (Pdata + dx);
            double dgdx = (g2 - g1) / dx; // find derivative at Pdata

            double Pdata3 = Pdata - g1 / dgdx; // find new evaluation point
            if (Pdata3 < 2.0 * Pdata || Pdata3 > 0 || Pdata3 > 0.5 * Pdata) // step is too large
                Pdata3 = Pdata + 10.0 * g1;

            double g3 = logMarginalLikelihoodSmoothed(v, delta, Pdata3) - Pdata3;

            // Try to do a Newton's method step
            if (Math.abs(g3) <= Math.abs(g2) && ((g3 > 0) || (Math.abs(dgdx) > 0.01))) {
                deltaP = Pdata3 - Pdata;
                Pdata = Pdata3;
            }  // otherwise try to go 10 times as far as one step
            else if (Math.abs(g2) <= Math.abs(g1)) {
                Pdata2 += g2;
                deltaP = Pdata2 - Pdata;
                Pdata = Pdata2;
            }  // otherwise go just one step
            else {
                deltaP = g1;
                Pdata += g1;
            }

            iterations++;

            if (iterations > 400) { // todo make class adjustable by acessor/setter
                System.err.println("Probabilities are not converging!!!"); // todo should throw exception
                return LogTricks.logZero;
            }
        }
        return Pdata;
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return ML_ANALYSIS;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            try {
                Reader reader;

                String fileName = xo.getStringAttribute(FILE_NAME);
                try {
                    File file = new File(fileName);
                    String name = file.getName();
                    String parent = file.getParent();

                    if (!file.isAbsolute()) {
                        parent = System.getProperty("user.dir");
                    }

                    reader = new FileReader(new File(parent, name));
                } catch (FileNotFoundException fnfe) {
                    throw new XMLParseException("File '" + fileName + "' can not be opened for " + getParserName() + " element.");
                }

                XMLObject cxo = (XMLObject) xo.getChild(COLUMN_NAME);
                String likelihoodName = cxo.getStringAttribute(Attribute.NAME);

                Trace trace = Trace.Utils.loadTrace(reader, likelihoodName);
                reader.close();
                if (trace == null)
                    throw new XMLParseException("Column '" + likelihoodName + "' can not be found for " + getParserName() + " element.");

                int maxState = Trace.Utils.getMaximumState(trace);

                int burnin = -1;
                if (xo.hasAttribute(BURN_IN)) {
                    // leaving the burnin attribute off will result in 10% being used
                    burnin = xo.getIntegerAttribute(BURN_IN);
                }

                if (burnin == -1) {
                    burnin = maxState / 10;
                }

                if (burnin < 0 || burnin >= maxState) {
                    burnin = maxState / 10;
                    System.out.println("WARNING: Burn-in larger than total number of states - using to 10%");
                }

                boolean harmonicOnly = false;
                if (cxo.hasAttribute(ONLY_HARMONIC))
                    harmonicOnly = cxo.getBooleanAttribute(ONLY_HARMONIC);

                boolean doBootstrap = true;
                if (cxo.hasAttribute(DO_BOOTSTRAP))
                    doBootstrap = cxo.getBooleanAttribute(DO_BOOTSTRAP);

                int bootstrapLength = 1000;
                if (cxo.hasAttribute(BOOTSTRAP_LENGTH))
                    bootstrapLength = cxo.getIntegerAttribute(BOOTSTRAP_LENGTH);


                MarginalLikelihoodAnalysis analysis = new MarginalLikelihoodAnalysis(trace,
                        burnin, harmonicOnly, doBootstrap, bootstrapLength);

                System.out.println(analysis.toString());

                return analysis;

            } catch (java.io.IOException ioe) {
                throw new XMLParseException(ioe.getMessage());
            }
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "Performs a trace analysis. Estimates the mean of the various statistics in the given log file.";
        }

        public Class getReturnType() {
            return MarginalLikelihoodAnalysis.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                new StringAttributeRule(FILE_NAME, "The name of a BEAST log file (can not include trees, which should be logged separately"),
                AttributeRule.newIntegerRule("burnIn", true)
                //, "The number of states (not sampled states, but actual states) that are discarded from the beginning of the trace before doing the analysis" ),
        };
    };
}
