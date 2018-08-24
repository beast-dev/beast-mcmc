package dr.evomodel.coalescent;

import dr.evolution.coalescent.DemographicFunction;
import dr.evolution.io.NewickImporter;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.util.Units;
import dr.evomodel.tree.TreeModel;
import dr.evomodelxml.coalescent.BNPRSamplingLikelihoodParser;
import dr.inference.model.*;
import dr.math.matrixAlgebra.Matrix;
import dr.util.Author;
import dr.util.Citable;
import dr.util.Citation;
import org.apache.commons.math.FunctionEvaluationException;
import org.apache.commons.math.MaxIterationsExceededException;
import org.apache.commons.math.analysis.UnivariateRealFunction;
import org.apache.commons.math.analysis.integration.RombergIntegrator;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;

/**
 * Created by mkarcher on 3/9/17.
 */
public class BNPRSamplingLikelihood extends AbstractModelLikelihood implements Citable {
    private boolean likelihoodKnown = false;
    private boolean storedLikelihoodKnown = false;
    private double logLikelihood;
    private double storedLogLikelihood;

    private boolean samplingTimesKnown = false;
    private boolean storedSamplingTimesKnown = false;
    private double[] samplingTimes;
    private double[] storedSamplingTimes;
    private double[] logPopSizes;
    private double[] storedLogPopSizes;

    private double[] epochWidths = null;
    private double[] midpoints = null;

    private Tree tree;
    private int numSamples;
    private Parameter betas;
    private DemographicModel population;
    private MatrixParameter covariates;
    private MatrixParameter powerCovariates;
    private Parameter powerBetas;
    private DemographicWrapper demo;
    private DemographicModel demoModel = null;
    private RombergIntegrator integrator;

    public BNPRSamplingLikelihood(Tree tree, Parameter betas, DemographicModel population) {
        this(tree, betas, population, null, null);
    }

    public BNPRSamplingLikelihood(Tree tree, Parameter betas, DemographicModel population, double[] epochWidths) {
        this(tree, betas, population, epochWidths, null);
    }

    public BNPRSamplingLikelihood(Tree tree, Parameter betas, DemographicModel population, double[] epochWidths,
                                  MatrixParameter covariates) {
        this(tree, betas, population, epochWidths, covariates, null, null);
    }

    public BNPRSamplingLikelihood(Tree tree, Parameter betas, DemographicModel population, double[] epochWidths,
                                  MatrixParameter covariates, MatrixParameter powerCovariates, Parameter powerBetas) {
        super(BNPRSamplingLikelihoodParser.SAMPLING_LIKELIHOOD);

        this.tree = tree;
        this.betas = betas;
        this.population = population;
        this.covariates = covariates;
        this.powerCovariates = powerCovariates;
        this.powerBetas = powerBetas;

        this.likelihoodKnown = false;
        this.samplingTimesKnown = false;

        this.demo = new DemographicWrapper(population.getDemographicFunction());
        this.demoModel = population;
        this.integrator = new RombergIntegrator();

        this.numSamples = tree.getExternalNodeCount();
        samplingTimes = new double[numSamples];
        storedSamplingTimes = new double[numSamples];
        logPopSizes = new double[numSamples];
        storedLogPopSizes = new double[numSamples];

//        this.covariates = covariates;

        if (epochWidths != null) {
            setEpochs(epochWidths);
        }

        if (tree instanceof TreeModel) {
            addModel((TreeModel) tree);
        }

        addModel(demoModel);

        if (this.betas != null) {
            addVariable(this.betas);
        }

        if (this.powerBetas != null) {
            addVariable(this.powerBetas);
        }

        if (this.covariates != null) {
            addVariable(this.covariates);
        }

        if (this.powerCovariates != null) {
            addVariable(this.powerCovariates);
        }

        setupSamplingTimes();
    }

    public void setEpochs(double[] epochWidths) {
        this.epochWidths = epochWidths;

        if (epochWidths == null) {
            this.midpoints = null;
        } else {
            this.midpoints = new double[epochWidths.length + 1];
            double start = 0;
            for (int i = 0; i < epochWidths.length; i++) {
                this.midpoints[i] = start + epochWidths[i] / 2.0;
                start += epochWidths[i];
            }
            this.midpoints[epochWidths.length] = start + epochWidths[epochWidths.length - 1] / 2.0;
        }
    }

    private void setupSamplingTimes() {
//        this.numSamples = tree.getExternalNodeCount();
//        samplingTimes = new double[numSamples];
//        logPopSizes = new double[numSamples];

        for (int i = 0; i < numSamples; i++) {
            NodeRef node = tree.getExternalNode(i);
            samplingTimes[i] = tree.getNodeHeight(node);
        }

        Arrays.sort(samplingTimes);

        samplingTimesKnown = true;
    }

    public double[] getSamplingTimes() {
        if (!samplingTimesKnown) {
            setupSamplingTimes();
        }

        return samplingTimes;
    }

    @Override
    public Model getModel() {
        return this;
    }

    @Override
    public double getLogLikelihood() {
        if (!samplingTimesKnown) {
            setupSamplingTimes();
        }

        if (!likelihoodKnown) {
            logLikelihood = calculateLogLikelihood();
            likelihoodKnown = true;
        }

        return logLikelihood;
    }

    public double calculateLogLikelihood() {
        double minSample = samplingTimes[0];
        double maxSample = samplingTimes[numSamples-1];

        double beta0 = betas.getParameterValue(0);
        double beta1 = betas.getParameterValue(1);

        double logLik = numSamples * beta0; // Start with beta0 part of likelihood.

//        System.out.printf("beta0 part = %f\n", logLik);

        DemographicFunction f = population.getDemographicFunction();

        // Gather log-population sizes to facilitate likelihood calculation.
        // Separated to make cut-paste into a new setupLogPopSizes() function easy.
        for (int i = 0; i < numSamples; i++) {
            logPopSizes[i] = f.getLogDemographic(samplingTimes[i]);
        }

        double eventLogLik = 0.0;
        // Calculate the event component of the inhomogeneous Poisson process log-likelihood
        for (int i = 0; i < numSamples; i++) {
            eventLogLik += beta1 * logPopSizes[i];

            if (this.powerCovariates != null) {
                for (int j = 0; j < powerBetas.getDimension(); j++) {
                    eventLogLik += powerBetas.getParameterValue(j) * evaluatePiecewise(this.epochWidths, powerCovariates.getColumnValues(j), samplingTimes[i]) * logPopSizes[i];
                }
            }

            if (this.covariates != null) {
                for (int j = 2; j < betas.getDimension(); j++) {
//                    eventLogLik += betas.getParameterValue(j) * covariates.getParameterValue(i, j - 2);
                    eventLogLik += betas.getParameterValue(j) * evaluatePiecewise(this.epochWidths, covariates.getColumnValues(j-2), samplingTimes[i]);
                }
            }
        }

//        System.out.printf("event part = %f\n", eventLogLik);

        logLik += eventLogLik;

        double integralLogLik = 0.0;

//        if (epochWidths == null) {
//            demo.setF(f); // Possibly unnecessary if the DemographicFunction reference stays the same after model changes
//            PowerFunctional pf = new PowerFunctional(demo, beta1);
//
//            // Calculate the integral component of the inhomogeneous Poisson process log-likelihood
//            try {
//                integralLogLik -= Math.exp(beta0) * integrator.integrate(pf, minSample, maxSample);
//            } catch (MaxIterationsExceededException e) {
//                throw new RuntimeException(e);
//            } catch (FunctionEvaluationException e) {
//                throw new RuntimeException(e);
//            }
//        } else {
//            double[] fBeta = new double[this.midpoints.length];
//            for (int i = 0; i < fBeta.length; i++) {
//                fBeta[i] = Math.exp(f.getLogDemographic(this.midpoints[i]) * beta1);
//
//                if (this.covariates != null) {
//                    for (int j = 2; j < betas.getDimension(); j++) {
//                        fBeta[i] *= Math.exp(betas.getParameterValue(j) * this.covariates.getParameterValue(j - 2, i));
//                    }
//                }
//            }
//            integralLogLik -= Math.exp(beta0) * integratePiecewise(this.epochWidths, fBeta, minSample, maxSample);
//        }

        double[] fBeta = new double[this.midpoints.length];
        double fBetaLog;
        for (int i = 0; i < fBeta.length; i++) {
            fBetaLog = beta1 * f.getLogDemographic(this.midpoints[i]);

            if (this.powerCovariates != null) {
                for (int j = 0; j < powerBetas.getDimension(); j++) {
                    fBetaLog += powerBetas.getParameterValue(j) * this.powerCovariates.getParameterValue(i, j) * f.getLogDemographic(this.midpoints[i]);
                }
            }

            if (this.covariates != null) {
                for (int j = 2; j < betas.getDimension(); j++) {
                    fBetaLog += betas.getParameterValue(j) * this.covariates.getParameterValue(i, j - 2);
                }
            }

            fBeta[i] = Math.exp(fBetaLog);
        }
        integralLogLik -= Math.exp(beta0) * integratePiecewise(this.epochWidths, fBeta, minSample, maxSample);


//        System.out.printf("integral part = %f\n", integralLogLik);
        logLik += integralLogLik;

        return logLik;
    }

    private static double evaluatePiecewise(double[] widths, double[] heights, double t) {
        double[] widths2 = new double[widths.length + 1];
        System.arraycopy(widths, 0, widths2, 0, widths.length);
        widths2[widths.length] = Double.POSITIVE_INFINITY;

//        double result = 0.0;
        int i = 0;
        while (t > widths2[i]) {
            t -= widths2[i];
            i++;
        }

        return heights[i];
    }

    private static double integratePiecewise(double[] widths, double[] heights, double start, double end) {
        double[] widths2 = new double[widths.length + 1];
        System.arraycopy(widths, 0, widths2, 0, widths.length);
        widths2[widths.length] = Double.POSITIVE_INFINITY;

        double result = 0.0;
        int i = 0;
        while (start > widths2[i]) {
            start -= widths2[i];
            end -= widths2[i];
            i++;
        }

        if (end < widths2[i]) {
            result = (end - start) * heights[i];
        } else {
            result += (widths2[i] - start) * heights[i];
            end -= widths2[i];
            i++;

            while (end > widths2[i]) {
                result += widths2[i] * heights[i];
                end -= widths2[i];
                i++;
            }

            if (end > 0.0) {
                result += end * heights[i];
            }
        }

        return result;
    }

    @Override
    public void makeDirty() {
        likelihoodKnown = false;
        samplingTimesKnown = false;
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
        makeDirty();
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        makeDirty();
    }

    @Override
    protected void storeState() {
        // copy the sampling times into the storedSamplingTimes
        System.arraycopy(samplingTimes, 0, storedSamplingTimes, 0, samplingTimes.length);

        storedSamplingTimesKnown = samplingTimesKnown;
        storedLikelihoodKnown = likelihoodKnown;
        storedLogLikelihood = logLikelihood;

        System.arraycopy(logPopSizes, 0, storedLogPopSizes, 0, logPopSizes.length);
    }

    @Override
    protected void restoreState() {
        // copy the sampling times into the storedSamplingTimes
        System.arraycopy(storedSamplingTimes, 0, samplingTimes, 0, storedSamplingTimes.length);

        samplingTimesKnown = storedSamplingTimesKnown;
        likelihoodKnown = storedLikelihoodKnown;
        logLikelihood = storedLogLikelihood;

        System.arraycopy(storedLogPopSizes, 0, logPopSizes, 0, storedLogPopSizes.length);
    }

    @Override
    protected void acceptState() {
        // nothing to do
    }

    public String toString() {
        return Double.toString(logLikelihood);
    }

    @Override
    public Citation.Category getCategory() {
        return Citation.Category.PRIOR_MODELS;
    }

    @Override
    public String getDescription() {
        return "Bayesian non-parametric preferential sampling";
    }

    @Override
    public List<Citation> getCitations() {
        return Collections.singletonList(CITATION);
    }

    public static Citation CITATION = new Citation(
            new Author[]{
                    new Author("MD", "Karcher"),
                    new Author("MA", "Suchard"),
                    new Author("G", "Dudas"),
                    new Author( "T", "Bedford"),
                    new Author("VN","Minin")
            },
//            "Estimating effective population size changes from preferentially sampled genetic sequences",
            Citation.Status.IN_PREPARATION
    );

    private class PowerFunctional implements UnivariateRealFunction {
        private double power;
        private UnivariateRealFunction f;

        public PowerFunctional(UnivariateRealFunction f, double power) {
            this.f = f;
            this.power = power;
        }

        @Override
        public double value(double t) throws FunctionEvaluationException {
            return Math.pow(f.value(t), power);
        }

        public double getPower() {
            return power;
        }

        public void setPower(double power) {
            this.power = power;
        }

        public UnivariateRealFunction getF() {
            return f;
        }

        public void setF(UnivariateRealFunction f) {
            this.f = f;
        }
    }

    private class DemographicWrapper implements UnivariateRealFunction {
        private DemographicFunction f;

        public DemographicWrapper(DemographicFunction f) {
            this.f = f;
        }

        public double value(double t) {
            return f.getDemographic(t);
        }

        public void setF(DemographicFunction f) {
            this.f = f;
        }
    }

    public static void main(String[] args) throws Exception {
//        DistributionLikelihood lik = new DistributionLikelihood(new NormalDistribution(0, 1));
//        Parameter par = new Parameter.Default(new double[]{0.0, 1.0});
//        lik.addData(par);
//        System.out.println("Likelihood = " + lik.getLogLikelihood());

        NewickImporter importer = new NewickImporter("(((t93_0:0.1823217653,t94_0:0.08702282874):0.853414946,t95_0:0.5591826007):16.46923208,((t96_0:4.730080227,(t97_0:3.316833,((t98_0:1.117259784,((t84_0:0.02918861288,((((t77_0:0.6269405877,t78_0:0.1883673111):0.1035487414,(((((t58_0:0.366605532,(t59_0:0.2902223473,t60_0:0.09330035051):0.05814249819):0.2687774288,(t61_0:0.03909732613,((((t31_0:0.4925976152,(t1_0:1.543243834,(t2_0:0.3786992111,t3_0:0.286843186):0.9369069124):1.84806174):1.388380979,(((((t4_0:1.57060332,t5_0:1.565578151):0.9340116213,(t6_0:0.6129297388,t7_0:0.5252059248):1.883951667):1.298791879,(t8_0:3.070277733,((t9_0:2.351705269,t10_0:2.318654989):0.03364809509,t11_0:1.986711821):0.6207718581):0.558736892):0.1439273897,((((t12_0:0.3565815494,t13_0:0.2420866743):1.027508236,t14_0:1.128497322):0.1217011969,t15_0:1.211013892):0.5793893699,((t16_0:0.4746827902,t17_0:0.4609704827):1.049116131,t18_0:1.497407856):0.1530705075):1.000179818):0.0008681768441,(t19_0:2.204632212,((t20_0:1.857557731,(t21_0:1.005514944,t22_0:0.9291487706):0.6230458901):0.1109842512,t23_0:1.582286036):0.1992671825):0.4433270097):0.3343580536):0.6454071651,(((t24_0:0.815171935,t25_0:0.6918964743):0.2147852872,(t26_0:0.1228387289,t27_0:0.05280374254):0.5882964168):2.025084596,(t28_0:0.2705161616,t29_0:0.2642467551):2.230464942):0.1103642803):6.832809449,(((t30_0:0.9906878227,t32_0:0.1856034526):1.311909106,t33_0:0.3458898346):6.842038256,t34_0:4.353598435):0.2887890018):0.414867716):0.05823058913):0.3120716399,(t35_0:1.623697863,((t36_0:1.029267152,t37_0:0.9940333973):0.2630810985,t38_0:1.077299678):0.2899006266):2.324485801):1.915577302,(t39_0:1.748484212,t40_0:1.720294263):3.606542332):0.6429878862,(t41_0:5.650860772,((t42_0:0.8033637796,t43_0:0.3493332566):3.763664014,(t44_0:0.8778820702,(t45_0:0.5483027215,t46_0:0.4562980629):0.1491010534):3.197472308):0.8000617359):0.2857177615):0.4150148533):0.7376635867,(((t47_0:0.3446780656,t48_0:0.3130649187):3.810586451,(((((t49_0:0.9994210908,t50_0:0.8953272412):0.5705160932,(t51_0:0.7994882815,t52_0:0.7164389879):0.6268954443):0.34752079,t53_0:1.598405188):0.2494469601,(t54_0:1.508633144,t55_0:1.473482952):0.05756852234):0.6596007127,(t56_0:1.076711525,t57_0:1.070680706):1.008030326):1.241809834):0.7439478609,((t62_0:0.1842275742,t63_0:0.1339199725):1.294985145,(t64_0:1.040418148,t65_0:0.7556259274):0.3408415617):1.72091551):0.5986199428):2.461710024,(((t66_0:3.108881856,((t67_0:2.000699305,((t68_0:1.382538601,t69_0:1.246844479):0.249533307,t70_0:1.366941817):0.176205828):0.5708242903,(t71_0:0.9490874836,t72_0:0.7691535035):0.8078785817):0.475066817):0.06683807466,((t73_0:0.7617022066,t74_0:0.7482687036):0.1047558487,t75_0:0.7677691921):0.8692021164):0.1592569402,t76_0:1.722131386):2.262393917):2.572913635):2.595279242,(t79_0:3.600210891,(t80_0:0.7070684518,t81_0:0.3253173205):2.684283986):2.510617438):0.5734740762):1.936674132,t82_0:6.475760804):0.2353584188):1.216419756):0.9696225718,(((t83_0:0.4830428417,t85_0:0.214507183):4.223119782,(((t86_0:1.283242236,t87_0:1.278786271):1.710456,t88_0:2.295501437):0.6485241096,(t89_0:1.619761915,t90_0:1.575819491):1.230444387):0.7762895827):1.668960635,(t91_0:4.403654576,t92_0:4.367902389):0.7940511151):1.277021867):11.16032274)");
        Tree tree = importer.importNextTree();
        Parameter betas = new Parameter.Default(new double[] {2.430198, 1.121209, 18.963120});
        Parameter N0 = new Parameter.Default(new double[] {0.957477459020635, 2.43878447702397, 1.16178606305887, 2.47803465457009, 3.59425875085121, 5.99484118979476, 8.45281123199342, 5.32731462686442, 2.14530385537217, 1.32838970512892, 1.57682364008761, 2.06977561883114, 0.975318782259623, 0.729019262672561, 0.528122712339147, 1.73826551939609, 1.59115234812228, 0.991767268014354, 1.19544144384689, 4.2431612093277, 9.68837537548171, 21.7809958450516, 72.9542238498992, 67.4709298097791, 17.5930802539968, 17.9641115502897, 20.2230985692042, 18.5751760066478, 15.316828425915, 53.4144639187011, 20.6589204576208, 11.4695515805863, 25.1358546005526, 41.1222703999958, 36.0749583093113, 43.5582733013315, 191.509377345225, 60.2774566116462, 13.6461065289052, 29.0150510257721, 11.8060344563918, 14.0771159483723, 17.2536924114547, 15.3833722037801, 46.48760509832, 41.881067569086, 27.6941667792898, 55.4922963887114, 253.986328605972, 305.304267851333, 243.190770701089, 144.095172595153, 172.765640994598, 174.182113284746, 454.771916114004, 514.94845166669, 266.355413273477, 2039.89326970745, 972.728248184194, 305.635176732534, 227.134976673311, 143.297859235514, 74.5162610231541, 45.3481762741632, 32.1527697118794, 7.6456190280646, 12.1418530975389, 9.48239531687115, 7.37481430050462, 4.49432342829254, 14.6077380251641, 13.5858667415758, 11.2001057413539, 4.12650576804321, 3.83168028086946, 5.56705828329053, 7.86450990980968, 5.6493737041167, 4.76307833501312, 6.45317992085796, 9.78615502142029, 13.8294443012538, 22.5127519861243, 28.9164141530354, 85.3528130029434, 34.8658623017608, 47.1574370857575, 40.2205642235078, 57.50963654004, 65.3719036601996, 44.8775886163862, 17.6951867195527, 9.07043638586854, 18.663530289391, 15.6125065593161, 2.12354416497279, 4.96867966440677, 3.63537271117703, 2.09626570198215, 1.03533797751651});
        double[] epochLengths = new double[] {0.404760913356317, 0.404760913356317, 0.404760913356317, 0.404760913356317, 0.404760913356317, 0.404760913356317, 0.404760913356317, 0.404760913356317, 0.404760913356317, 0.404760913356317, 0.404760913356317, 0.404760913356317, 0.404760913356317, 0.404760913356317, 0.404760913356317, 0.404760913356317, 0.404760913356317, 0.404760913356317, 0.404760913356317, 0.404760913356317, 0.404760913356317, 0.404760913356317, 0.404760913356317, 0.404760913356317, 0.404760913356317, 0.404760913356317, 0.404760913356317, 0.404760913356317, 0.404760913356317, 0.404760913356317, 0.404760913356317, 0.404760913356317, 0.404760913356317, 0.404760913356317, 0.404760913356317, 0.404760913356317, 0.404760913356317, 0.404760913356317, 0.404760913356317, 0.404760913356317, 0.404760913356317, 0.404760913356317, 0.404760913356317, 0.404760913356317, 0.404760913356317, 0.404760913356317, 0.404760913356317, 0.404760913356317, 0.404760913356317, 0.404760913356317, 0.404760913356317, 0.404760913356317, 0.404760913356317, 0.404760913356317, 0.404760913356317, 0.404760913356317, 0.404760913356317, 0.404760913356317, 0.404760913356317, 0.404760913356317, 0.404760913356317, 0.404760913356317, 0.404760913356317, 0.404760913356317, 0.404760913356317, 0.404760913356317, 0.404760913356317, 0.404760913356317, 0.404760913356317, 0.404760913356317, 0.404760913356317, 0.404760913356317, 0.404760913356317, 0.404760913356317, 0.404760913356317, 0.404760913356317, 0.404760913356317, 0.404760913356317, 0.404760913356317, 0.404760913356317, 0.404760913356317, 0.404760913356317, 0.404760913356317, 0.404760913356317, 0.404760913356317, 0.404760913356317, 0.404760913356317, 0.404760913356317, 0.404760913356317, 0.404760913356317, 0.404760913356317, 0.404760913356317, 0.404760913356317, 0.404760913356317, 0.404760913356317, 0.404760913356317, 0.404760913356317, 0.404760913356317, 0.404760913356317};
        Parameter covariate1 = new Parameter.Default(new double[]{-0.00505951141695401, -0.015178534250862, -0.0252975570847698, -0.0354165799186778, -0.0455356027525857, -0.0556546255864937, -0.0657736484204016, -0.0758926712543096, -0.0860116940882175, -0.0961307169221254, -0.106249739756033, -0.116368762589941, -0.126487785423849, -0.136606808257757, -0.146725831091665, -0.156844853925573, -0.166963876759481, -0.177082899593389, -0.187201922427297, -0.197320945261205, -0.207439968095113, -0.217558990929021, -0.227678013762929, -0.237797036596836, -0.247916059430744, -0.258035082264652, -0.26815410509856, -0.278273127932468, -0.288392150766376, -0.298511173600284, -0.308630196434192, -0.3187492192681, -0.328868242102008, -0.338987264935916, -0.349106287769824, -0.359225310603732, -0.36934433343764, -0.379463356271548, -0.389582379105455, -0.399701401939363, -0.409820424773271, -0.419939447607179, -0.430058470441087, -0.440177493274995, -0.450296516108903, -0.460415538942811, -0.470534561776719, -0.480653584610627, -0.490772607444535, -0.500891630278443, -0.511010653112351, -0.521129675946259, -0.531248698780167, -0.541367721614075, -0.551486744447982, -0.56160576728189, -0.571724790115798, -0.581843812949706, -0.591962835783614, -0.602081858617522, -0.61220088145143, -0.622319904285338, -0.632438927119246, -0.642557949953154, -0.652676972787062, -0.66279599562097, -0.672915018454878, -0.683034041288786, -0.693153064122694, -0.703272086956601, -0.713391109790509, -0.723510132624417, -0.733629155458325, -0.743748178292233, -0.753867201126141, -0.763986223960049, -0.774105246793957, -0.784224269627865, -0.794343292461773, -0.804462315295681, -0.814581338129589, -0.824700360963497, -0.834819383797405, -0.844938406631313, -0.85505742946522, -0.865176452299128, -0.875295475133036, -0.885414497966944, -0.895533520800852, -0.90565254363476, -0.915771566468668, -0.925890589302576, -0.936009612136484, -0.946128634970392, -0.9562476578043, -0.966366680638208, -0.976485703472116, -0.986604726306024, -0.996723749139932, -1.00684277197384});
        MatrixParameter covariates = new MatrixParameter("Covariates", new Parameter[]{covariate1});
//        Parameter N0 = new Parameter.Default(new double[] {1.0, 0.5, 1.0, 1.0, 10.0});
//        double[] epochLengths = new double[] {0.35, 0.35, 0.35, 0.35};

//        int numGridPoints = 2;

        DemographicModel Ne = new PiecewisePopulationModel("Ne(t)", N0, epochLengths, false, Units.Type.DAYS);

//        System.out.println("Integration: " + integratePiecewise(epochLengths, N0.getParameterValues(), 3.0, 14.0));

//        System.out.println("Got here A");

//        BNPRSamplingLikelihood samplingLikelihood1 = new BNPRSamplingLikelihood(tree, betas, Ne);
        BNPRSamplingLikelihood samplingLikelihood = new BNPRSamplingLikelihood(tree, betas, Ne, epochLengths, covariates);

//        System.out.println("Sampling times: " + Arrays.toString(samplingLikelihood.getSamplingTimes()));
//
        System.out.println("Sampling likelihood: " + samplingLikelihood.getLogLikelihood());
//
//        samplingLikelihood.makeDirty();
//
//        System.out.println("Sampling likelihood 2: " + samplingLikelihood2.getLogLikelihood());
//
//        betas.setParameterValue(0, 0.0);
//
//        System.out.println("Sampling likelihood: " + samplingLikelihood.getLogLikelihood());
//
//        samplingLikelihood.makeDirty();
//
//        System.out.println("Sampling likelihood: " + samplingLikelihood.getLogLikelihood());
//
//        System.out.println("Midpoints: " + Arrays.toString(samplingLikelihood.midpoints));

//        Parameter row1 = new Parameter.Default(new double[]{1.0, 2.0, 3.0});
//        MatrixParameter mat = new MatrixParameter("Matrix");
//        mat.addParameter(row1);

//        System.out.println(mat.getParameterValue(1,0));
    }
}
