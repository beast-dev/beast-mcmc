package dr.evomodel.coalescent;

import dr.evolution.tree.Tree;
import dr.evomodelxml.coalescent.BNPRLikelihoodParser;
import dr.inference.model.MatrixParameter;
import dr.inference.model.Parameter;
import dr.math.MathUtils;

import java.util.Arrays;
import java.util.List;

/**
 * Created by mkarcher on 9/15/16.
 */
public class BNPRLikelihood extends GMRFSkyrideLikelihood {

    protected double cutOff;
    protected int numGridPoints;
    protected double[] gridPoints;

    protected double precAlpha;
    protected double precBeta;

    protected Parameter samplingBetas;

    protected double logSamplingLikelihood;
    protected double storedLogSamplingLikelihood;

    protected double[] samplingTimes;
    protected boolean samplingTimesKnown;

    public BNPRLikelihood(List<Tree> treeList,
                          Parameter popParameter,
                          Parameter groupParameter,
                          Parameter precParameter,
                          Parameter lambda,
                          //Parameter beta,
                          MatrixParameter dMatrix,
                          boolean timeAwareSmoothing,
                          //boolean rescaleByRootHeight,
                          double cutOff,
                          int numGridPoints) {

        super(BNPRLikelihoodParser.BNPR_LIKELIHOOD);

        this.popSizeParameter = popParameter;
        this.groupSizeParameter = groupParameter;
        this.precisionParameter = precParameter;
        this.lambdaParameter = lambda;
        //this.betaParameter = beta;
        this.dMatrix = dMatrix;
        this.timeAwareSmoothing = timeAwareSmoothing;
        //this.rescaleByRootHeight = rescaleByRootHeight;

        this.cutOff = cutOff;
        this.numGridPoints = numGridPoints;
        setupGridPoints();

        addVariable(popSizeParameter);
        addVariable(precisionParameter);
        addVariable(lambdaParameter);
        //if (betaParameter != null) {
        //    addVariable(betaParameter);
        //}

        setTree(treeList);

        int correctFieldLength = getCorrectFieldLength();

        if (popSizeParameter.getDimension() <= 1) {
            // popSize dimension hasn't been set yet, set it here:
            popSizeParameter.setDimension(correctFieldLength);
        }

        fieldLength = popSizeParameter.getDimension();
        if (correctFieldLength != fieldLength) {
            throw new IllegalArgumentException("Population size parameter should have length " + correctFieldLength);
        }

        // Field length must be set by this point
        wrapSetupIntervals();
        coalescentIntervals = new double[fieldLength];
        storedCoalescentIntervals = new double[fieldLength];
        sufficientStatistics = new double[fieldLength];
        storedSufficientStatistics = new double[fieldLength];

        setupGMRFWeights();

        addStatistic(new DeltaStatistic());

        initializationReport();

        /* Force all entries in groupSizeParameter = 1 for compatibility with Tracer */
        if (groupSizeParameter != null) {
            for (int i = 0; i < groupSizeParameter.getDimension(); i++)
                groupSizeParameter.setParameterValue(i, 1.0);
        }
    }

    protected int getCorrectFieldLength() {
        return numGridPoints + 1;
    }

    protected void setupGridPoints() {
        if (gridPoints == null) {
            gridPoints = new double[numGridPoints];
        } else {
            Arrays.fill(gridPoints, 0);
        }

        for (int pt = 0; pt < numGridPoints; pt++) {
            gridPoints[pt] = (pt + 1) * (cutOff / numGridPoints);
        }
    }

    protected double[] getSamplingTimes() {
        if (!samplingTimesKnown) {
            Tree tree = this.tree;
            int n = tree.getExternalNodeCount();
            double[] nodeHeights = new double[n];
            double maxHeight = 0;
            samplingTimes = new double[n];

            for (int i = 0; i < n; i++) {
                nodeHeights[i] = tree.getNodeHeight(tree.getExternalNode(i));
                if (nodeHeights[i] > maxHeight) {
                    maxHeight = nodeHeights[i];
                }
            }

            for (int i = 0; i < n; i++) {
                samplingTimes[i] = maxHeight - nodeHeights[i];
            }
        }

        return samplingTimes;
    }

    private double[] diff(double[] doubles) {
        double[] result = new double[doubles.length - 1];

        for (int i = 0; i < doubles.length - 1; i++) {
            result[i] = doubles[i+1] - doubles[i];
        }

        return result;
    }

    private int[] bin(double[] data, double[] grid) {
        int[] result = new int[grid.length - 1];

        for (int i = 0; i < grid.length - 1; i++) {
            for (int j = 0; j < data.length; j++) {
                if (data[j] >= grid[i] && data[j] < grid[i + 1]) {
                    result[i]++;
                }
            }
        }
        return result;
    }

    private int[] binNA(int[] data) {
        int[] result = new int[data.length];
        boolean encountered = false;

        for (int i = data.length - 1; i >= 0; i--) {
            if (encountered || data[i] > 0) {
                encountered = true;
                result[i] = data[i];
            } else {
                result[i] = -1;
            }
        }

        return result;
    }

    public double calculateLogSamplingLikelihood() {

        if (!intervalsKnown) {
            // intervalsKnown -> false when handleModelChanged event occurs in super.
            wrapSetupIntervals();
            setupSufficientStatistics();
            intervalsKnown = true;
        }

        double[] currentGamma = popSizeParameter.getParameterValues();
        double[] currentBetas = samplingBetas.getParameterValues();

        int[] binnedSamples = binNA(bin(samplingTimes, gridPoints));
        double[] gridDiff = diff(gridPoints);

        double beta0 = currentBetas[0];
        double beta1 = currentBetas[1];

        double currentLike = beta0 * samplingTimes.length;
        for (int i = 0; binnedSamples[i] != -1 && i < fieldLength; i++) {
            currentLike += beta1 * currentGamma[i] * binnedSamples[i] - gridDiff[i] * Math.exp(beta0 + beta1 * currentGamma[i]);
        }

        return currentLike;
    }

    public double calculateLogSamplingLikelihoodSubGamma(double[] gamma) {

        if (!intervalsKnown) {
            // intervalsKnown -> false when handleModelChanged event occurs in super.
            wrapSetupIntervals();
            setupSufficientStatistics();
            intervalsKnown = true;
        }

        double[] currentGamma = gamma;
        double[] currentBetas = samplingBetas.getParameterValues();

        int[] binnedSamples = binNA(bin(samplingTimes, gridPoints));
        double[] gridDiff = diff(gridPoints);

        double beta0 = currentBetas[0];
        double beta1 = currentBetas[1];

        double currentLike = beta0 * samplingTimes.length;
        for (int i = 0; binnedSamples[i] != -1 && i < fieldLength; i++) {
            currentLike += beta1 * currentGamma[i] * binnedSamples[i] - gridDiff[i] * Math.exp(beta0 + beta1 * currentGamma[i]);
        }

        return currentLike;
    }

    public double calculateLogSamplingLikelihoodSubBetas(double[] betas) {

        if (!intervalsKnown) {
            // intervalsKnown -> false when handleModelChanged event occurs in super.
            wrapSetupIntervals();
            setupSufficientStatistics();
            intervalsKnown = true;
        }

        double[] currentGamma = popSizeParameter.getParameterValues();

        int[] binnedSamples = binNA(bin(samplingTimes, gridPoints));
        double[] gridDiff = diff(gridPoints);

        double beta0 = betas[0];
        double beta1 = betas[1];

        double currentLike = beta0 * samplingTimes.length;
        for (int i = 0; binnedSamples[i] != -1 && i < fieldLength; i++) {
            currentLike += beta1 * currentGamma[i] * binnedSamples[i] - gridDiff[i] * Math.exp(beta0 + beta1 * currentGamma[i]);
        }

        return currentLike;
    }

    public double getLogLikelihood() {
        if (!likelihoodKnown) {
            logLikelihood = calculateLogCoalescentLikelihood();
            logFieldLikelihood = calculateLogFieldLikelihood();
            logSamplingLikelihood = calculateLogSamplingLikelihood();
            likelihoodKnown = true;
        }

        return logLikelihood + logFieldLikelihood + logSamplingLikelihood;
    }

    public double calculateLogCoalescentLikelihoodSubGamma(double[] gamma) {

        if (!intervalsKnown) {
            // intervalsKnown -> false when handleModelChanged event occurs in super.
            wrapSetupIntervals();
            setupSufficientStatistics();
            intervalsKnown = true;
        }

        // Matrix operations taken from block update sampler to calculate data likelihood and field prior

        double currentLike = 0;
        double[] currentGamma = gamma;

        for (int i = 0; i < fieldLength; i++) {
            currentLike += -numCoalEvents[i] * currentGamma[i] + ploidySums[i] - sufficientStatistics[i] * Math.exp(-currentGamma[i]);
        }

        return currentLike;
    }

    public double calculateLogFieldLikelihoodSubGamma(double[] gamma) {

    }

    public double getLogLikelihoodSubGamma(double[] gamma) {
        if (!likelihoodKnown) {
            logLikelihood = calculateLogCoalescentLikelihoodSubGamma(gamma);
            logFieldLikelihood = calculateLogFieldLikelihoodSubGamma(gamma);
            logSamplingLikelihood = calculateLogSamplingLikelihoodSubGamma(gamma);
            likelihoodKnown = true;
        }

        return logLikelihood + logFieldLikelihood + logSamplingLikelihood;
    }

    public int getNumGridPoints() {
        return this.gridPoints.length;
    }

    public double getPrecAlpha() {
        return precAlpha;
    }

    public void setPrecAlpha(double precAlpha) {
        this.precAlpha = precAlpha;
    }

    public double getPrecBeta() {
        return precBeta;
    }

    public void setPrecBeta(double precBeta) {
        this.precBeta = precBeta;
    }

    public Parameter getSamplingBetas() {
        return samplingBetas;
    }
}
