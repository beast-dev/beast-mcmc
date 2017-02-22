package dr.evomodel.coalescent;

import dr.evolution.coalescent.IntervalType;
import dr.evolution.coalescent.TreeIntervals;
import dr.evolution.io.NewickImporter;
import dr.evolution.tree.Tree;
import dr.evomodel.coalescent.operators.BNPRBlockUpdateOperator;
import dr.evomodel.tree.TreeModel;
import dr.evomodelxml.coalescent.BNPRLikelihoodParser;
import dr.inference.model.MatrixParameter;
import dr.inference.model.Parameter;
import no.uib.cipr.matrix.BandCholesky;
import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.SymmTridiagMatrix;
import no.uib.cipr.matrix.UpperSPDBandMatrix;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Created by mkarcher on 9/15/16.
 */
public class BNPRLikelihood extends GMRFSkyrideLikelihood {
    public static final double NUGGET = 0.0001;

    protected double cutOff;
    protected int numGridPoints;
    protected int oldFieldLength;
    protected double[] numCoalEvents;
    protected double[] storedNumCoalEvents;
    protected double[] gridPoints;
    protected double[] midPoints;


    protected double precAlpha;
    protected double precBeta;

    protected Parameter samplingBetas;

    protected double logSamplingLikelihood;
    protected double storedLogSamplingLikelihood;

    protected boolean samplingAware;
    protected double[] samplingTimes;
    protected boolean samplingTimesKnown;

    public BNPRLikelihood(List<Tree> treeList,
                          Parameter popParameter,
                          Parameter groupParameter,
                          Parameter precParameter,
                          Parameter lambda,
                          Parameter beta,
                          MatrixParameter dMatrix,
                          boolean timeAwareSmoothing,
                          boolean rescaleByRootHeight,
                          double cutOff,
                          int numGridPoints,
                          boolean samplingAware) {

        super(BNPRLikelihoodParser.BNPR_LIKELIHOOD);

        this.popSizeParameter = popParameter;
        this.groupSizeParameter = groupParameter;
        this.precisionParameter = precParameter;
        this.lambdaParameter = lambda;
        this.betaParameter = beta;
        this.dMatrix = dMatrix;
        this.timeAwareSmoothing = timeAwareSmoothing;
        this.rescaleByRootHeight = rescaleByRootHeight;

        this.cutOff = cutOff;
        this.numGridPoints = numGridPoints;
        this.samplingAware = samplingAware;

        int correctFieldLength = getCorrectFieldLength();

        if (popSizeParameter.getDimension() <= 1) {
            // popSize dimension hasn't been set yet, set it here:
            popSizeParameter.setDimension(correctFieldLength);
        }

        fieldLength = popSizeParameter.getDimension();
        if (correctFieldLength != fieldLength) {
            throw new IllegalArgumentException("Population size parameter should have length " + correctFieldLength);
        }

        setupGridPoints();
        setupMidPoints();

        addVariable(popSizeParameter);
        addVariable(precisionParameter);
        addVariable(lambdaParameter);
        //if (betaParameter != null) {
        //    addVariable(betaParameter);
        //}

        setTree(treeList);

        oldFieldLength = super.getCorrectFieldLength();

        samplingTimes = getSamplingTimes();
        samplingTimesKnown = true;

        // Field length must be set by this point
        wrapSetupIntervals();
        coalescentIntervals = new double[oldFieldLength];
        storedCoalescentIntervals = new double[oldFieldLength];
        sufficientStatistics = new double[fieldLength];
        storedSufficientStatistics = new double[fieldLength];
        numCoalEvents = new double[fieldLength];
        storedNumCoalEvents = new double[fieldLength];

        setupGPWeights();

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

    protected void setTree(List<Tree> treeList) {
        if (treeList.size() != 1) {
            throw new RuntimeException("BNPRLikelihood only implemented for one tree");
        }
        this.tree = treeList.get(0);
        this.treesSet = null;
        if (tree instanceof TreeModel) {
            addModel((TreeModel) tree);
        }
    }

    public void initializationReport() {
        System.out.println("Creating a BNPR model");
        System.out.println("\tPopulation sizes: " + popSizeParameter.getDimension());
    }

    public void wrapSetupIntervals() {
        setupIntervals();
    }

    protected void setupGridPoints() {
        if (gridPoints == null) {
            gridPoints = new double[numGridPoints];
        } else {
            Arrays.fill(gridPoints, 0);
        }

        for (int pt = 0; pt < numGridPoints; pt++) {
            gridPoints[pt] = pt * cutOff / (numGridPoints - 1);
        }
    }

    private void setupMidPoints() {
        if (gridPoints == null) {
            setupGridPoints();
        }

        midPoints = new double[fieldLength];

        midPoints[0] = 0;
        for (int i = 1; i < fieldLength - 1; i++) {
            midPoints[i] = (gridPoints[i-1] + gridPoints[i]) / 2.0;
        }
        midPoints[fieldLength - 1] = gridPoints[numGridPoints - 1];
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

    //***********************
    // Calculate Likelihood
    //***********************

    public double getLogLikelihood() {
        if (!likelihoodKnown) {
            logLikelihood = calculateLogCoalescentLikelihood();
            logFieldLikelihood = calculateLogFieldLikelihood();

            if (samplingAware)
                logSamplingLikelihood = calculateLogSamplingLikelihood();

            likelihoodKnown = true;
        }

        double ll = logLikelihood + logFieldLikelihood;
        if (samplingAware)
            ll += logSamplingLikelihood;

        return ll;
    }

    protected double calculateLogCoalescentLikelihood() {
//        makeIntervalsKnown();

        // Matrix operations taken from block update sampler to calculate data likelihood and field prior

        double currentLike = 0;
        double[] currentGamma = popSizeParameter.getParameterValues();

        for (int i = 0; i < fieldLength; i++) {
            System.out.println("Sufficient statistic " + i + ": " + sufficientStatistics[i]);
            currentLike += -currentGamma[i] - sufficientStatistics[i] * Math.exp(-currentGamma[i]);
        }

        return currentLike;// + LogNormalDistribution.logPdf(Math.exp(popSizeParameter.getParameterValue(coalescentIntervals.length - 1)), mu, sigma);
    }

    public double getLogLikelihoodSubGamma(double[] gamma) {
        double logLikelihood = calculateLogCoalescentLikelihoodSubGamma(gamma);
        double logFieldLikelihood = calculateLogFieldLikelihoodSubGamma(gamma);
        double logSamplingLikelihood = calculateLogSamplingLikelihoodSubGamma(gamma);

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
            currentLike += -numCoalEvents[i] * currentGamma[i] - sufficientStatistics[i] * Math.exp(-currentGamma[i]);
        }

        return currentLike;
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

    public double calculateLogFieldLikelihood() {

        if (!intervalsKnown) {
            //intervalsKnown -> false when handleModelChanged event occurs in super.
            wrapSetupIntervals();
            setupSufficientStatistics();
            intervalsKnown = true;
        }

        DenseVector diagonal1 = new DenseVector(fieldLength);
        DenseVector currentGamma = new DenseVector(popSizeParameter.getParameterValues());

        //updateGammaWithCovariates(currentGamma);

        //double currentLike = handleMissingValues();
        double currentLike = 0.0;

        SymmTridiagMatrix currentQ = getScaledWeightMatrix(precisionParameter.getParameterValue(0), lambdaParameter.getParameterValue(0));
        currentQ.mult(currentGamma, diagonal1);

        currentLike += 0.5 * (fieldLength - 1) * Math.log(precisionParameter.getParameterValue(0)) - 0.5 * currentGamma.dot(diagonal1);
        if (lambdaParameter.getParameterValue(0) == 1) {
            currentLike -= (fieldLength - 1) / 2.0 * LOG_TWO_TIMES_PI;
        } else {
            currentLike -= fieldLength / 2.0 * LOG_TWO_TIMES_PI;
        }

        return currentLike;
    }

    // TODO: Potentially should be in a BNPRHelper class
    public double calculateLogFieldLikelihoodSubGamma(double[] gamma) {

        if (!intervalsKnown) {
            //intervalsKnown -> false when handleModelChanged event occurs in super.
            wrapSetupIntervals();
            setupSufficientStatistics();
            intervalsKnown = true;
        }

        DenseVector diagonal1 = new DenseVector(fieldLength);
        DenseVector currentGamma = new DenseVector(gamma);

        //updateGammaWithCovariates(currentGamma);

        //double currentLike = handleMissingValues();
        double currentLike = 0.0;

        SymmTridiagMatrix currentQ = getScaledWeightMatrix(precisionParameter.getParameterValue(0), lambdaParameter.getParameterValue(0));
        currentQ.mult(currentGamma, diagonal1);

        currentLike += 0.5 * (fieldLength - 1) * Math.log(precisionParameter.getParameterValue(0)) - 0.5 * currentGamma.dot(diagonal1);
        if (lambdaParameter.getParameterValue(0) == 1) {
            currentLike -= (fieldLength - 1) / 2.0 * LOG_TWO_TIMES_PI;
        } else {
            currentLike -= fieldLength / 2.0 * LOG_TWO_TIMES_PI;
        }

        return currentLike;
    }

    protected void setupSufficientStatistics() {
        int index = 0;

        double length = 0;
        double weight = 0;
        for (int i = 0; i < getIntervalCount(); i++) {
            length += getInterval(i);
            weight += getInterval(i) * getLineageCount(i) * (getLineageCount(i) - 1);
            if (getIntervalType(i) == CoalescentEventType.COALESCENT) {
                coalescentIntervals[index] = length;
                sufficientStatistics[index] = weight / 2.0;
                index++;
                length = 0;
                weight = 0;
            }
        }
    }

    protected void setupGPWeights() { // TODO: Switch to Gaussian Process weights

        setupSufficientStatistics();

        //Set up the weight Matrix
        double[] offdiag = new double[fieldLength - 1];
        double[] diag = new double[fieldLength];

        //First set up the offdiagonal entries;

        for (int i = 0; i < fieldLength - 1; i++) {
            offdiag[i] = -1.0 / (midPoints[i+1] - midPoints[i]);
        }

        //Then set up the diagonal entries;
        for (int i = 1; i < fieldLength - 1; i++)
            diag[i] = -(offdiag[i] + offdiag[i - 1]);

        //Take care of the endpoints
        diag[0] = -offdiag[0] + NUGGET;
        diag[fieldLength - 1] = -offdiag[fieldLength - 2];

        weightMatrix = new SymmTridiagMatrix(diag, offdiag);
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

    // Main function for testing

    public static void main(String[] args) throws Exception {
        NewickImporter importer = new NewickImporter("((((5:0.5,1:0.2):0.5,0:1):0.2,2:0.8):0.2,3:1.4)");
        Tree tree = importer.importNextTree();
        List<Tree> treeList = wrapTree(tree);

        System.out.println("Got here A");

        BNPRLikelihood bnprLikelihood = new BNPRLikelihood(treeList, new Parameter.Default(1.0),
                null, new Parameter.Default(1.0),
                new Parameter.Default(1.0), new Parameter.Default(1.0), null,
                true, true, 1.4, 5, false);

        System.out.println("Got here B");

        double ll = bnprLikelihood.calculateLogCoalescentLikelihood();
        System.out.println("Log-likelihood: " + ll);

        System.out.println("Grid: " + Arrays.toString(bnprLikelihood.gridPoints));
        System.out.println("Midpoints: " + Arrays.toString(bnprLikelihood.midPoints));

        System.out.println("Got here C");

//        System.out.println("Q matrix:");
//
//        SymmTridiagMatrix Qmat = bnprLikelihood.getScaledWeightMatrix(1);
//        System.out.println(Qmat);
//
//        System.out.println("Got here D");
//
//        BandCholesky chol = BandCholesky.factorize(new UpperSPDBandMatrix(Qmat, 1));
//        DenseVector randNorm = BNPRBlockUpdateOperator.getMultiNormal(new DenseVector(bnprLikelihood.midPoints.length), chol);
//
//        System.out.println("Randnorm: " + Arrays.toString(randNorm.getData()));

        System.out.println("Got here E");

        int intervalCount = bnprLikelihood.getIntervalCount();
        System.out.println("Interval count: " + intervalCount);
        for (int i = 0; i < intervalCount; i++) {
            System.out.printf("Num lineages %d: %d\n", i, bnprLikelihood.getLineageCount(i));
            System.out.printf("Interval %d: %.2f\n", i, bnprLikelihood.getInterval(i));
        }

        System.out.println("Got here F");

        System.out.println("Tree root: " + bnprLikelihood.tree.getNodeHeight(bnprLikelihood.tree.getRoot()));

    }

    private static List<Tree> wrapTree(Tree tree) {
        List<Tree> treeList = new ArrayList<Tree>();
        treeList.add(tree);
        return treeList;
    }

    public class GriddedTreeIntervals extends TreeIntervals {
        protected double[] grid;
        protected IntervalType[] eventTypes;

        public GriddedTreeIntervals() {
            super();
        }

        public GriddedTreeIntervals(Tree tree) {
            super(tree);
        }

        public void setGrid(double[] grid) {
            this.grid = grid.clone();
            // TODO: Propagate new grid to intervals, event types, etc.
        }
    }
}
