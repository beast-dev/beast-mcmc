package dr.evomodel.coalescent;

import dr.evolution.coalescent.IntervalList;
import dr.evolution.coalescent.IntervalType;
import dr.evolution.io.NewickImporter;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeUtils;
import dr.evolution.util.Units;
import dr.evomodel.tree.TreeModel;
import dr.evomodelxml.coalescent.BNPRLikelihoodParser;
import dr.inference.model.Parameter;
import dr.util.Author;
import dr.util.Citable;
import dr.util.Citation;
import dr.util.HeapSort;
import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.SymmTridiagMatrix;

import java.util.*;

/**
 * Created by mkarcher on 9/15/16.
 */
public class BNPRLikelihood extends AbstractCoalescentLikelihood implements CoalescentIntervalProvider, Citable {
    public static final double NUGGET = 0.0001;

    protected Parameter popSizeParameter;
    protected Parameter groupSizeParameter;
    protected Parameter precisionParameter;

    protected Tree tree;

    protected int fieldLength;
    protected double[] coalescentIntervals;
    protected double[] storedCoalescentIntervals;
    protected double[] sufficientStatistics;
    protected double[] storedSufficientStatistics;

    protected double cutOff;
    protected int numGridPoints;
    protected int oldFieldLength;
    protected double[] numCoalEvents;
    protected double[] storedNumCoalEvents;
    protected double[] gridPoints;
    protected double[] midPoints;
    protected GriddedTreeIntervals gti;

    protected double precAlpha;
    protected double precBeta;

    protected Parameter samplingBetas;

    protected double logSamplingLikelihood;
    protected double storedLogSamplingLikelihood;

    protected double logFieldLikelihood;
    protected double storedLogFieldLikelihood;

    protected SymmTridiagMatrix weightMatrix;
    protected SymmTridiagMatrix storedWeightMatrix;

    protected boolean intervalsKnown;
    protected boolean samplingAware;
    protected double[] samplingTimes;
    protected boolean samplingTimesKnown;

    public BNPRLikelihood(Tree tree,
                          Parameter popParameter,
                          Parameter groupParameter,
                          Parameter precParameter,
                          double cutOff,
                          int numGridPoints,
                          boolean samplingAware) throws TreeUtils.MissingTaxonException {

        super(BNPRLikelihoodParser.BNPR_LIKELIHOOD, tree, null, null);

        this.popSizeParameter = popParameter;
        this.groupSizeParameter = groupParameter;
        this.precisionParameter = precParameter;

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
        //if (betaParameter != null) {
        //    addVariable(betaParameter);
        //}

        setTree(tree);

        oldFieldLength = getOldFieldLength(); // TODO: Refactor

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

    protected int getOldFieldLength() {
        return tree.getExternalNodeCount() - 1;
    }

    protected void setTree(Tree tree) {
        this.tree = tree;
    }

    public void initializationReport() {
        System.out.println("Creating a BNPR model");
        System.out.println("\tPopulation sizes: " + popSizeParameter.getDimension());
    }

    public void wrapSetupIntervals() {
//        setupIntervals();
        this.gti = new GriddedTreeIntervals(this.tree, this.gridPoints);
        this.gti.calculateIntervals();
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
            samplingTimes = new double[n];

            for (int i = 0; i < n; i++) {
                samplingTimes[i] = tree.getNodeHeight(tree.getExternalNode(i));
            }

            HeapSort.sort(samplingTimes);
            samplingTimesKnown = true;
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


    protected static double[] repeatGamma(double[] gamma, int[] indices) {
        double[] result = new double[indices.length];

        for (int i = 0; i < indices.length; i++) {
            result[i] = gamma[indices[i]];
        }

        return result;
    }

    //***********************
    // Calculate Likelihood
    //***********************

//    public double getLogLikelihood() {
//        if (!likelihoodKnown) {
//            logLikelihood = calculateLogCoalescentLikelihood();
//            logFieldLikelihood = calculateLogFieldLikelihood();
//
//            if (samplingAware)
//                logSamplingLikelihood = calculateLogSamplingLikelihood();
//
//            likelihoodKnown = true;
//        }
//
//        double ll = logLikelihood + logFieldLikelihood;
//        if (samplingAware)
//            ll += logSamplingLikelihood;
//
//        return ll;
//    }

    @Override
    public double calculateLogLikelihood() {
        return 0;
    }

    public double getLogLikelihoodSubGamma(double[] gamma) {
        double logLikelihood = calculateLogCoalescentLikelihoodSubGamma(gamma);
        double logFieldLikelihood = calculateLogFieldLikelihoodSubGamma(gamma);
        double logSamplingLikelihood = 0;

        if (samplingAware) logSamplingLikelihood = calculateLogSamplingLikelihoodSubGamma(gamma);

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
        double[] repeatedGamma = repeatGamma(currentGamma, gti.getIndices());

        for (int i = 0; i < gti.getIntervalCount(); i++) {
//            System.out.println("Interval "+i);
            double intervalLength = gti.getInterval(i);
//            System.out.println("Interval length: " + intervalLength);
            IntervalTypeN intervalType = gti.getIntervalTypeN(i);
//            System.out.println("Interval type : " + intervalType.type() + " " + intervalType.N());
            int lineages = gti.getLineageCount(i);
//            System.out.println("Number of lineages: " + lineages);
//            System.out.println("Coalescent pressure: " + lineages * (lineages-1) / 2);

            currentLike += -intervalLength * (lineages * (lineages-1) / 2) * Math.exp(-repeatedGamma[i]);
            if (intervalType.type() == IntervalType.COALESCENT) {
                currentLike += -repeatedGamma[i] * intervalType.N();
            }
        }

        return currentLike;// + LogNormalDistribution.logPdf(Math.exp(popSizeParameter.getParameterValue(coalescentIntervals.length - 1)), mu, sigma);
    }


    public double calculateLogCoalescentLikelihoodSubGamma(double[] gamma) {

        double currentLike = 0;
        double[] currentGamma = gamma;
        double[] repeatedGamma = repeatGamma(currentGamma, gti.getIndices());

        for (int i = 0; i < gti.getIntervalCount(); i++) {
//            System.out.println("Interval "+i);
            double intervalLength = gti.getInterval(i);
//            System.out.println("Interval length: " + intervalLength);
            IntervalTypeN intervalType = gti.getIntervalTypeN(i);
//            System.out.println("Interval type : " + intervalType.type() + " " + intervalType.N());
            int lineages = gti.getLineageCount(i);
//            System.out.println("Number of lineages: " + lineages);
//            System.out.println("Coalescent pressure: " + lineages * (lineages-1) / 2);

            currentLike += -intervalLength * (lineages * (lineages-1) / 2) * Math.exp(-repeatedGamma[i]);
            if (intervalType.type() == IntervalType.COALESCENT) {
                currentLike += -repeatedGamma[i] * intervalType.N();
            }
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

    public double calculateLogFieldLikelihood() { // TODO: Convert to GP prior

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

        SymmTridiagMatrix currentQ = getScaledWeightMatrix(precisionParameter.getParameterValue(0));
        currentQ.mult(currentGamma, diagonal1);

        currentLike += 0.5 * (fieldLength - 1) * Math.log(precisionParameter.getParameterValue(0)) - 0.5 * currentGamma.dot(diagonal1);
//        if (lambdaParameter.getParameterValue(0) == 1) {
//            currentLike -= (fieldLength - 1) / 2.0 * LOG_TWO_TIMES_PI;
//        } else {
//            currentLike -= fieldLength / 2.0 * LOG_TWO_TIMES_PI;
//        }

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

        SymmTridiagMatrix currentQ = getScaledWeightMatrix(precisionParameter.getParameterValue(0));
        currentQ.mult(currentGamma, diagonal1);

        currentLike += 0.5 * (fieldLength - 1) * Math.log(precisionParameter.getParameterValue(0)) - 0.5 * currentGamma.dot(diagonal1);
//        if (lambdaParameter.getParameterValue(0) == 1) {
//            currentLike -= (fieldLength - 1) / 2.0 * LOG_TWO_TIMES_PI;
//        } else {
//            currentLike -= fieldLength / 2.0 * LOG_TWO_TIMES_PI;
//        }

        return currentLike;
    }

    public Parameter getPrecisionParameter() {
        return precisionParameter;
    }

    public Parameter getPopSizeParameter() {
        return popSizeParameter;
    }

    protected void setupSufficientStatistics() {
        int index = 0;

        double length = 0;
        double weight = 0;
        for (int i = 0; i < gti.getIntervalCount(); i++) {
            length += gti.getInterval(i);
            weight += gti.getInterval(i) * gti.getLineageCount(i) * (gti.getLineageCount(i) - 1);
            if (gti.getIntervalType(i) == IntervalType.COALESCENT) {
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

    public SymmTridiagMatrix getScaledWeightMatrix(double precision) {
        SymmTridiagMatrix a = weightMatrix.copy();
        for (int i = 0; i < a.numRows() - 1; i++) {
            a.set(i, i, a.get(i, i) * precision);
            a.set(i + 1, i, a.get(i + 1, i) * precision);
        }
        a.set(fieldLength - 1, fieldLength - 1, a.get(fieldLength - 1, fieldLength - 1) * precision);
        return a;
    }

    public SymmTridiagMatrix getStoredScaledWeightMatrix(double precision) {
        SymmTridiagMatrix a = storedWeightMatrix.copy();
        for (int i = 0; i < a.numRows() - 1; i++) {
            a.set(i, i, a.get(i, i) * precision);
            a.set(i + 1, i, a.get(i + 1, i) * precision);
        }
        a.set(fieldLength - 1, fieldLength - 1, a.get(fieldLength - 1, fieldLength - 1) * precision);
        return a;
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

    // From CoalescentIntervalProvider
    public int getNumberOfCoalescentEvents() {
        return tree.getExternalNodeCount() - 1;
    }

    // From CoalescentIntervalProvider
    public double getCoalescentEventsStatisticValue(int i) {
        return sufficientStatistics[i];
    }

    @Override
    public Citation.Category getCategory() {
        return Citation.Category.TREE_PRIORS;
    }

    @Override
    public String getDescription() {
        return "Bayesian nonparametric phylodynamic reconstruction coalescent";
    }

    @Override
    public List<Citation> getCitations() {
        return Collections.singletonList(CITATION);
    }

    public static Citation CITATION = new Citation(
            new Author[]{
                    new Author("MD", "Karcher"),
                    new Author("JA", "Palacios"),
                    new Author("T", "Bedford"),
                    new Author("MA", "Suchard"),
                    new Author("VN", "Minin")
            },
            "Quantifying and mitigating the effect of preferential sampling on phylodynamic inference",
            2016,
            "PLoS Comput Biol",
            12, "e1004789",
            "10.1371/journal.pcbi.1004789"
    );

    // Main function for testing

    public static void main(String[] args) throws Exception {
        NewickImporter importer = new NewickImporter("((((5:0.5,1:0.2):0.5,0:1):0.2,2:0.8):0.2,3:1.4)");
        Tree tree = importer.importNextTree();
//        List<Tree> treeList = wrapTree(tree);

        System.out.println("Got here A");

        BNPRLikelihood bnprLikelihood = new BNPRLikelihood(tree, new Parameter.Default(new double[] {10.0, 4.7, 2.1, 4.8, 1.9, 1.59}),
                null, new Parameter.Default(1.0),
                1.4, 5, false);
//        BNPRLikelihood bnprLikelihood = new BNPRLikelihood(treeList, new Parameter.Default(1.0),
//                null, new Parameter.Default(1.0),
//                new Parameter.Default(1.0), new Parameter.Default(1.0), null,
//                true, true, 1.4, 5, false);


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

//        System.out.println("Got here E");
//
//        int intervalCount = bnprLikelihood.getIntervalCount();
//        System.out.println("Interval count: " + intervalCount);
//        for (int i = 0; i < intervalCount; i++) {
//            System.out.printf("Num lineages %d: %d\n", i, bnprLikelihood.getLineageCount(i));
//            System.out.printf("Interval %d: %.2f\n", i, bnprLikelihood.getInterval(i));
//        }
//
//        System.out.println("Got here F");
//
//        System.out.println("Tree root: " + bnprLikelihood.tree.getNodeHeight(bnprLikelihood.tree.getRoot()));
//
//        System.out.println("Got here G");

        GriddedTreeIntervals gti = new GriddedTreeIntervals(bnprLikelihood.tree, bnprLikelihood.gridPoints);
//        System.out.println("eventMap: " + gti.eventMap.keySet().toString());
        System.out.println("Got here H");

        System.out.println("Interval count: " + gti.getIntervalCount());
        System.out.println("Sample count: " + gti.getSampleCount());
        System.out.println("First interval length: " + gti.getInterval(0));
        System.out.println("Got here I");

        System.out.println("Lineages during first interval: " + gti.getLineageCount(0));
        System.out.println("Got here J");

        System.out.println("Last interval length: " + gti.getInterval(gti.getIntervalCount()-1));
        System.out.println("Got here K");

        System.out.println("toString: " + gti.toString());

        System.out.println("Got here L");

        System.out.println("Sampling times: " + Arrays.toString(bnprLikelihood.getSamplingTimes()));

        System.out.println("Got here M");

        System.out.println("Repetition index: " + Arrays.toString(gti.getIndices()));
        System.out.println("Repeated gamma: " + Arrays.toString(repeatGamma(bnprLikelihood.getPopSizeParameter().getParameterValues(), gti.getIndices())));
    }

    private static List<Tree> wrapTree(Tree tree) {
        List<Tree> treeList = new ArrayList<Tree>();
        treeList.add(tree);
        return treeList;
    }

    @Override
    public Type getUnits() {
        return null; // TODO: Implement
    }

    @Override
    public void setUnits(Type units) {
        // TODO: Implement
    }

    public static class GriddedTreeIntervals implements IntervalList {
        protected Tree tree = null;
        protected double[] grid;
        protected TreeMap<Double, IntervalTypeN> eventMap;
        protected Set<Double> times;
        protected double[] nodeTimes;
        protected int[] childCounts;

        protected boolean intervalsKnown = false;
        protected boolean storedIntervalsKnown;

        protected double multifurcationLimit = -1.0;
        protected double simultaneousSamplingLimit = 1e-9;

        public GriddedTreeIntervals() {

        }

        public GriddedTreeIntervals(Tree tree, double[] grid) {
            setTree(tree);
            setGrid(grid);
        }

        /**
         * Sets the tree for which intervals are obtained
         */
        public void setTree(Tree tree) {
            this.tree = tree;
            intervalsKnown = false;
        }

        public Tree getTree() {
            return this.tree;
        }

        public void setGrid(double[] grid) {
            this.grid = grid.clone();
            intervalsKnown = false;
        }

        public double[] getGrid() {
            return this.grid;
        }

        protected void calculateIntervals() {
            eventMap = new TreeMap<Double, IntervalTypeN>();
            this.times = eventMap.keySet();

            int n = this.tree.getNodeCount();
            this.nodeTimes = new double[n];
            this.childCounts = new int[n];
            collectTimes(this.tree, this.nodeTimes, this.childCounts);

            for (int i = 0; i < grid.length; i++) {
                eventMap.put(new Double(grid[i]), new IntervalTypeN(IntervalType.NOTHING));
            }

            for (int i = 0; i < n; i++) {
//                System.out.printf("Node %d has %d children at time %f\n", i, childCounts[i], nodeTimes[i]);
                Double currentTime = new Double(nodeTimes[i]);
                Map.Entry<Double, IntervalTypeN> floorEntry = eventMap.floorEntry(currentTime);

                if (childCounts[i] == 0) {
                    if (floorEntry.getValue().type() == IntervalType.SAMPLE && Math.abs(nodeTimes[i] - floorEntry.getKey().doubleValue()) <= this.simultaneousSamplingLimit) {
                        int N = floorEntry.getValue().N();
                        eventMap.put(currentTime, new IntervalTypeN(IntervalType.SAMPLE, N + 1));
                    } else {
                        eventMap.put(currentTime, new IntervalTypeN(IntervalType.SAMPLE));
                    }
                } else {
                    if (floorEntry.getValue().type() == IntervalType.COALESCENT && Math.abs(nodeTimes[i] - floorEntry.getKey().doubleValue()) <= this.multifurcationLimit) {
                        int N = floorEntry.getValue().N();
                        eventMap.put(currentTime, new IntervalTypeN(IntervalType.COALESCENT, N + 1));
                    } else {
                        eventMap.put(currentTime, new IntervalTypeN(IntervalType.COALESCENT));
                    }
                }
            }

            intervalsKnown = true;
        }

        public int[] getIndices() {
            int[] result = new int[getIntervalCount()];
            int currentPlace = 0;

            for (int i = 1; i < this.grid.length; i++) {
                int n = this.eventMap.subMap(grid[i-1], grid[i]).size();

                for (int j = 0; j < n; j++) {
                    result[currentPlace++] = i;
                }
            }

            return result;
        }

        /**
         * Specifies that the intervals are unknown (e.g., the tree has changed).
         */
        public void setIntervalsUnknown() {
            intervalsKnown = false;
        }

        /**
         * Sets the limit for which adjacent events are merged.
         *
         * @param multifurcationLimit A value of 0 means merge addition of leafs (terminal nodes) when possible but
         *                            return each coalescense as a separate event.
         */
        public void setMultifurcationLimit(double multifurcationLimit) {
            this.multifurcationLimit = multifurcationLimit;
            intervalsKnown = false;
        }

        public void setSimultaneousSamplingLimit(double simultaneousSamplingLimit) {
            this.simultaneousSamplingLimit = simultaneousSamplingLimit;
            intervalsKnown = false;
        }

        public int getSampleCount() {
            return tree.getExternalNodeCount();
        }

        public int getIntervalCount() {
            if (!intervalsKnown) {
                calculateIntervals();
            }
            return this.times.size() - 1;
        }

        public double getInterval(int i) {
            if (!intervalsKnown) {
                calculateIntervals();
            }
            if (i >= getIntervalCount()) {
                throw new IllegalArgumentException();
            }
            Double[] timesArray = new Double[times.size()];
            times.toArray(timesArray);
            return timesArray[i+1] - timesArray[i];
        }

        /**
         * Returns the number of uncoalesced lineages within this interval.
         * Required for s-coalescents, where new lineages are added as
         * earlier samples are come across.
         */
        public int getLineageCount(int i) {
            if (!intervalsKnown) {
                calculateIntervals();
            }
            if (i >= getIntervalCount()) throw new IllegalArgumentException();

            int lineages = 0;
            Double time = (Double) this.times.toArray()[i+1];
            NavigableMap<Double, IntervalTypeN> submap = this.eventMap.headMap(time, false);
//            System.out.println("getLineageCount: number of times less than " + time + ": " + submap.size());

            for (Double timeI : submap.keySet()) {
//                System.out.println("getLineageCount: loop: time: " + timeI + ", N: " + submap.get(timeI).N());
                if (submap.get(timeI).type() == IntervalType.SAMPLE) {
                    lineages += submap.get(timeI).N();
                } else if (submap.get(timeI).type() == IntervalType.COALESCENT) {
                    lineages -= submap.get(timeI).N();
                }
            }

            return lineages;
        }

        /**
         * Returns the number coalescent events in an interval
         */
        public int getCoalescentEvents(int i) {
            if (!intervalsKnown) {
                calculateIntervals();
            }
            if (i >= getIntervalCount()) throw new IllegalArgumentException();
            Double time = (Double) this.times.toArray()[i];
            IntervalTypeN typeN = this.eventMap.get(time);
            int result = 0;

            if (typeN.type() == IntervalType.COALESCENT) {
                result += typeN.N();
            }

            return result;
        }

        /**
         * Returns the type of interval observed.
         */
        public IntervalType getIntervalType(int i) {
            if (!intervalsKnown) {
                calculateIntervals();
            }
            if (i >= getIntervalCount()) throw new IllegalArgumentException();
            Double time = (Double) this.times.toArray()[i+1];
            IntervalTypeN typeN = this.eventMap.get(time);
            return typeN.type();
        }

        /**
         * Returns the number of events at the end of the interval
         */
        public int getIntervalN(int i) {
            if (!intervalsKnown) {
                calculateIntervals();
            }
            if (i >= getIntervalCount()) throw new IllegalArgumentException();
            Double time = (Double) this.times.toArray()[i+1];
            IntervalTypeN typeN = this.eventMap.get(time);
            return typeN.N();
        }

        /**
         * Returns the type of interval observed.
         */
        public IntervalTypeN getIntervalTypeN(int i) {
            if (!intervalsKnown) {
                calculateIntervals();
            }
            if (i >= getIntervalCount()) throw new IllegalArgumentException();
            Double time = (Double) this.times.toArray()[i+1];
            IntervalTypeN typeN = this.eventMap.get(time);
            return typeN;
        }

        /**
         * get the total height of the genealogy represented by these
         * intervals.
         */
        public double getTotalDuration() {

            if (!intervalsKnown) {
                calculateIntervals();
            }
            return this.eventMap.lastKey().doubleValue();
        }

        /**
         * Checks whether this set of coalescent intervals is fully resolved
         * (i.e. whether is has exactly one coalescent event in each
         * subsequent interval)
         */
        public boolean isBinaryCoalescent() {
            if (!intervalsKnown) {
                calculateIntervals();
            }
            for (int i = 0; i < getIntervalCount(); i++) {
                if (getCoalescentEvents(i) > 0) {
                    if (getCoalescentEvents(i) != 1) return false;
                }
            }

            return true;
        }

        /**
         * Checks whether this set of coalescent intervals coalescent only
         * (i.e. whether is has exactly one or more coalescent event in each
         * subsequent interval)
         */
        public boolean isCoalescentOnly() {
            if (!intervalsKnown) {
                calculateIntervals();
            }
            for (int i = 0; i < getIntervalCount(); i++) {
                if (getCoalescentEvents(i) < 1) return false;
            }

            return true;
        }

        /**
         * Returns the time of the start of an interval
         *
         * @param i which interval
         * @return start time
         */
        public double getIntervalTime(int i) {
            if (!intervalsKnown) {
                calculateIntervals();
            }
            if (i >= getIntervalCount()) throw new IllegalArgumentException();
            return ((Double) this.times.toArray()[i+1]).doubleValue();
        }

        /**
         * extract coalescent times and tip information into array times from tree.
         */
        protected static void collectTimes(Tree tree, double[] times, int[] childCounts) {

            for (int i = 0; i < tree.getNodeCount(); i++) {
                NodeRef node = tree.getNode(i);
                times[i] = tree.getNodeHeight(node);
                childCounts[i] = tree.getChildCount(node);
            }
        }

        /**
         * Return the units that this tree is expressed in.
         */
        public final Type getUnits() {
            return tree.getUnits();
        }

        /**
         * Sets the units that this tree is expressed in.
         */
        public final void setUnits(Type units) {
            throw new IllegalArgumentException("Can't set interval's units");
        }

        /**
         * Extra functionality to store and restore values for caching
         */
        public void storeState() { // TODO: Implement
            if (intervalsKnown) {


            }

            storedIntervalsKnown = intervalsKnown;
        }

        public void restoreState() { // TODO: Implement
            intervalsKnown = storedIntervalsKnown;

            if (intervalsKnown) {

            }
        }

        public String toString() {
            StringBuilder sb = new StringBuilder();
            for (int i = 0; i < getIntervalCount(); i++) {
                sb.append("[ ");
                sb.append(getInterval(i));
                sb.append(": ");
                sb.append(getIntervalTime(i));
                sb.append(": ");
                sb.append(getLineageCount(i));
                sb.append(" ]");
            }
            return sb.toString();
        }
    }

    public static class IntervalTypeN {
        private final IntervalType type;
        private final int N;

        public IntervalTypeN(IntervalType type) {
            this.type = type;
            this.N = 1;
        }

        public IntervalTypeN(IntervalType type, int N) {
            this.type = type;
            this.N = N;
        }

        public IntervalType type() {
            return this.type;
        }

        public int N() {
            return this.N;
        }
    }
}
