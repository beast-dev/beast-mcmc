package dr.app.tools;

import dr.evolution.alignment.SimpleAlignment;
import dr.evolution.alignment.SitePatterns;
import dr.evolution.datatype.DataType;
import dr.evolution.datatype.GeneralDataType;
import dr.evolution.sequence.Sequence;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.SimpleNode;
import dr.evolution.tree.SimpleTree;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxon;
import dr.evomodel.branchmodel.RewardsAwareBranchModel;
import dr.evomodel.branchmodel.RewardsAwareBranchModelRewardRateGradient;
import dr.evomodel.branchmodel.RewardsAwareCtmcSimulation;
import dr.evomodel.branchratemodel.ArbitraryBranchRates;
import dr.evomodel.branchratemodel.DefaultBranchRateModel;
import dr.evomodel.branchratemodel.RewardRates;
import dr.evomodel.branchratemodel.RewardsAwareMixtureBranchRates;
import dr.evomodel.siteratemodel.GammaSiteRateModel;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.LogAdditiveCtmcRateProvider;
import dr.evomodel.substmodel.LogRateSubstitutionModel;
import dr.evomodel.tree.DefaultTreeModel;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.tree.TreeParameterModel;
import dr.evomodel.treedatalikelihood.DiscreteDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.PreOrderSettings;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.discrete.discretetreedataLikelihood.representations.RewardsAwarePartialsRepresentation;
import dr.evomodel.treelikelihood.PartialsRescalingScheme;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.CompoundParameter;
import dr.inference.model.Parameter;
import dr.inference.operators.RewardsMixtureIndicatorAndAtomIndicesOperator;
import dr.math.MathUtils;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Independent-likelihood simulation and operator validation driver for the
 * rewards-aware CTMC machinery.
 *
 * Usage:
 * java -cp beast.jar dr.app.tools.RewardsAwareIndependentOperatorValidation \
 *   outputDir tipCount stateCount scenario sweeps burnIn logEvery seed
 */
public final class RewardsAwareIndependentOperatorValidation {

    private static final double BRANCH_HEIGHT_UNIT = 0.18;
    private static final double REWARD_RATE_HMC_STEP_SIZE =
            getDoubleProperty("rewardValidation.rewardRateHmcStepSize", 0.01);
    private static final int REWARD_RATE_HMC_LEAPFROG_STEPS =
            Integer.getInteger("rewardValidation.rewardRateHmcLeapfrogSteps", 5);
    private static final int REWARD_RATE_HMC_STEPS =
            Integer.getInteger("rewardValidation.rewardRateHmcSteps", 1);
    private static final double CTS_RANDOM_WALK_SCALE =
            getDoubleProperty("rewardValidation.ctsRandomWalkScale", 0.05);
    private static final int CTS_MOVES_PER_SWEEP =
            Integer.getInteger("rewardValidation.ctsMovesPerSweep", 0);

    private RewardsAwareIndependentOperatorValidation() {
        // executable utility
    }

    public static void main(final String[] args) throws Exception {
        if (args.length != 8) {
            throw new IllegalArgumentException(
                    "Expected 8 arguments: outputDir tipCount stateCount scenario sweeps burnIn logEvery seed");
        }

        final File outputDir = new File(args[0]);
        final int tipCount = Integer.parseInt(args[1]);
        final int stateCount = Integer.parseInt(args[2]);
        final Scenario scenario = Scenario.fromId(Integer.parseInt(args[3]));
        final int sweeps = Integer.parseInt(args[4]);
        final int burnIn = Integer.parseInt(args[5]);
        final int logEvery = Integer.parseInt(args[6]);
        final long seed = Long.parseLong(args[7]);

        validateArguments(outputDir, tipCount, stateCount, sweeps, burnIn, logEvery);
        MathUtils.setSeed(seed);

        final ValidationRun run = new ValidationRun(outputDir, tipCount, stateCount, scenario,
                sweeps, burnIn, logEvery, seed);
        run.execute();
    }

    private static void validateArguments(final File outputDir,
                                          final int tipCount,
                                          final int stateCount,
                                          final int sweeps,
                                          final int burnIn,
                                          final int logEvery) throws IOException {
        if (tipCount < 2) {
            throw new IllegalArgumentException("tipCount must be at least 2: " + tipCount);
        }
        if (stateCount < 3) {
            throw new IllegalArgumentException("stateCount must be at least 3: " + stateCount);
        }
        if (sweeps <= 0) {
            throw new IllegalArgumentException("sweeps must be positive: " + sweeps);
        }
        if (burnIn < 0 || burnIn >= sweeps) {
            throw new IllegalArgumentException("burnIn must be in [0, sweeps): " + burnIn);
        }
        if (logEvery <= 0) {
            throw new IllegalArgumentException("logEvery must be positive: " + logEvery);
        }
        if (!outputDir.exists() && !outputDir.mkdirs()) {
            throw new IOException("Could not create output directory: " + outputDir);
        }
    }

    private static double getDoubleProperty(final String name, final double defaultValue) {
        final String raw = System.getProperty(name);
        return raw == null ? defaultValue : Double.parseDouble(raw);
    }

    private enum Scenario {
        REWARD_VALUES(1, "reward_values", true, false, false, false),
        REWARD_ORDER(2, "reward_order", false, true, false, false),
        ALL_REWARD_RATES(3, "all_reward_rates", true, true, false, false),
        TOTAL_REWARD_CTS(4, "total_reward_cts", false, false, true, false),
        TOTAL_REWARD_INDICATOR(5, "total_reward_indicator", false, false, false, true),
        TOTAL_REWARD_CTS_INDICATOR(6, "total_reward_cts_indicator", false, false, true, true),
        EVERYTHING(7, "everything", true, true, true, true);

        final int id;
        final String label;
        final boolean updateRewardValues;
        final boolean updateRewardOrder;
        final boolean updateCtsRewards;
        final boolean updateIndicators;

        Scenario(final int id,
                 final String label,
                 final boolean updateRewardValues,
                 final boolean updateRewardOrder,
                 final boolean updateCtsRewards,
                 final boolean updateIndicators) {
            this.id = id;
            this.label = label;
            this.updateRewardValues = updateRewardValues;
            this.updateRewardOrder = updateRewardOrder;
            this.updateCtsRewards = updateCtsRewards;
            this.updateIndicators = updateIndicators;
        }

        static Scenario fromId(final int id) {
            for (Scenario scenario : values()) {
                if (scenario.id == id) {
                    return scenario;
                }
            }
            throw new IllegalArgumentException("Unknown scenario id: " + id);
        }
    }

    private static final class ValidationRun {

        private final File outputDir;
        private final int tipCount;
        private final int stateCount;
        private final Scenario scenario;
        private final int sweeps;
        private final int burnIn;
        private final int logEvery;
        private final long seed;

        private final TreeModel tree;
        private final DataType dataType;
        private final double[] rootFrequencies;
        private final double[] truthRewardValues;
        private final int[] truthMapping;
        private final double[] truthRewardByState;
        private final RewardsAwareCtmcSimulation simulation;
        private final double[] trueCts;
        private final int[] trueIndicators;
        private final int[] trueAtoms;

        private final Parameter ctsRewards;
        private final Parameter indicator;
        private final Parameter atomIndices;
        private final RewardRates rewardRates;
        private final RewardsAwareBranchModel rewardsAwareBranchModel;
        private final TreeDataLikelihood independentLikelihood;

        private final RewardsMixtureIndicatorAndAtomIndicesOperator indicatorOperator;
        private final RewardRateHmcKernel rewardRateHmc;
        private final CtsRewardRandomWalk ctsRandomWalk;
        private final RewardOrderOperator rewardOrderOperator;

        private int sampleCount = 0;

        ValidationRun(final File outputDir,
                      final int tipCount,
                      final int stateCount,
                      final Scenario scenario,
                      final int sweeps,
                      final int burnIn,
                      final int logEvery,
                      final long seed) {
            this.outputDir = outputDir;
            this.tipCount = tipCount;
            this.stateCount = stateCount;
            this.scenario = scenario;
            this.sweeps = sweeps;
            this.burnIn = burnIn;
            this.logEvery = logEvery;
            this.seed = seed;

            this.tree = createBalancedTree(tipCount);
            this.dataType = createDataType(stateCount);
            this.rootFrequencies = uniform(stateCount);
            this.truthRewardValues = createTruthRewardValues(stateCount);
            this.truthMapping = createTruthMapping(stateCount);
            this.truthRewardByState = valuesByState(truthRewardValues, truthMapping);

            final LogRateSubstitutionModel truthModel =
                    createIndependentRewardSubstitutionModel("truth.independent.model", dataType,
                            rootFrequencies, stateCount);
            final double[] qMatrix = new double[stateCount * stateCount];
            truthModel.getInfinitesimalMatrix(qMatrix);

            this.simulation = RewardsAwareCtmcSimulation.withRootFrequencies(
                    "independent.operator.validation.simulation",
                    stateCount,
                    truthRewardByState,
                    tree,
                    qMatrix,
                    rootFrequencies
            );

            this.trueCts = simulation.getBranchRewardProportions();
            this.trueIndicators = simulation.getIndicators();
            this.trueAtoms = simulation.getAtomIndices();

            this.rewardRates = createInferenceRewardRates(
                    "inference.rewardRates",
                    truthRewardValues,
                    truthMapping,
                    scenario.updateRewardValues,
                    scenario.updateRewardOrder
            );
            this.ctsRewards = new Parameter.Default("inference.ctsRewards",
                    initializeCtsRewards(trueCts, scenario.updateCtsRewards));
            this.indicator = new Parameter.Default("inference.indicator",
                    initializeIndicators(trueIndicators, scenario.updateIndicators));
            this.atomIndices = new Parameter.Default("inference.atomIndices",
                    initializeAtoms(trueAtoms, scenario.updateIndicators, stateCount));

            final RewardsAwareMixtureBranchRates branchRates = new RewardsAwareMixtureBranchRates(
                    tree,
                    ctsRewards,
                    indicator,
                    atomIndices,
                    rewardRates,
                    new ArbitraryBranchRates.BranchRateTransform.None(),
                    false,
                    TreeParameterModel.Type.WITHOUT_ROOT
            );

            final LogRateSubstitutionModel inferenceModel =
                    createIndependentRewardSubstitutionModel("inference.independent.model", dataType,
                            rootFrequencies, stateCount);

            this.rewardsAwareBranchModel = new RewardsAwareBranchModel(
                    tree,
                    inferenceModel,
                    rewardRates,
                    indicator,
                    branchRates,
                    atomIndices,
                    false
            );

            final SimpleAlignment alignment = simulation.getGeneratedTipsAlignment();
            final SitePatterns patterns = new SitePatterns(alignment, null, 0, -1, 1, true);
            final GammaSiteRateModel siteRateModel =
                    new GammaSiteRateModel("independent.operator.validation.siteModel");
            final DiscreteDataLikelihoodDelegate delegate = new DiscreteDataLikelihoodDelegate(
                    tree,
                    patterns,
                    rewardsAwareBranchModel,
                    siteRateModel,
                    false,
                    false,
                    PartialsRescalingScheme.NONE,
                    false,
                    new PreOrderSettings(true, false, false, false),
                    new RewardsAwarePartialsRepresentation(rewardsAwareBranchModel),
                    DiscreteDataLikelihoodDelegate.PartialTransform.IDENTITY,
                    DiscreteDataLikelihoodDelegate.PartialTransform.IDENTITY,
                    null,
                    null
            );
            this.independentLikelihood =
                    new TreeDataLikelihood(delegate, tree, new DefaultBranchRateModel());

            this.indicatorOperator = new RewardsMixtureIndicatorAndAtomIndicesOperator(
                    indicator,
                    atomIndices,
                    rewardsAwareBranchModel,
                    independentLikelihood,
                    1.0 / Math.max(1, trueCts.length),
                    false,
                    false,
                    2,
                    0.5,
                    1.0
            );
            this.rewardRateHmc = createRewardRateHmcKernel(independentLikelihood,
                    rewardsAwareBranchModel, rewardRates);
            this.ctsRandomWalk = new CtsRewardRandomWalk(ctsRewards, independentLikelihood,
                    Math.max(1, CTS_MOVES_PER_SWEEP == 0 ? trueCts.length : CTS_MOVES_PER_SWEEP),
                    CTS_RANDOM_WALK_SCALE);
            this.rewardOrderOperator = new RewardOrderOperator(rewardRates.getValues(),
                    rewardRates.getStateIndices(), independentLikelihood);
        }

        void execute() throws IOException {
            writeRunMetadata();
            writeTruth();
            final File logFile = new File(outputDir, "independent_operator_validation.log");
            final double initialLogTarget = logTarget(independentLikelihood);

            try (PrintWriter log = new PrintWriter(logFile)) {
                writeLogHeader(log);
                for (int sweep = 0; sweep <= sweeps; sweep++) {
                    if (sweep > 0) {
                        step();
                    }
                    if (sweep >= burnIn && ((sweep - burnIn) % logEvery == 0 || sweep == sweeps)) {
                        writeLogRow(log, sweep);
                        sampleCount++;
                    }
                    if (sweep > 0 && sweep % Math.max(logEvery, 1000) == 0) {
                        log.flush();
                    }
                }
            }

            writeSummary(initialLogTarget, logTarget(independentLikelihood));
        }

        private void step() {
            if (scenario.updateIndicators) {
                indicatorOperator.doOperation();
            }
            if (scenario.updateCtsRewards) {
                ctsRandomWalk.doOperation();
            }
            if (scenario.updateRewardOrder) {
                rewardOrderOperator.doOperation();
            }
            if (scenario.updateRewardValues) {
                for (int h = 0; h < REWARD_RATE_HMC_STEPS; h++) {
                    rewardRateHmc.doOperation();
                }
            }
        }

        private void writeRunMetadata() throws IOException {
            try (PrintWriter out = new PrintWriter(new File(outputDir, "run_metadata.tsv"))) {
                out.println("field\tvalue");
                out.println("seed\t" + seed);
                out.println("tipCount\t" + tipCount);
                out.println("stateCount\t" + stateCount);
                out.println("scenarioId\t" + scenario.id);
                out.println("scenario\t" + scenario.label);
                out.println("sweeps\t" + sweeps);
                out.println("burnIn\t" + burnIn);
                out.println("logEvery\t" + logEvery);
                out.println("branchCount\t" + trueCts.length);
                out.println("updateRewardValues\t" + scenario.updateRewardValues);
                out.println("updateRewardOrder\t" + scenario.updateRewardOrder);
                out.println("updateCtsRewards\t" + scenario.updateCtsRewards);
                out.println("updateIndicators\t" + scenario.updateIndicators);
            }
        }

        private void writeTruth() throws IOException {
            try (PrintWriter out = new PrintWriter(new File(outputDir, "truth.tsv"))) {
                out.println("type\tindex\tvalue");
                for (int i = 0; i < truthRewardValues.length; i++) {
                    out.println("rewardValue\t" + i + "\t" + truthRewardValues[i]);
                }
                for (int i = 0; i < truthMapping.length; i++) {
                    out.println("rewardMapping\t" + i + "\t" + truthMapping[i]);
                }
                for (int i = 0; i < truthRewardByState.length; i++) {
                    out.println("rewardByState\t" + i + "\t" + truthRewardByState[i]);
                }
                for (int i = 0; i < trueCts.length; i++) {
                    out.println("ctsReward\t" + i + "\t" + trueCts[i]);
                    out.println("indicator\t" + i + "\t" + trueIndicators[i]);
                    out.println("atom\t" + i + "\t" + trueAtoms[i]);
                }
            }
            writeNewick(new File(outputDir, "tree.nwk"), tree);
            writeAlignment(new File(outputDir, "independent_tips.fasta"), simulation.getGeneratedTipsAlignment());
        }

        private void writeLogHeader(final PrintWriter out) {
            out.print("sample");
            out.print('\t');
            out.print("logLikelihood");
            for (int i = 0; i < stateCount; i++) {
                out.print('\t');
                out.print("rewardValue.");
                out.print(i);
            }
            for (int i = 0; i < stateCount; i++) {
                out.print('\t');
                out.print("rewardMapping.");
                out.print(i);
            }
            out.print('\t');
            out.print("ctsRMSE");
            out.print('\t');
            out.print("indicatorAccuracy");
            out.print('\t');
            out.print("indicatorBrier");
            for (int i = 0; i < indicator.getDimension(); i++) {
                out.print('\t');
                out.print("indicator.");
                out.print(i);
            }
            for (int i = 0; i < ctsRewards.getDimension(); i++) {
                out.print('\t');
                out.print("ctsReward.");
                out.print(i);
            }
            out.println();
        }

        private void writeLogRow(final PrintWriter out, final int sweep) {
            out.print(sweep);
            out.print('\t');
            out.print(logTarget(independentLikelihood));
            for (int i = 0; i < rewardRates.getValues().getDimension(); i++) {
                out.print('\t');
                out.print(rewardRates.getValues().getParameterValue(i));
            }
            for (int i = 0; i < rewardRates.getStateIndices().getDimension(); i++) {
                out.print('\t');
                out.print(rewardRates.getStateIndices().getParameterValue(i));
            }
            out.print('\t');
            out.print(ctsRmse());
            out.print('\t');
            out.print(indicatorAccuracy());
            out.print('\t');
            out.print(indicatorBrier());
            for (int i = 0; i < indicator.getDimension(); i++) {
                out.print('\t');
                out.print(indicator.getParameterValue(i));
            }
            for (int i = 0; i < ctsRewards.getDimension(); i++) {
                out.print('\t');
                out.print(ctsRewards.getParameterValue(i));
            }
            out.println();
        }

        private void writeSummary(final double initialLogTarget,
                                  final double finalLogTarget) throws IOException {
            try (PrintWriter out = new PrintWriter(new File(outputDir, "summary.tsv"))) {
                out.println("metric\tvalue");
                out.println("initialLogLikelihood\t" + initialLogTarget);
                out.println("finalLogLikelihood\t" + finalLogTarget);
                out.println("sampleCount\t" + sampleCount);
                out.println("finalCtsRMSE\t" + ctsRmse());
                out.println("finalIndicatorAccuracy\t" + indicatorAccuracy());
                out.println("finalIndicatorBrier\t" + indicatorBrier());
                out.println("rewardRateHmcAttemptCount\t" + rewardRateHmc.getAttemptCount());
                out.println("rewardRateHmcAcceptCount\t" + rewardRateHmc.getAcceptCount());
                out.println("rewardRateHmcAcceptanceRate\t" +
                        ratio(rewardRateHmc.getAcceptCount(), rewardRateHmc.getAttemptCount()));
                out.println("ctsAttemptCount\t" + ctsRandomWalk.getAttemptCount());
                out.println("ctsAcceptCount\t" + ctsRandomWalk.getAcceptCount());
                out.println("ctsAcceptanceRate\t" +
                        ratio(ctsRandomWalk.getAcceptCount(), ctsRandomWalk.getAttemptCount()));
                out.println("rewardOrderAttemptCount\t" + rewardOrderOperator.getAttemptCount());
                out.println("rewardOrderAcceptCount\t" + rewardOrderOperator.getAcceptCount());
                out.println("rewardOrderAcceptanceRate\t" +
                        ratio(rewardOrderOperator.getAcceptCount(), rewardOrderOperator.getAttemptCount()));
            }
        }

        private double ctsRmse() {
            double sum = 0.0;
            for (int i = 0; i < trueCts.length; i++) {
                final double error = ctsRewards.getParameterValue(i) - trueCts[i];
                sum += error * error;
            }
            return Math.sqrt(sum / trueCts.length);
        }

        private double indicatorAccuracy() {
            int correct = 0;
            for (int i = 0; i < trueIndicators.length; i++) {
                final int value = (int) Math.round(indicator.getParameterValue(i));
                if (value == trueIndicators[i]) {
                    correct++;
                }
            }
            return correct / (double) trueIndicators.length;
        }

        private double indicatorBrier() {
            double sum = 0.0;
            for (int i = 0; i < trueIndicators.length; i++) {
                final double error = indicator.getParameterValue(i) - trueIndicators[i];
                sum += error * error;
            }
            return sum / trueIndicators.length;
        }
    }

    private static TreeModel createBalancedTree(final int tipCount) {
        final SimpleNode root = buildBalancedSubtree(0, tipCount, heightForTipCount(tipCount));
        final SimpleTree simpleTree = new SimpleTree(root);
        return new DefaultTreeModel("balancedTree." + tipCount, simpleTree);
    }

    private static SimpleNode buildBalancedSubtree(final int fromInclusive,
                                                   final int toExclusive,
                                                   final double height) {
        final int count = toExclusive - fromInclusive;
        final SimpleNode node = new SimpleNode();
        node.setHeight(height);
        if (count == 1) {
            node.setTaxon(new Taxon("taxon" + String.format("%03d", fromInclusive + 1)));
            node.setHeight(0.0);
            return node;
        }
        final int leftCount = count / 2;
        final int split = fromInclusive + leftCount;
        node.addChild(buildBalancedSubtree(fromInclusive, split, heightForTipCount(leftCount)));
        node.addChild(buildBalancedSubtree(split, toExclusive, heightForTipCount(count - leftCount)));
        return node;
    }

    private static double heightForTipCount(final int tipCount) {
        int levels = 0;
        int n = 1;
        while (n < tipCount) {
            n <<= 1;
            levels++;
        }
        return levels * BRANCH_HEIGHT_UNIT;
    }

    private static DataType createDataType(final int stateCount) {
        final String[] codes = new String[stateCount];
        for (int i = 0; i < stateCount; i++) {
            codes[i] = stateCode(i);
        }
        final GeneralDataType dataType = new GeneralDataType(codes);
        dataType.setId("rewardStateDataType." + stateCount);
        return dataType;
    }

    private static String stateCode(final int state) {
        if (state < 10) {
            return Integer.toString(state);
        }
        return "s" + state;
    }

    private static double[] uniform(final int dimension) {
        final double[] values = new double[dimension];
        Arrays.fill(values, 1.0 / dimension);
        return values;
    }

    private static double[] createTruthRewardValues(final int stateCount) {
        final double[] values = new double[stateCount];
        values[0] = 0.0;
        values[1] = 1.0;
        for (int i = 2; i < stateCount; i++) {
            values[i] = (i - 1.0) / (stateCount - 1.0);
        }
        return values;
    }

    private static int[] createTruthMapping(final int stateCount) {
        final int[] mapping = new int[stateCount];
        for (int i = 0; i < stateCount; i++) {
            mapping[i] = i;
        }
        if (stateCount > 1) {
            swap(mapping, 0, 1);
        }
        if (stateCount > 3) {
            swap(mapping, 2, stateCount - 1);
        }
        return mapping;
    }

    private static void swap(final int[] values, final int i, final int j) {
        final int tmp = values[i];
        values[i] = values[j];
        values[j] = tmp;
    }

    private static double[] valuesByState(final double[] rewardValues, final int[] mapping) {
        final double[] byState = new double[mapping.length];
        for (int i = 0; i < mapping.length; i++) {
            byState[i] = rewardValues[mapping[i]];
        }
        return byState;
    }

    private static double[] initializeCtsRewards(final double[] trueCts, final boolean update) {
        if (!update) {
            return Arrays.copyOf(trueCts, trueCts.length);
        }
        final double[] values = new double[trueCts.length];
        Arrays.fill(values, 0.5);
        return values;
    }

    private static double[] initializeIndicators(final int[] trueIndicators, final boolean update) {
        final double[] values = new double[trueIndicators.length];
        if (!update) {
            for (int i = 0; i < values.length; i++) {
                values[i] = trueIndicators[i];
            }
            return values;
        }
        Arrays.fill(values, 0.0);
        return values;
    }

    private static double[] initializeAtoms(final int[] trueAtoms,
                                            final boolean update,
                                            final int stateCount) {
        final double[] values = new double[trueAtoms.length];
        if (!update) {
            for (int i = 0; i < values.length; i++) {
                values[i] = trueAtoms[i];
            }
            return values;
        }
        for (int i = 0; i < values.length; i++) {
            values[i] = i % stateCount;
        }
        return values;
    }

    private static RewardRates createInferenceRewardRates(final String id,
                                                          final double[] truthValues,
                                                          final int[] truthMapping,
                                                          final boolean updateValues,
                                                          final boolean updateOrder) {
        final Parameter fixedValues = new Parameter.Default(id + ".fixed", new double[]{0.0, 1.0});
        final double[] internalValues = new double[Math.max(0, truthValues.length - 2)];
        for (int i = 0; i < internalValues.length; i++) {
            final double truth = truthValues[i + 2];
            internalValues[i] = updateValues ? clamp01(0.35 + 0.25 * (i % 2) + 0.1 * truth) : truth;
        }
        final Parameter varyingValues = new Parameter.Default(id + ".varying", internalValues);
        varyingValues.addBounds(new Parameter.DefaultBounds(1.0, 0.0, varyingValues.getDimension()));

        final CompoundParameter values = new CompoundParameter(id + ".values");
        values.addParameter(fixedValues);
        values.addParameter(varyingValues);

        final double[] mapping = new double[truthMapping.length];
        if (updateOrder) {
            for (int i = 0; i < mapping.length; i++) {
                mapping[i] = i;
            }
        } else {
            for (int i = 0; i < mapping.length; i++) {
                mapping[i] = truthMapping[i];
            }
        }

        return new RewardRates(
                values,
                fixedValues,
                varyingValues,
                new Parameter.Default(id + ".mapping", mapping)
        );
    }

    private static LogRateSubstitutionModel createIndependentRewardSubstitutionModel(
            final String id,
            final DataType dataType,
            final double[] rootFrequencies,
            final int stateCount) {

        final double[] logRates = new double[stateCount * (stateCount - 1)];
        int index = 0;
        for (int from = 0; from < stateCount; from++) {
            for (int to = 0; to < stateCount; to++) {
                if (from == to) {
                    continue;
                }
                final double rate = 0.20 + (1 + ((from + 2 * to) % stateCount)) / (double) stateCount;
                logRates[index++] = Math.log(rate);
            }
        }

        final FrequencyModel frequencyModel = new FrequencyModel(dataType, rootFrequencies);
        final LogRateSubstitutionModel model = new LogRateSubstitutionModel(
                id,
                dataType,
                frequencyModel,
                new LogAdditiveCtmcRateProvider.DataAugmented.Basic(
                        id + ".rateProvider",
                        new Parameter.Default(id + ".logRates", logRates)
                )
        );
        model.setNormalization(true);
        model.setScaleRatesByFrequencies(false);
        return model;
    }

    private static RewardRateHmcKernel createRewardRateHmcKernel(
            final TreeDataLikelihood independentLikelihood,
            final RewardsAwareBranchModel rewardsAwareBranchModel,
            final RewardRates rewardRates) {

        final List<GradientWrtParameterProvider> gradients =
                new ArrayList<GradientWrtParameterProvider>();
        gradients.add(new RewardsAwareBranchModelRewardRateGradient(
                independentLikelihood,
                rewardsAwareBranchModel,
                rewardRates.getValues(),
                rewardRates.getVaryingValues(),
                null));

        return new RewardRateHmcKernel(
                rewardRates.getVaryingValues(),
                gradients,
                new LogTargetEvaluator() {
                    @Override
                    public double evaluate() {
                        return logTarget(independentLikelihood);
                    }
                },
                0.0,
                1.0,
                REWARD_RATE_HMC_STEP_SIZE,
                REWARD_RATE_HMC_LEAPFROG_STEPS);
    }

    private static double logTarget(final TreeDataLikelihood likelihood) {
        likelihood.makeDirty();
        final double logLikelihood = likelihood.getLogLikelihood();
        if (Double.isNaN(logLikelihood) || logLikelihood == Double.POSITIVE_INFINITY) {
            throw new IllegalStateException("Invalid log likelihood: " + logLikelihood);
        }
        return logLikelihood;
    }

    private static double clamp01(final double value) {
        return Math.max(0.0, Math.min(1.0, value));
    }

    private static double ratio(final int numerator, final int denominator) {
        return denominator == 0 ? Double.NaN : numerator / (double) denominator;
    }

    private static void writeNewick(final File file, final TreeModel tree) throws IOException {
        try (PrintWriter out = new PrintWriter(file)) {
            out.println(toNewick(tree, tree.getRoot()) + ";");
        }
    }

    private static String toNewick(final Tree tree, final NodeRef node) {
        final StringBuilder sb = new StringBuilder();
        if (tree.isExternal(node)) {
            sb.append(tree.getNodeTaxon(node).getId());
        } else {
            sb.append('(');
            for (int i = 0; i < tree.getChildCount(node); i++) {
                if (i > 0) {
                    sb.append(',');
                }
                sb.append(toNewick(tree, tree.getChild(node, i)));
            }
            sb.append(')');
        }
        if (!tree.isRoot(node)) {
            sb.append(':').append(tree.getBranchLength(node));
        }
        return sb.toString();
    }

    private static void writeAlignment(final File file, final SimpleAlignment alignment) throws IOException {
        try (PrintWriter out = new PrintWriter(file)) {
            for (int i = 0; i < alignment.getSequenceCount(); i++) {
                final Sequence sequence = alignment.getSequence(i);
                out.println('>' + sequence.getTaxon().getId());
                out.println(sequence.getSequenceString());
            }
        }
    }

    private interface LogTargetEvaluator {
        double evaluate();
    }

    private static final class RewardOrderOperator {
        private final Parameter rewardValues;
        private final Parameter mapping;
        private final TreeDataLikelihood likelihood;
        private int attemptCount = 0;
        private int acceptCount = 0;

        private RewardOrderOperator(final Parameter rewardValues,
                                    final Parameter mapping,
                                    final TreeDataLikelihood likelihood) {
            this.rewardValues = rewardValues;
            this.mapping = mapping;
            this.likelihood = likelihood;
        }

        private void doOperation() {
            attemptCount++;
            final double oldLogTarget = logTarget(likelihood);
            final double[] oldMapping = mapping.getParameterValues();
            final int k = rewardValues.getDimension();
            final int extreme = MathUtils.nextBoolean() ? 1 : 0;
            final int other = k == 2 ? 1 - extreme : 2 + MathUtils.nextInt(k - 2);
            remapRewardLabelsBySwap(extreme, other);
            mapping.fireParameterChangedEvent();
            final double newLogTarget = logTarget(likelihood);
            if (Double.isFinite(newLogTarget) &&
                    Math.log(MathUtils.nextDouble()) < newLogTarget - oldLogTarget) {
                acceptCount++;
            } else {
                restore(mapping, oldMapping);
            }
        }

        private void remapRewardLabelsBySwap(final int i, final int j) {
            for (int s = 0; s < mapping.getDimension(); s++) {
                final int label = (int) Math.round(mapping.getParameterValue(s));
                if (label == i) {
                    mapping.setParameterValueQuietly(s, j);
                } else if (label == j) {
                    mapping.setParameterValueQuietly(s, i);
                }
            }
        }

        private int getAttemptCount() {
            return attemptCount;
        }

        private int getAcceptCount() {
            return acceptCount;
        }
    }

    private static final class CtsRewardRandomWalk {
        private final Parameter ctsRewards;
        private final TreeDataLikelihood likelihood;
        private final int movesPerSweep;
        private final double scale;
        private int attemptCount = 0;
        private int acceptCount = 0;

        private CtsRewardRandomWalk(final Parameter ctsRewards,
                                    final TreeDataLikelihood likelihood,
                                    final int movesPerSweep,
                                    final double scale) {
            this.ctsRewards = ctsRewards;
            this.likelihood = likelihood;
            this.movesPerSweep = movesPerSweep;
            this.scale = scale;
        }

        private void doOperation() {
            for (int i = 0; i < movesPerSweep; i++) {
                attemptCount++;
                final int index = MathUtils.nextInt(ctsRewards.getDimension());
                final double oldValue = ctsRewards.getParameterValue(index);
                final double oldLogTarget = logTarget(likelihood);
                final double newValue = reflect01(oldValue + scale * MathUtils.nextGaussian());
                ctsRewards.setParameterValue(index, newValue);
                final double newLogTarget = logTarget(likelihood);
                if (Double.isFinite(newLogTarget) &&
                        Math.log(MathUtils.nextDouble()) < newLogTarget - oldLogTarget) {
                    acceptCount++;
                } else {
                    ctsRewards.setParameterValue(index, oldValue);
                }
            }
        }

        private static double reflect01(final double raw) {
            double value = raw;
            while (value < 0.0 || value > 1.0) {
                if (value < 0.0) {
                    value = -value;
                }
                if (value > 1.0) {
                    value = 2.0 - value;
                }
            }
            return value;
        }

        private int getAttemptCount() {
            return attemptCount;
        }

        private int getAcceptCount() {
            return acceptCount;
        }
    }

    private static final class RewardRateHmcKernel {

        private static final double TRANSFORM_EPSILON = 1.0e-12;

        private final Parameter rewardRateVaryingValues;
        private final List<GradientWrtParameterProvider> gradients;
        private final LogTargetEvaluator target;
        private final double lower;
        private final double upper;
        private final double width;
        private final double stepSize;
        private final int leapfrogSteps;

        private int attemptCount = 0;
        private int acceptCount = 0;

        private RewardRateHmcKernel(final Parameter rewardRateVaryingValues,
                                    final List<GradientWrtParameterProvider> gradients,
                                    final LogTargetEvaluator target,
                                    final double lower,
                                    final double upper,
                                    final double stepSize,
                                    final int leapfrogSteps) {
            this.rewardRateVaryingValues = rewardRateVaryingValues;
            this.gradients = gradients;
            this.target = target;
            this.lower = lower;
            this.upper = upper;
            this.width = upper - lower;
            this.stepSize = stepSize;
            this.leapfrogSteps = leapfrogSteps;
        }

        private void doOperation() {
            if (rewardRateVaryingValues.getDimension() == 0) {
                return;
            }
            attemptCount++;
            final double[] oldValues = rewardRateVaryingValues.getParameterValues();
            final double[] theta = valuesToTheta(oldValues);
            final double oldLogDensity = transformedLogDensity(theta);
            if (!Double.isFinite(oldLogDensity)) {
                restoreValues(oldValues);
                return;
            }

            final double[] momentum = new double[theta.length];
            for (int i = 0; i < momentum.length; i++) {
                momentum[i] = MathUtils.nextGaussian();
            }
            final double oldKineticEnergy = kineticEnergy(momentum);

            double[] gradient = transformedGradient(theta);
            if (!allFinite(gradient)) {
                restoreValues(oldValues);
                return;
            }
            for (int i = 0; i < momentum.length; i++) {
                momentum[i] += 0.5 * stepSize * gradient[i];
            }

            for (int step = 0; step < leapfrogSteps; step++) {
                for (int i = 0; i < theta.length; i++) {
                    theta[i] += stepSize * momentum[i];
                }
                gradient = transformedGradient(theta);
                if (!allFinite(gradient)) {
                    restoreValues(oldValues);
                    return;
                }
                final double scale = step == leapfrogSteps - 1 ? 0.5 : 1.0;
                for (int i = 0; i < momentum.length; i++) {
                    momentum[i] += scale * stepSize * gradient[i];
                }
            }

            final double newLogDensity = transformedLogDensity(theta);
            final double newKineticEnergy = kineticEnergy(momentum);
            final double logAcceptanceProbability =
                    newLogDensity - oldLogDensity + oldKineticEnergy - newKineticEnergy;
            if (Double.isFinite(logAcceptanceProbability) &&
                    Math.log(MathUtils.nextDouble()) < logAcceptanceProbability) {
                acceptCount++;
            } else {
                restoreValues(oldValues);
            }
        }

        private double transformedLogDensity(final double[] theta) {
            setValuesFromTheta(theta);
            final double logLikelihood = target.evaluate();
            if (!Double.isFinite(logLikelihood)) {
                return Double.NEGATIVE_INFINITY;
            }
            double logJacobian = 0.0;
            for (double value : theta) {
                final double s = logistic(value);
                logJacobian += Math.log(width) + Math.log(s) + Math.log1p(-s);
            }
            return logLikelihood + logJacobian;
        }

        private double[] transformedGradient(final double[] theta) {
            setValuesFromTheta(theta);
            final double logLikelihood = target.evaluate();
            if (!Double.isFinite(logLikelihood)) {
                final double[] invalid = new double[theta.length];
                Arrays.fill(invalid, Double.NaN);
                return invalid;
            }
            final double[] valueGradient = new double[theta.length];
            for (GradientWrtParameterProvider gradientProvider : gradients) {
                final double[] contribution = gradientProvider.getGradientLogDensity();
                if (contribution.length != valueGradient.length) {
                    throw new IllegalStateException("Reward-rate gradient dimension mismatch.");
                }
                for (int i = 0; i < valueGradient.length; i++) {
                    valueGradient[i] += contribution[i];
                }
            }
            final double[] thetaGradient = new double[theta.length];
            for (int i = 0; i < theta.length; i++) {
                final double s = logistic(theta[i]);
                final double jacobian = width * s * (1.0 - s);
                thetaGradient[i] = valueGradient[i] * jacobian + 1.0 - 2.0 * s;
            }
            return thetaGradient;
        }

        private void setValuesFromTheta(final double[] theta) {
            for (int i = 0; i < theta.length; i++) {
                final double s = logistic(theta[i]);
                rewardRateVaryingValues.setParameterValueQuietly(i, lower + width * s);
            }
            rewardRateVaryingValues.fireParameterChangedEvent();
        }

        private double[] valuesToTheta(final double[] values) {
            final double[] theta = new double[values.length];
            for (int i = 0; i < values.length; i++) {
                double s = (values[i] - lower) / width;
                s = Math.max(TRANSFORM_EPSILON, Math.min(1.0 - TRANSFORM_EPSILON, s));
                theta[i] = Math.log(s) - Math.log1p(-s);
            }
            return theta;
        }

        private void restoreValues(final double[] values) {
            restore(rewardRateVaryingValues, values);
        }

        private static double logistic(final double x) {
            if (x >= 0.0) {
                final double z = Math.exp(-x);
                return 1.0 / (1.0 + z);
            }
            final double z = Math.exp(x);
            return z / (1.0 + z);
        }

        private static double kineticEnergy(final double[] momentum) {
            double sum = 0.0;
            for (double value : momentum) {
                sum += value * value;
            }
            return 0.5 * sum;
        }

        private static boolean allFinite(final double[] values) {
            for (double value : values) {
                if (!Double.isFinite(value)) {
                    return false;
                }
            }
            return true;
        }

        private int getAttemptCount() {
            return attemptCount;
        }

        private int getAcceptCount() {
            return acceptCount;
        }
    }

    private static void restore(final Parameter parameter, final double[] values) {
        for (int i = 0; i < values.length; i++) {
            parameter.setParameterValueQuietly(i, values[i]);
        }
        parameter.fireParameterChangedEvent();
    }
}
