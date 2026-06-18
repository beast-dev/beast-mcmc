package dr.evomodel.branchmodel;

import dr.evolution.alignment.SimpleAlignment;
import dr.evolution.datatype.DataType;
import dr.evolution.datatype.GeneralDataType;
import dr.evolution.sequence.DelimitedSequence;
import dr.evolution.sequence.Sequence;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxon;
import dr.inference.markovjumps.StateChange;
import dr.inference.markovjumps.StateHistory;
import dr.inference.model.Parameter;
import dr.math.MathUtils;
import dr.xml.Reportable;

import java.util.Arrays;

/*
 * @author Filippo Monti
 */
public final class RewardsAwareCtmcSimulation implements Reportable {

    private static final double Q_TOLERANCE = 1.0e-8;

    private final String id;
    private final int stateCount;
    private final Tree tree;
    private final double[] qMatrix;
    private final double[] rewardRates;
    private final double[] rootFrequencies;
    private final int rootState;

    private final int branchCount;
    private final int[] nodeNumberByBranchIndex;
    private final int[] branchIndexByNodeNumber;
    private final int[] nodeStates;

    private final double[] branchTimes;
    private final double[] branchRewardTotals;
    private final double[] branchRewardProportions;
    private final double[][] branchDwellTimes;
    private final int[][] branchTransitionCounts;
    private final int[] branchJumpCounts;
    private final int[] branchStartStates;
    private final int[] branchEndStates;
    private final int[] indicators;
    private final int[] atomIndices;
    private final int[] tipStates;
    private final int[] tipNodeNumbers;
    private final String[] tipTaxonIds;
    private final SimpleAlignment generatedTips;
    private final Summary summary;

    private final Parameter branchRewardTotalsParameter;
    private final Parameter ctsRewardsParameter;
    private final Parameter branchDwellTimesParameter;
    private final Parameter indicatorParameter;
    private final Parameter atomIndicesParameter;
    private final Parameter rewardRatesParameter;

    public RewardsAwareCtmcSimulation(final int stateCount,
                                      final double[] rewardRates,
                                      final Tree tree,
                                      final double[] qMatrix) {
        this(null, stateCount, rewardRates, tree, qMatrix, null, -1);
    }

    public RewardsAwareCtmcSimulation(final String id,
                                      final int stateCount,
                                      final double[] rewardRates,
                                      final Tree tree,
                                      final double[] qMatrix,
                                      final double[] rootFrequencies,
                                      final int rootState) {
        if (tree == null) {
            throw new IllegalArgumentException("tree must be non-null");
        }
        if (stateCount <= 0) {
            throw new IllegalArgumentException("stateCount must be positive: " + stateCount);
        }
        validateRewardRates(rewardRates, stateCount);
        validateQMatrix(qMatrix, stateCount);
        if (rootFrequencies != null && rootState >= 0) {
            throw new IllegalArgumentException("Specify either rootFrequencies or rootState, not both.");
        }

        this.id = id == null ? "rewardsAwareCtmcSimulation" : id;
        this.stateCount = stateCount;
        this.tree = tree;
        this.qMatrix = Arrays.copyOf(qMatrix, qMatrix.length);
        this.rewardRates = Arrays.copyOf(rewardRates, rewardRates.length);
        this.rootState = rootState >= 0 ? validateRootState(rootState, stateCount) : -1;
        this.rootFrequencies = rootFrequencies == null
                ? (this.rootState >= 0 ? pointMass(this.rootState, stateCount) : stationaryFrequencies(qMatrix, stateCount))
                : normalizeFrequencies(rootFrequencies, stateCount);

        this.branchCount = tree.getNodeCount() - 1;
        this.nodeNumberByBranchIndex = new int[branchCount];
        Arrays.fill(nodeNumberByBranchIndex, -1);
        this.branchIndexByNodeNumber = new int[tree.getNodeCount()];
        Arrays.fill(branchIndexByNodeNumber, -1);
        initializeBranchMappings();

        this.nodeStates = new int[tree.getNodeCount()];
        Arrays.fill(nodeStates, -1);

        this.branchTimes = new double[branchCount];
        this.branchRewardTotals = new double[branchCount];
        this.branchRewardProportions = new double[branchCount];
        this.branchDwellTimes = new double[branchCount][stateCount];
        this.branchTransitionCounts = new int[branchCount][stateCount * stateCount];
        this.branchJumpCounts = new int[branchCount];
        this.branchStartStates = new int[branchCount];
        this.branchEndStates = new int[branchCount];
        this.indicators = new int[branchCount];
        this.atomIndices = new int[branchCount];

        simulate();

        this.tipStates = new int[tree.getExternalNodeCount()];
        this.tipNodeNumbers = new int[tree.getExternalNodeCount()];
        this.tipTaxonIds = new String[tree.getExternalNodeCount()];
        collectTipStates();
        this.generatedTips = buildGeneratedTips();
        this.summary = summarize();

        this.branchRewardTotalsParameter =
                new Parameter.Default(this.id + ".branchRewardTotals", branchRewardTotals);
        this.ctsRewardsParameter =
                new Parameter.Default(this.id + ".branchRewardProportions", branchRewardProportions);
        this.branchDwellTimesParameter =
                new Parameter.Default(this.id + ".branchDwellTimes", flatten(branchDwellTimes));
        this.indicatorParameter =
                new Parameter.Default(this.id + ".indicator", toDouble(indicators));
        this.atomIndicesParameter =
                new Parameter.Default(this.id + ".atomIndices", toDouble(atomIndices));
        this.rewardRatesParameter =
                new Parameter.Default(this.id + ".rewardRates", this.rewardRates);
    }

    public static RewardsAwareCtmcSimulation withFixedRootState(final String id,
                                                                final int stateCount,
                                                                final double[] rewardRates,
                                                                final Tree tree,
                                                                final double[] qMatrix,
                                                                final int rootState) {
        return new RewardsAwareCtmcSimulation(id, stateCount, rewardRates, tree, qMatrix, null, rootState);
    }

    public static RewardsAwareCtmcSimulation withRootFrequencies(final String id,
                                                                 final int stateCount,
                                                                 final double[] rewardRates,
                                                                 final Tree tree,
                                                                 final double[] qMatrix,
                                                                 final double[] rootFrequencies) {
        return new RewardsAwareCtmcSimulation(id, stateCount, rewardRates, tree, qMatrix, rootFrequencies, -1);
    }

    public int getStateCount() {
        return stateCount;
    }

    public Tree getTree() {
        return tree;
    }

    public double[] getQMatrix() {
        return Arrays.copyOf(qMatrix, qMatrix.length);
    }

    public double[] getRewardRates() {
        return Arrays.copyOf(rewardRates, rewardRates.length);
    }

    public double[] getRootFrequencies() {
        return Arrays.copyOf(rootFrequencies, rootFrequencies.length);
    }

    public int getRootState() {
        return rootState < 0 ? nodeStates[tree.getRoot().getNumber()] : rootState;
    }

    public double[] getBranchRewardTotals() {
        return Arrays.copyOf(branchRewardTotals, branchRewardTotals.length);
    }

    public double[] getBranchRewardProportions() {
        return Arrays.copyOf(branchRewardProportions, branchRewardProportions.length);
    }

    public int[] getIndicators() {
        return Arrays.copyOf(indicators, indicators.length);
    }

    public int[] getAtomIndices() {
        return Arrays.copyOf(atomIndices, atomIndices.length);
    }

    public int[] getGeneratedTipStates() {
        return Arrays.copyOf(tipStates, tipStates.length);
    }

    public int[] getNodeStates() {
        return Arrays.copyOf(nodeStates, nodeStates.length);
    }

    public int[] getBranchJumpCounts() {
        return Arrays.copyOf(branchJumpCounts, branchJumpCounts.length);
    }

    public double[][] getBranchDwellTimes() {
        return copy(branchDwellTimes);
    }

    public int[][] getBranchTransitionCounts() {
        final int[][] out = new int[branchTransitionCounts.length][];
        for (int i = 0; i < out.length; i++) {
            out[i] = Arrays.copyOf(branchTransitionCounts[i], branchTransitionCounts[i].length);
        }
        return out;
    }

    public int[] getBranchStartStates() {
        return Arrays.copyOf(branchStartStates, branchStartStates.length);
    }

    public int[] getBranchEndStates() {
        return Arrays.copyOf(branchEndStates, branchEndStates.length);
    }

    public int getNodeNumberForBranchIndex(final int branchIndex) {
        if (branchIndex < 0 || branchIndex >= branchCount) {
            throw new IllegalArgumentException("branchIndex out of range: " + branchIndex);
        }
        return nodeNumberByBranchIndex[branchIndex];
    }

    public int getBranchIndexForNodeNumber(final int nodeNumber) {
        if (nodeNumber < 0 || nodeNumber >= branchIndexByNodeNumber.length) {
            throw new IllegalArgumentException("nodeNumber out of range: " + nodeNumber);
        }
        return branchIndexByNodeNumber[nodeNumber];
    }

    public SimpleAlignment getGeneratedTipsAlignment() {
        return generatedTips;
    }

    public Summary getSummary() {
        return summary;
    }

    public Parameter getBranchRewardTotalsParameter() {
        return branchRewardTotalsParameter;
    }

    public Parameter getCtsRewardsParameter() {
        return ctsRewardsParameter;
    }

    public Parameter getBranchDwellTimesParameter() {
        return branchDwellTimesParameter;
    }

    public Parameter getIndicatorParameter() {
        return indicatorParameter;
    }

    public Parameter getAtomIndicesParameter() {
        return atomIndicesParameter;
    }

    public Parameter getRewardRatesParameter() {
        return rewardRatesParameter;
    }

    public void copyBranchRewardTotalsInto(final Parameter parameter) {
        copyIntoParameter(parameter, branchRewardTotals);
    }

    public void copyBranchRewardProportionsInto(final Parameter parameter) {
        copyIntoParameter(parameter, branchRewardProportions);
    }

    public void copyBranchDwellTimesInto(final Parameter parameter) {
        copyIntoParameter(parameter, flatten(branchDwellTimes));
    }

    public void copyIndicatorsInto(final Parameter parameter) {
        copyIntoParameter(parameter, toDouble(indicators));
    }

    public void copyAtomIndicesInto(final Parameter parameter) {
        copyIntoParameter(parameter, toDouble(atomIndices));
    }

    public static void copyIntoParameter(final Parameter parameter, final double[] values) {
        if (parameter == null) {
            return;
        }
        if (parameter.getDimension() != values.length) {
            throw new IllegalArgumentException(
                    "Parameter " + parameter.getParameterName() + " has dimension " +
                            parameter.getDimension() + ", expected " + values.length
            );
        }
        for (int i = 0; i < values.length; i++) {
            parameter.setParameterValueQuietly(i, values[i]);
        }
        parameter.fireParameterChangedEvent();
    }

    private void initializeBranchMappings() {
        for (int i = 0; i < tree.getNodeCount(); i++) {
            final NodeRef node = tree.getNode(i);
            if (tree.isRoot(node)) {
                continue;
            }
            final int nodeNumber = node.getNumber();
            final int branchIndex = getWithoutRootParameterIndex(nodeNumber);
            nodeNumberByBranchIndex[branchIndex] = nodeNumber;
            branchIndexByNodeNumber[nodeNumber] = branchIndex;
        }
        for (int i = 0; i < nodeNumberByBranchIndex.length; i++) {
            if (nodeNumberByBranchIndex[i] < 0) {
                throw new IllegalStateException("No node number mapped to branch index " + i);
            }
        }
    }

    private int getWithoutRootParameterIndex(final int nodeNumber) {
        final int rootNumber = tree.getRoot().getNumber();
        return nodeNumber > rootNumber ? nodeNumber - 1 : nodeNumber;
    }

    private void simulate() {
        final int simulatedRootState = rootState >= 0 ? rootState : MathUtils.randomChoicePDF(rootFrequencies);
        final NodeRef root = tree.getRoot();
        nodeStates[root.getNumber()] = simulatedRootState;
        traverse(root, simulatedRootState);
    }

    private void traverse(final NodeRef parent, final int parentState) {
        final int childCount = tree.getChildCount(parent);
        for (int i = 0; i < childCount; i++) {
            final NodeRef child = tree.getChild(parent, i);
            final int branchIndex = branchIndexByNodeNumber[child.getNumber()];
            final double branchTime = tree.getBranchLength(child);
            if (branchTime < -Q_TOLERANCE) {
                throw new IllegalArgumentException("Negative branch length for node " + child.getNumber() + ": " + branchTime);
            }
            final double nonNegativeBranchTime = Math.max(0.0, branchTime);
            final StateHistory history = simulateAlongBranch(parentState, nonNegativeBranchTime);
            processBranch(branchIndex, nonNegativeBranchTime, history);

            final int childState = history.getEndingState();
            nodeStates[child.getNumber()] = childState;
            traverse(child, childState);
        }
    }

    private StateHistory simulateAlongBranch(final int startingState, final double branchTime) {
        if (branchTime == 0.0 || -qMatrix[startingState * stateCount + startingState] <= 0.0) {
            final StateHistory history = new StateHistory(0.0, startingState, stateCount);
            history.addEndingState(new StateChange(branchTime, startingState));
            return history;
        }
        return StateHistory.simulateUnconditionalOnEndingState(
                0.0,
                startingState,
                branchTime,
                qMatrix,
                stateCount
        );
    }

    private void processBranch(final int branchIndex,
                               final double branchTime,
                               final StateHistory history) {
        branchTimes[branchIndex] = branchTime;
        branchStartStates[branchIndex] = history.getStartingState();
        branchEndStates[branchIndex] = history.getEndingState();
        branchJumpCounts[branchIndex] = history.getNumberOfJumps();

        final double[] waitingTimes = history.getWaitingTimes();
        final int[] jumpCounts = history.getJumpCounts();
        System.arraycopy(waitingTimes, 0, branchDwellTimes[branchIndex], 0, stateCount);
        System.arraycopy(jumpCounts, 0, branchTransitionCounts[branchIndex], 0, stateCount * stateCount);

        double totalReward = 0.0;
        for (int state = 0; state < stateCount; state++) {
            totalReward += waitingTimes[state] * rewardRates[state];
        }
        branchRewardTotals[branchIndex] = totalReward;
        branchRewardProportions[branchIndex] =
                branchTime > 0.0 ? totalReward / branchTime : rewardRates[history.getStartingState()];

        if (history.getNumberOfJumps() == 0) {
            indicators[branchIndex] = 1;
            atomIndices[branchIndex] = history.getStartingState();
        } else {
            indicators[branchIndex] = 0;
            atomIndices[branchIndex] = history.getStartingState();
        }
    }

    private void collectTipStates() {
        for (int i = 0; i < tree.getExternalNodeCount(); i++) {
            final NodeRef tip = tree.getExternalNode(i);
            tipNodeNumbers[i] = tip.getNumber();
            tipStates[i] = nodeStates[tip.getNumber()];
            final Taxon taxon = tree.getNodeTaxon(tip);
            tipTaxonIds[i] = taxon == null ? "node" + tip.getNumber() : taxon.getId();
        }
    }

    private SimpleAlignment buildGeneratedTips() {
        final String[] codes = new String[stateCount];
        for (int i = 0; i < stateCount; i++) {
            codes[i] = Integer.toString(i);
        }
        final GeneralDataType dataType = new GeneralDataType(codes);
        dataType.setId(id + ".dataType");

        final SimpleAlignment alignment = new SimpleAlignment();
        alignment.setDataType(dataType);
        for (int i = 0; i < tipStates.length; i++) {
            final Taxon taxon = new Taxon(tipTaxonIds[i]);
            final String sequenceString = dataType.getCode(tipStates[i]);
            final Sequence sequence = dataType.isDelimited()
                    ? new DelimitedSequence(taxon, sequenceString, dataType)
                    : new Sequence(taxon, sequenceString);
            sequence.setDataType(dataType);
            alignment.addSequence(sequence);
        }
        alignment.setId(id + ".generatedTips");
        return alignment;
    }

    private Summary summarize() {
        final double[] totalDwellTimes = new double[stateCount];
        final int[] transitionCounts = new int[stateCount * stateCount];
        final int[] tipStateCounts = new int[stateCount];

        int totalJumps = 0;
        int atomicBranchCount = 0;
        int continuousBranchCount = 0;
        double totalTreeTime = 0.0;
        double totalReward = 0.0;
        double minRewardProportion = Double.POSITIVE_INFINITY;
        double maxRewardProportion = Double.NEGATIVE_INFINITY;

        for (int branch = 0; branch < branchCount; branch++) {
            totalJumps += branchJumpCounts[branch];
            if (indicators[branch] == 1) {
                atomicBranchCount++;
            } else {
                continuousBranchCount++;
            }
            totalTreeTime += branchTimes[branch];
            totalReward += branchRewardTotals[branch];
            minRewardProportion = Math.min(minRewardProportion, branchRewardProportions[branch]);
            maxRewardProportion = Math.max(maxRewardProportion, branchRewardProportions[branch]);

            for (int state = 0; state < stateCount; state++) {
                totalDwellTimes[state] += branchDwellTimes[branch][state];
            }
            for (int i = 0; i < transitionCounts.length; i++) {
                transitionCounts[i] += branchTransitionCounts[branch][i];
            }
        }

        for (int tipState : tipStates) {
            tipStateCounts[tipState]++;
        }

        final double meanRewardProportion =
                totalTreeTime > 0.0 ? totalReward / totalTreeTime : Double.NaN;

        return new Summary(
                totalJumps,
                atomicBranchCount,
                continuousBranchCount,
                totalTreeTime,
                totalReward,
                meanRewardProportion,
                minRewardProportion,
                maxRewardProportion,
                totalDwellTimes,
                transitionCounts,
                tipStateCounts
        );
    }

    @Override
    public String getReport() {
        final StringBuilder sb = new StringBuilder();
        sb.append("RewardsAwareCtmcSimulation ").append(id).append('\n');
        sb.append(summary).append('\n');
        sb.append("branchRewardTotals=").append(Arrays.toString(branchRewardTotals)).append('\n');
        sb.append("branchRewardProportions=").append(Arrays.toString(branchRewardProportions)).append('\n');
        sb.append("branchDwellTimes=").append(Arrays.toString(flatten(branchDwellTimes))).append('\n');
        sb.append("indicators=").append(Arrays.toString(indicators)).append('\n');
        sb.append("atomIndices=").append(Arrays.toString(atomIndices)).append('\n');
        sb.append("tipStates=").append(Arrays.toString(tipStates)).append('\n');
        return sb.toString();
    }

    private static void validateRewardRates(final double[] rewardRates, final int stateCount) {
        if (rewardRates == null || rewardRates.length != stateCount) {
            throw new IllegalArgumentException("rewardRates must have length stateCount=" + stateCount);
        }
        for (int i = 0; i < rewardRates.length; i++) {
            if (!Double.isFinite(rewardRates[i])) {
                throw new IllegalArgumentException("rewardRates[" + i + "] is not finite: " + rewardRates[i]);
            }
        }
    }

    private static void validateQMatrix(final double[] qMatrix, final int stateCount) {
        if (qMatrix == null || qMatrix.length != stateCount * stateCount) {
            throw new IllegalArgumentException("qMatrix must have length stateCount^2=" + (stateCount * stateCount));
        }
        for (int row = 0; row < stateCount; row++) {
            double rowSum = 0.0;
            double offDiagonalSum = 0.0;
            for (int col = 0; col < stateCount; col++) {
                final double value = qMatrix[row * stateCount + col];
                if (!Double.isFinite(value)) {
                    throw new IllegalArgumentException("qMatrix[" + row + "," + col + "] is not finite: " + value);
                }
                rowSum += value;
                if (row == col) {
                    if (value > Q_TOLERANCE) {
                        throw new IllegalArgumentException("qMatrix diagonal must be non-positive at row " + row);
                    }
                } else {
                    if (value < -Q_TOLERANCE) {
                        throw new IllegalArgumentException("qMatrix off-diagonal must be non-negative at [" +
                                row + "," + col + "]");
                    }
                    offDiagonalSum += Math.max(0.0, value);
                }
            }
            if (Math.abs(rowSum) > Q_TOLERANCE) {
                throw new IllegalArgumentException("qMatrix row " + row + " must sum to zero; found " + rowSum);
            }
            final double diagonal = qMatrix[row * stateCount + row];
            if (Math.abs(diagonal + offDiagonalSum) > Q_TOLERANCE) {
                throw new IllegalArgumentException("qMatrix diagonal at row " + row +
                        " must equal negative off-diagonal sum.");
            }
        }
    }

    private static int validateRootState(final int rootState, final int stateCount) {
        if (rootState < 0 || rootState >= stateCount) {
            throw new IllegalArgumentException("rootState out of range: " + rootState);
        }
        return rootState;
    }

    private static double[] normalizeFrequencies(final double[] frequencies, final int stateCount) {
        if (frequencies.length != stateCount) {
            throw new IllegalArgumentException("rootFrequencies must have length stateCount=" + stateCount);
        }
        final double[] out = Arrays.copyOf(frequencies, frequencies.length);
        double sum = 0.0;
        for (int i = 0; i < out.length; i++) {
            if (!Double.isFinite(out[i]) || out[i] < 0.0) {
                throw new IllegalArgumentException("Invalid root frequency at state " + i + ": " + out[i]);
            }
            sum += out[i];
        }
        if (!(sum > 0.0)) {
            throw new IllegalArgumentException("rootFrequencies must have positive sum");
        }
        for (int i = 0; i < out.length; i++) {
            out[i] /= sum;
        }
        return out;
    }

    private static double[] pointMass(final int state, final int stateCount) {
        final double[] frequencies = new double[stateCount];
        frequencies[state] = 1.0;
        return frequencies;
    }

    private static double[] stationaryFrequencies(final double[] qMatrix, final int stateCount) {
        final double[][] a = new double[stateCount][stateCount];
        final double[] b = new double[stateCount];

        for (int row = 0; row < stateCount - 1; row++) {
            for (int col = 0; col < stateCount; col++) {
                a[row][col] = qMatrix[col * stateCount + row];
            }
        }
        for (int col = 0; col < stateCount; col++) {
            a[stateCount - 1][col] = 1.0;
        }
        b[stateCount - 1] = 1.0;

        final double[] pi = solveLinearSystem(a, b);
        for (int i = 0; i < pi.length; i++) {
            if (pi[i] < 0.0 && pi[i] > -1.0e-10) {
                pi[i] = 0.0;
            }
            if (!Double.isFinite(pi[i]) || pi[i] < 0.0) {
                throw new IllegalArgumentException(
                        "Could not derive a valid stationary root distribution from qMatrix: " +
                                Arrays.toString(pi)
                );
            }
        }
        return normalizeFrequencies(pi, stateCount);
    }

    private static double[] solveLinearSystem(final double[][] a, final double[] b) {
        final int n = b.length;
        final double[][] m = new double[n][n + 1];
        for (int i = 0; i < n; i++) {
            System.arraycopy(a[i], 0, m[i], 0, n);
            m[i][n] = b[i];
        }

        for (int col = 0; col < n; col++) {
            int pivot = col;
            double pivotAbs = Math.abs(m[col][col]);
            for (int row = col + 1; row < n; row++) {
                final double abs = Math.abs(m[row][col]);
                if (abs > pivotAbs) {
                    pivot = row;
                    pivotAbs = abs;
                }
            }
            if (pivotAbs < 1.0e-12) {
                throw new IllegalArgumentException("Singular linear system while solving stationary distribution");
            }
            if (pivot != col) {
                final double[] tmp = m[col];
                m[col] = m[pivot];
                m[pivot] = tmp;
            }

            final double invPivot = 1.0 / m[col][col];
            for (int j = col; j <= n; j++) {
                m[col][j] *= invPivot;
            }
            for (int row = 0; row < n; row++) {
                if (row == col) {
                    continue;
                }
                final double factor = m[row][col];
                if (factor == 0.0) {
                    continue;
                }
                for (int j = col; j <= n; j++) {
                    m[row][j] -= factor * m[col][j];
                }
            }
        }

        final double[] x = new double[n];
        for (int i = 0; i < n; i++) {
            x[i] = m[i][n];
        }
        return x;
    }

    private static double[][] copy(final double[][] values) {
        final double[][] out = new double[values.length][];
        for (int i = 0; i < values.length; i++) {
            out[i] = Arrays.copyOf(values[i], values[i].length);
        }
        return out;
    }

    private static double[] flatten(final double[][] values) {
        int length = 0;
        for (double[] row : values) {
            length += row.length;
        }
        final double[] out = new double[length];
        int offset = 0;
        for (double[] row : values) {
            System.arraycopy(row, 0, out, offset, row.length);
            offset += row.length;
        }
        return out;
    }

    private static double[] toDouble(final int[] values) {
        final double[] out = new double[values.length];
        for (int i = 0; i < values.length; i++) {
            out[i] = values[i];
        }
        return out;
    }

    public static final class Summary {
        public final int totalJumps;
        public final int atomicBranchCount;
        public final int continuousBranchCount;
        public final double totalTreeTime;
        public final double totalReward;
        public final double meanRewardProportion;
        public final double minRewardProportion;
        public final double maxRewardProportion;
        public final double[] dwellTimesByState;
        public final int[] transitionCounts;
        public final int[] tipStateCounts;

        private Summary(final int totalJumps,
                        final int atomicBranchCount,
                        final int continuousBranchCount,
                        final double totalTreeTime,
                        final double totalReward,
                        final double meanRewardProportion,
                        final double minRewardProportion,
                        final double maxRewardProportion,
                        final double[] dwellTimesByState,
                        final int[] transitionCounts,
                        final int[] tipStateCounts) {
            this.totalJumps = totalJumps;
            this.atomicBranchCount = atomicBranchCount;
            this.continuousBranchCount = continuousBranchCount;
            this.totalTreeTime = totalTreeTime;
            this.totalReward = totalReward;
            this.meanRewardProportion = meanRewardProportion;
            this.minRewardProportion = minRewardProportion;
            this.maxRewardProportion = maxRewardProportion;
            this.dwellTimesByState = Arrays.copyOf(dwellTimesByState, dwellTimesByState.length);
            this.transitionCounts = Arrays.copyOf(transitionCounts, transitionCounts.length);
            this.tipStateCounts = Arrays.copyOf(tipStateCounts, tipStateCounts.length);
        }

        @Override
        public String toString() {
            return "Summary{" +
                    "totalJumps=" + totalJumps +
                    ", atomicBranchCount=" + atomicBranchCount +
                    ", continuousBranchCount=" + continuousBranchCount +
                    ", totalTreeTime=" + totalTreeTime +
                    ", totalReward=" + totalReward +
                    ", meanRewardProportion=" + meanRewardProportion +
                    ", minRewardProportion=" + minRewardProportion +
                    ", maxRewardProportion=" + maxRewardProportion +
                    ", dwellTimesByState=" + Arrays.toString(dwellTimesByState) +
                    ", transitionCounts=" + Arrays.toString(transitionCounts) +
                    ", tipStateCounts=" + Arrays.toString(tipStateCounts) +
                    '}';
        }
    }
}
