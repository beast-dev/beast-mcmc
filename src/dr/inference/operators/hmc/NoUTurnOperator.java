package dr.inference.operators.hmc;

import dr.inference.hmc.ReversibleHMCProvider;
import dr.inference.operators.GibbsOperator;
import dr.inference.operators.SimpleMCMCOperator;
import dr.math.MathUtils;
import dr.math.matrixAlgebra.WrappedVector;

import java.util.Arrays;

public class NoUTurnOperator extends SimpleMCMCOperator implements GibbsOperator {

    class Options {
        private double logProbErrorTol = 100.0;
        private int findMax = 100;
        private int maxHeight = 10;
    }

    private final Options options = new Options();

    public NoUTurnOperator(ReversibleHMCProvider hmcProvider,
                           boolean adaptiveStepsize,
                           double weight) {

        this.hmcProvider = hmcProvider;
        this.adaptiveStepsize = adaptiveStepsize;
        setWeight(weight);
    }

    @Override
    public String getOperatorName() {
        return "No-U-Turn Operator";
    }

    @Override
    public double doOperation() {

        final double[] initialPosition = hmcProvider.getInitialPosition();

        if (stepSizeInformation == null) {
            stepSizeInformation = findReasonableStepSize(initialPosition,
                    hmcProvider.getGradientProvider().getGradientLogDensity(), hmcProvider.getStepSize());
        }

        double[] position = takeOneStep(getCount() + 1, initialPosition);

        hmcProvider.setParameter(position);
        return 0;
    }

    private double[] takeOneStep(long m, double[] initialPosition) {

        double[] endPosition = Arrays.copyOf(initialPosition, initialPosition.length);
        final WrappedVector initialMomentum = hmcProvider.drawMomentum();

        final double initialJointDensity = hmcProvider.getJointProbability(initialMomentum);

        double logSliceU = Math.log(MathUtils.nextDouble()) + initialJointDensity;

        TreeState trajectoryTree = new TreeState(initialPosition, initialMomentum.getBuffer(),
                hmcProvider.getGradientProvider().getGradientLogDensity(), 1, true);

        int height = 0;

        while (trajectoryTree.flagContinue) {
            double[] tmp = updateTrajectoryTree(trajectoryTree, height, logSliceU, initialJointDensity);
            if (tmp != null) {
                endPosition = tmp;
            }

            height++;

            if (height > options.maxHeight) {
                trajectoryTree.flagContinue = false;
            }
        }
        if (adaptiveStepsize) {
            stepSizeInformation.update(m, trajectoryTree.cumAcceptProb, trajectoryTree.numAcceptProbStates);
        }
        return endPosition;
    }

    private double[] updateTrajectoryTree(TreeState trajectoryTree, int depth, double logSliceU,
                                          double initialJointDensity) {

        double[] endPosition = null;

        final double uniform1 = MathUtils.nextDouble();
        int direction = (uniform1 < 0.5) ? -1 : 1;
        TreeState nextTrajectoryTree = buildTree(
                trajectoryTree.getPosition(direction), trajectoryTree.getMomentum(direction),
                trajectoryTree.getGradient(direction),
                direction, logSliceU, depth, stepSizeInformation.getStepSize(), initialJointDensity);

        if (nextTrajectoryTree.flagContinue) {

            final double uniform = MathUtils.nextDouble();
            final double acceptProb = (double) nextTrajectoryTree.numNodes / (double) trajectoryTree.numNodes;
            if (uniform < acceptProb) {
                endPosition = nextTrajectoryTree.getSample();
            }
        }

        trajectoryTree.mergeNextTree(nextTrajectoryTree, direction);

        return endPosition;
    }

    private TreeState buildTree(double[] position, double[] momentum, double[] gradient, int direction,
                                double logSliceU, int height, double stepSize, double initialJointDensity) {

        if (height == 0) {
            return buildBaseCase(position, momentum, gradient, direction, logSliceU, stepSize, initialJointDensity);
        } else {
            return buildRecursiveCase(position, momentum, gradient, direction, logSliceU, height, stepSize,
                    initialJointDensity);
        }
    }


    private TreeState buildBaseCase(double[] inPosition, double[] inMomentum, double[] inGradient, int direction,
                                    double logSliceU, double stepSize, double initialJointDensity) {

        // Make deep copy of position and momentum
        WrappedVector position = new WrappedVector.Raw(Arrays.copyOf(inPosition, inPosition.length));
        WrappedVector momentum = new WrappedVector.Raw(Arrays.copyOf(inMomentum, inMomentum.length));
        WrappedVector gradient = new WrappedVector.Raw(Arrays.copyOf(inGradient, inGradient.length));

        hmcProvider.setParameter(position.getBuffer());

        // "one reversibleHMC integral
        hmcProvider.reversiblePositionMomentumUpdate(position, momentum, gradient, direction, stepSize);

        double logJointProbAfter = hmcProvider.getJointProbability(momentum);

        final int numNodes = (logSliceU <= logJointProbAfter ? 1 : 0);

        final boolean flagContinue = (logSliceU < options.logProbErrorTol + logJointProbAfter);

        // Values for dual-averaging
        final double acceptProb = Math.min(1.0, Math.exp(logJointProbAfter - initialJointDensity));
        final int numAcceptProbStates = 1;

        hmcProvider.setParameter(inPosition);

        return new TreeState(position.getBuffer(), momentum.getBuffer(), gradient.getBuffer(), numNodes, flagContinue
                , acceptProb,
                numAcceptProbStates);
    }

    private TreeState buildRecursiveCase(double[] inPosition, double[] inMomentum, double[] gradient, int direction,
                                         double logSliceU, int height, double stepSize, double initialJointDensity) {

        TreeState subtree = buildTree(inPosition, inMomentum, gradient, direction, logSliceU,
                height - 1, // Recursion
                stepSize, initialJointDensity);

        if (subtree.flagContinue) {

            TreeState nextSubtree = buildTree(subtree.getPosition(direction), subtree.getMomentum(direction),
                    subtree.getGradient(direction), direction,
                    logSliceU, height - 1, stepSizeInformation.getStepSize(), initialJointDensity);

            subtree.mergeNextTree(nextSubtree, direction);

        }
        return subtree;
    }

    private static boolean computeStopCriterion(boolean flagContinue, TreeState state) {
        return computeStopCriterion(flagContinue,
                state.getPosition(1), state.getPosition(-1),
                state.getMomentum(1), state.getMomentum(-1));
    }

    private StepSize findReasonableStepSize(double[] initialPosition, double[] initialGradient,
                                            double forcedInitialStepSize) {

        if (forcedInitialStepSize != 0) {
            return new StepSize(forcedInitialStepSize);
        } else {
            double stepSize = 0.1;

            WrappedVector momentum = hmcProvider.drawMomentum();
            int count = 1;
            int dim = initialPosition.length;
            WrappedVector position = new WrappedVector.Raw(Arrays.copyOf(initialPosition, dim));
            WrappedVector gradient = new WrappedVector.Raw(Arrays.copyOf(initialGradient, dim));

            double probBefore = hmcProvider.getJointProbability(momentum);

            hmcProvider.reversiblePositionMomentumUpdate(position, momentum, gradient, 1, stepSize);

            double probAfter = hmcProvider.getJointProbability(momentum);

            double a = ((probAfter - probBefore) > Math.log(0.5) ? 1 : -1);

            double probRatio = Math.exp(probAfter - probBefore);

            while (Math.pow(probRatio, a) > Math.pow(2, -a)) {

                probBefore = probAfter;
                hmcProvider.reversiblePositionMomentumUpdate(position, momentum, gradient, 1, stepSize);

                probAfter = hmcProvider.getJointProbability(momentum);
                probRatio = Math.exp(probAfter - probBefore);

                stepSize = Math.pow(2, a) * stepSize;
                count++;

                if (count > options.findMax) {
                    throw new RuntimeException("Cannot find a reasonable step-size in " + options.findMax + " " +
                            "iterations");
                }
            }
            hmcProvider.setParameter(initialPosition);
            return new StepSize(stepSize);
        }
    }


    private static boolean computeStopCriterion(boolean flagContinue,
                                                double[] positionPlus, double[] positionMinus,
                                                double[] momentumPlus, double[] momentumMinus) {

        double[] positionDifference = subtractArray(positionPlus, positionMinus);

        return flagContinue &&
                getDotProduct(positionDifference, momentumMinus) >= 0 &&
                getDotProduct(positionDifference, momentumPlus) >= 0;
    }

    private static double getDotProduct(double[] x, double[] y) {

        assert (x.length == y.length);
        final int dim = x.length;

        double total = 0.0;
        for (int i = 0; i < dim; i++) {
            total += x[i] * y[i];
        }
        return total;
    }

    private static double[] subtractArray(double[] a, double[] b) {

        assert (a.length == b.length);
        final int dim = a.length;

        double[] result = new double[dim];
        for (int i = 0; i < dim; i++) {
            result[i] = a[i] - b[i];
        }

        return result;
    }

    private class TreeState {

        private TreeState(double[] position, double[] moment, double[] gradient,
                          int numNodes, boolean flagContinue) {
            this(position, moment, gradient, numNodes, flagContinue, 0.0, 0);
        }

        private TreeState(double[] position, double[] moment, double[] gradient,
                          int numNodes, boolean flagContinue,
                          double cumAcceptProb, int numAcceptProbStates) {
            this.position = new double[3][];
            this.momentum = new double[3][];
            this.gradient = new double[3][]; //todo: (for gradient) no need for 3 but 2? If changed to 2, getIndex should also be changed

            for (int i = 0; i < 3; ++i) {
                this.position[i] = position;
                this.momentum[i] = moment;
                this.gradient[i] = gradient;
            }

            // Recursion variables
            this.numNodes = numNodes;
            this.flagContinue = flagContinue;

            // Dual-averaging variables
            this.cumAcceptProb = cumAcceptProb;
            this.numAcceptProbStates = numAcceptProbStates;
        }

        private double[] getPosition(int direction) {
            return position[getIndex(direction)];
        }

        private double[] getMomentum(int direction) {
            return momentum[getIndex(direction)];
        }

        private double[] getGradient(int direction) {
            return gradient[getIndex(direction)];
        }

        private double[] getSample() {
            /*
            Returns a state chosen uniformly from the acceptable states along a hamiltonian dynamics trajectory tree.
            The sample is updated recursively while building trees.
            */
            return position[getIndex(0)];
        }

        private void setPosition(int direction, double[] position) {
            this.position[getIndex(direction)] = position;
        }

        private void setMomentum(int direction, double[] momentum) {
            this.momentum[getIndex(direction)] = momentum;
        }

        private void setGradient(int direction, double[] gradient) {
            this.gradient[getIndex(direction)] = gradient;
        }

        private void setSample(double[] position) { setPosition(0, position); }

        private int getIndex(int direction) { // valid directions: -1, 0, +1
            assert (direction >= -1 && direction <= 1);
            return direction + 1;
        }

        private void mergeNextTree(TreeState nextTree, int direction) {

            setPosition(direction, nextTree.getPosition(direction));
            setMomentum(direction, nextTree.getMomentum(direction));
            setGradient(direction, nextTree.getGradient(direction));

            updateSample(nextTree);

            numNodes += nextTree.numNodes;
            flagContinue = computeStopCriterion(nextTree.flagContinue, this);

            cumAcceptProb += nextTree.cumAcceptProb;
            numAcceptProbStates += nextTree.numAcceptProbStates;
        }

        private void updateSample(TreeState nextTree) {
            double uniform = MathUtils.nextDouble();
            if (nextTree.numNodes > 0
                    && uniform < ((double) nextTree.numNodes / (double) (numNodes + nextTree.numNodes))) {
                setSample(nextTree.getSample());
            }
        }

        final private double[][] position;
        final private double[][] momentum;
        final private double[][] gradient;

        private int numNodes;
        private boolean flagContinue;

        private double cumAcceptProb;
        private int numAcceptProbStates;
    }

    private ReversibleHMCProvider hmcProvider;
    private StepSize stepSizeInformation;
    private boolean adaptiveStepsize;
}

