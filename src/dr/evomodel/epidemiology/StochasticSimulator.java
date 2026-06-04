package dr.evomodel.epidemiology;

import dr.math.MathUtils;

public abstract class StochasticSimulator {

    protected CompartmentalModel compartmentalModel;
    protected double cutOff;
    protected int numGridPoints;
    protected int numSpecies;
    protected int numReactionChannels;
    protected int[][] vMatrix;
    protected double intervalWidth;

    public StochasticSimulator(CompartmentalModel compartmentalModel) {
        this.compartmentalModel = compartmentalModel;
        this.cutOff = compartmentalModel.cutOff;
        this.numGridPoints = compartmentalModel.numGridPoints;
        this.numSpecies = compartmentalModel.numSpecies;
        this.numReactionChannels = compartmentalModel.numReactionChannels;
        this.vMatrix = compartmentalModel.vMatrix;
        this.intervalWidth = cutOff/numGridPoints;
    }

    public abstract void simulateTrajectory();

    protected int sampleReactionChannel(double[] reactionInt, double reactionIntSum) {
        double r = MathUtils.nextDouble();
        double threshold = r*reactionIntSum;
        double cumulative = 0.0;

        for (int i = 0; i < reactionInt.length; i++) {
            cumulative = cumulative + reactionInt[i];
            if (threshold < cumulative) {
                return i;
            }
        }
        return reactionInt.length - 1;
    }

}
