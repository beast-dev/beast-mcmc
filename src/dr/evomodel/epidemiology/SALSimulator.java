package dr.evomodel.epidemiology;

public class SALSimulator extends TauLeapingSimulator {

    public SALSimulator(CompartmentalModel compartmentalModel, double epsilon, int criticalNumber) {
        super(compartmentalModel, epsilon, criticalNumber);
    }

    protected double[] getPoissonIntensities(double[] currentCounts, double[] reactionInt, double tau) {
        return compartmentalModel.getSALPoissonIntensities(currentCounts, reactionInt, tau);
    }
}
