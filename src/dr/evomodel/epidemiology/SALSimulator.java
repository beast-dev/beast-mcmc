package dr.evomodel.epidemiology;

import dr.inference.model.Parameter;

public class SALSimulator extends TauLeapingSimulator {

    public SALSimulator(CompartmentalModel compartmentalModel, double epsilon, int criticalNumber) {
        super(compartmentalModel, epsilon, criticalNumber);
    }

    public SALSimulator(CompartmentalModel compartmentalModel, double epsilon, int criticalNumber, Parameter elapsedTime) {
        this(compartmentalModel, epsilon, criticalNumber);
        this.elapsedTime = elapsedTime;
    }

    protected double[] getPoissonIntensities(double[] currentCounts, double[] reactionInt, double tau, double simTime) {
        return compartmentalModel.getSALPoissonIntensities(currentCounts, reactionInt, tau, simTime);
    }
}
