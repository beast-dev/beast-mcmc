package dr.evomodel.epidemiology;

import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;

import java.util.List;

public abstract class CompartmentalModel extends AbstractModel {

    protected List<Parameter> rateParameters;
    protected List<Parameter> compartmentCounts;
    protected Parameter originOne;
    protected Parameter originTwo;
    protected int numGridPoints;
    protected double cutOff;
    protected int numReactionChannels;
    // number of different kinds of species/particle types
    protected int numSpecies;
    // v matrix describes how count vector changes with reaction
    // row corresponds to species/particle type, column corresponds to reaction channel
    protected int[][] vMatrix;

    public CompartmentalModel(String name) {super(name);}

    protected abstract void setOriginTimeCompartmentCounts(int index);

    protected abstract int[] getHighestOrdersOfReactions();

    protected abstract void setDefaultCompartmentCounts(int index);

    protected abstract double[] getReactionIntensities(double[] compartmentCounts);

    //protected abstract double[] getMaxFiringTimes(double[] compartmentCounts, double[] reactionIntensities);

    //protected abstract double[] getTauLeapingPoissonIntensities(double[] currentCounts, double[] reactionInt, double tau);

    protected abstract double[] getSALPoissonIntensities(double[] currentCounts, double[] reactionInt, double tau);

    protected abstract double[] getUpdatedCompartmentCounts(double[] currentCounts, double[] countsNew);

    // Checks of all counts are nonnegative
    protected boolean hasMinimalCounts(double[] counts) {
        for (int i = 0; i < counts.length; i++) {
            if (counts[i] < 0) {
                return false;
            }
        }
        return true;
    }

    protected double getOldestOrigin() {
        return originOne.getParameterValue(0);
    }

    public double[] introduceSecondPathogen(
            //double previousTime,
            double simulationTime,
            double[] currentCounts) {
        return currentCounts;
    }

    protected void handleModelChangedEvent(Model model, Object object, int index) {
        // no intermediates need to be recalculated...
    }

    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        // no intermediates need to be recalculated...
    }

    protected void storeState() {
    } // no additional state needs storing

    protected void restoreState() {
    } // no additional state needs restoring

    protected void acceptState() {
    } // no additional state needs accepting
}
