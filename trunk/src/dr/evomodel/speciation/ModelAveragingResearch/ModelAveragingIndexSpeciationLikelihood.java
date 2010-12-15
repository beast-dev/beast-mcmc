package dr.evomodel.speciation.ModelAveragingResearch;

import dr.evolution.util.Units;
import dr.inference.model.*;

/**
 * @author Walter Xie
 */
public class ModelAveragingIndexSpeciationLikelihood extends AbstractModelLikelihood implements Units {


    public ModelAveragingIndexSpeciationLikelihood(String name, Variable<Integer> indexVariable, Variable<Double> maxIndexVariable) {

        super(name);

        this.indexVariable = indexVariable; // integer index parameter size = real size - 1
        for (int i=0; i<indexVariable.getSize(); i++) {
            indexVariable.setValue(i, 0);
        }
        indexVariable.addBounds(new Bounds.Staircase(indexVariable));
        addVariable(indexVariable);

        for (int i=0; i<maxIndexVariable.getSize(); i++) {
            maxIndexVariable.setValue(i, 0.0);
        }
        this.maxIndexVariable = maxIndexVariable;
        addVariable(maxIndexVariable);
    }

    // **************************************************************
    // ModelListener IMPLEMENTATION
    // **************************************************************

    protected final void handleModelChangedEvent(Model model, Object object, int index) {
        likelihoodKnown = false;
    }

    // **************************************************************
    // VariableListener IMPLEMENTATION
    // **************************************************************

    protected final void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        likelihoodKnown = false;
    } // No parameters to respond to

    // **************************************************************
    // Model IMPLEMENTATION
    // **************************************************************

    /**
     * Stores the precalculated state: likelihood
     */
    protected final void storeState() {
        storedLikelihoodKnown = likelihoodKnown;
        storedLogLikelihood = logLikelihood;
    }

    /**
     * Restores the precalculated state: computed likelihood
     */
    protected final void restoreState() {
        likelihoodKnown = storedLikelihoodKnown;
        logLikelihood = storedLogLikelihood;
    }

    protected final void acceptState() {
    } // nothing to do

    // **************************************************************
    // Likelihood IMPLEMENTATION
    // **************************************************************

    public final Model getModel() {
        return this;
    }

    public final double getLogLikelihood() {
        if (!likelihoodKnown) {
            logLikelihood = calculateLogLikelihood();
            likelihoodKnown = true;
        }
        return logLikelihood;
    }

    public final void makeDirty() {
        likelihoodKnown = false;
    }

    /**
     * Calculates the log likelihood of this set of coalescent intervals,
     * given a demographic model.
     *
     * @return the log likelihood
     */
    private double calculateLogLikelihood() {

        // Rule: index k cannot be appeared unless k-1 appeared

            if (!isValidate(indexVariable.getValues())) {
//                output("illegal index variable", indexVariable);
                return Double.NEGATIVE_INFINITY;
            }

        Double maxI = (double)(int) getMaxIndex(indexVariable.getValues());
        maxIndexVariable.setValue(0, maxI);

        return 0.0;
    }

    private boolean isValidate(Integer[] pattern) {
        // Rule: index k cannot be appeared unless k-1 appeared before it appears
        int[] indexFreq = new int[pattern.length];
        for (int i = 0; i < pattern.length; i++) {
            if (pattern[i] > 0) // not validate 0
                indexFreq[pattern[i] - 1] += 1; // integer index parameter size = real size - 1

            if (i > 0 && (pattern[i] - pattern[i - 1] > 1)) {
                for (int f = 0; f < i; f++) {
                    if (indexFreq[f] < 1) return false;
                }
            }
        }

        return true;
    }

    private int getMaxIndex(Integer[] pattern) {
        int max = 0;

        for (int p : pattern) {
            if (p > max) {
                max = p;
            }
        }

        return max;
    }

    private void output(String message, Variable<Integer> indexVariable) {
        System.out.print(message + ": ");
        for (int i = 0; i < indexVariable.getSize(); i++) {
            System.out.print(indexVariable.getValue(i) + "\t");
        }
        System.out.println();
    }

    // **************************************************************
    // Loggable IMPLEMENTATION
    // **************************************************************

    /**
     * @return the log columns.
     */
    public final dr.inference.loggers.LogColumn[] getColumns() {

        String columnName = getId();
        if (columnName == null) columnName = getModelName() + ".likelihood";

        return new dr.inference.loggers.LogColumn[]{
                new LikelihoodColumn(columnName)
        };
    }

    private final class LikelihoodColumn extends dr.inference.loggers.NumberColumn {
        public LikelihoodColumn(String label) {
            super(label);
        }

        public double getDoubleValue() {
            return getLogLikelihood();
        }
    }

    // **************************************************************
    // Units IMPLEMENTATION
    // **************************************************************

    /**
     * Sets the units these coalescent intervals are
     * measured in.
     */
    public final void setUnits(Type u) {

    }

    /**
     * Returns the units these coalescent intervals are
     * measured in.
     */
    public final Type getUnits() {
        return Type.YEARS;
    }

    // ****************************************************************
    // Private and protected stuff
    // ****************************************************************

    Variable<Double> maxIndexVariable = null;

    Variable<Integer> indexVariable = null;

    private double logLikelihood;
    private double storedLogLikelihood;
    private boolean likelihoodKnown = false;
    private boolean storedLikelihoodKnown = false;
}