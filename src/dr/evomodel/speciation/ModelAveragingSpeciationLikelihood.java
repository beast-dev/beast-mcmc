package dr.evomodel.speciation;

import dr.evolution.tree.Tree;
import dr.evolution.util.Units;
import dr.evomodelxml.speciation.SpeciationLikelihoodParser;
import dr.inference.model.*;

import java.util.List;

/**
 * @author Alexei Drummond
 */
public class ModelAveragingSpeciationLikelihood extends AbstractModelLikelihood implements Units {

    // PUBLIC STUFF
    /**
     * @param trees            the tree
     * @param speciationModels the model of speciation
     * @param id               a unique identifier for this likelihood
     */
    public ModelAveragingSpeciationLikelihood(List<Tree> trees, List<MaskableSpeciationModel> speciationModels, Variable<Integer> indexVariable, String id) {
        this(SpeciationLikelihoodParser.SPECIATION_LIKELIHOOD, trees, speciationModels, indexVariable);
        setId(id);
    }

    public ModelAveragingSpeciationLikelihood(String name, List<Tree> trees, List<MaskableSpeciationModel> speciationModels, Variable<Integer> indexParameter) {

        super(name);

        this.trees = trees;
        this.speciationModels = speciationModels;

        if (trees.size() != speciationModels.size()) {
            throw new IllegalArgumentException("The number of trees and the number of speciation models should be equal.");
        }

        for (Tree tree : trees) {
            if (tree instanceof Model) {
                addModel((Model) tree);
            }
        }
        for (SpeciationModel speciationModel : speciationModels) {
            if (speciationModel != null) {
                addModel(speciationModel);
            }
        }


        if (indexParameter.getSize() != trees.size()) {
            throw new IllegalArgumentException("Index parameter must be same size as the number of trees.");
        }
        indexParameter.addBounds(new Bounds.Staircase(indexParameter));
        addVariable(indexParameter);

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

        double logL = 0;
        for (int i = 0; i < trees.size(); i++) {

            MaskableSpeciationModel model = speciationModels.get(i);
            SpeciationModel mask = speciationModels.get(indexVariable.getValue(i));
            if (model != mask) {
                model.mask(mask);
            } else {
                model.unmask();
            }

            logL += model.calculateTreeLogLikelihood(trees.get(i));
        }

        return logL;
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
        for (SpeciationModel speciationModel : speciationModels) {
            speciationModel.setUnits(u);
        }
    }

    /**
     * Returns the units these coalescent intervals are
     * measured in.
     */
    public final Type getUnits() {
        return speciationModels.get(0).getUnits();
    }

    // ****************************************************************
    // Private and protected stuff
    // ****************************************************************

    /**
     * The speciation models.
     */
    List<MaskableSpeciationModel> speciationModels = null;

    /**
     * The trees.
     */
    List<Tree> trees = null;

    Variable<Integer> indexVariable = null;

    private double logLikelihood;
    private double storedLogLikelihood;
    private boolean likelihoodKnown = false;
    private boolean storedLikelihoodKnown = false;
}