package dr.evomodel.continuous;

import dr.evomodel.tree.TreeModel;
import dr.inference.model.AbstractModelLikelihood;
import dr.inference.model.CompoundParameter;
import dr.inference.model.Model;
import dr.inference.model.Variable;

import java.util.List;

/**
 * @author Andrew Rambaut
 * @author Marc Suchard
 * @version $Id$
 */
public class AntigenicTraitLikelihood extends AbstractModelLikelihood {

    public final static String ANTIGENIC_TRAIT_LIKELIHOOD = "antigenicTraitLikelihood";

    public AntigenicTraitLikelihood(
            String traitName, TreeModel treeModel,
            MultivariateDiffusionModel diffusionModel,
            CompoundParameter tipTraitParameter,
            CompoundParameter virusTraitParameter,
            CompoundParameter serumTraitParameter,
            double[][] antigenicAssayTable,
            String[] virusNames,
            String[] seraNames) {

        super(ANTIGENIC_TRAIT_LIKELIHOOD);
        int tipCount = tipTraitParameter.getDimension();

        // the total number of viruses is the number of rows in the table
        int virusCount = antigenicAssayTable.length;
        // the number of sera is the number of columns
        int serumCount = antigenicAssayTable[0].length;

    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Variable.ChangeType type) {
    }

    @Override
    protected void storeState() {
    }

    @Override
    protected void restoreState() {
    }

    @Override
    protected void acceptState() {
    }

    public Model getModel() {
        return null;
    }

    public double getLogLikelihood() {
        return 0;
    }

    public void makeDirty() {
    }
}
