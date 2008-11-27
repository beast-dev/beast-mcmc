package dr.evomodel.treelikelihood;

import dr.evolution.alignment.PatternList;
import dr.evolution.tree.Tree;
import dr.evolution.util.TaxonList;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id$
 */
public abstract class TipPartialsModel extends AbstractModel {

    /**
     * @param name Model Name
     */
    public TipPartialsModel(String name, TaxonList includeTaxa, TaxonList excludeTaxa) {
        super(name);

        this.includeTaxa = includeTaxa;
        this.excludeTaxa = excludeTaxa;
    }

    public final void setTree(Tree tree) {
        this.tree = tree;

        int extNodeCount = tree.getExternalNodeCount();

        excluded = new boolean[extNodeCount];
        if (includeTaxa != null) {
            for (int i = 0; i < extNodeCount; i++) {
                if (includeTaxa.getTaxonIndex(tree.getNodeTaxon(tree.getExternalNode(i))) == -1) {
                    excluded[i] = true;
                }
            }
        }

        if (excludeTaxa != null) {
            for (int i = 0; i < extNodeCount; i++) {
                if (excludeTaxa.getTaxonIndex(tree.getNodeTaxon(tree.getExternalNode(i))) != -1) {
                    excluded[i] = true;
                }
            }

        }

        states = new int[extNodeCount][];
    }

    public final void setStates(PatternList patternList, int sequenceIndex, int nodeIndex) {
        if (patternCount == 0) {
            patternCount = patternList.getPatternCount();
            stateCount = patternList.getDataType().getStateCount();
        }
        if (this.states[nodeIndex] == null) {
            this.states[nodeIndex] = new int[patternCount];
        }

        for (int i = 0; i < patternCount; i++) {
            this.states[nodeIndex][i] = patternList.getPatternState(sequenceIndex, i);
        }
    }

    protected void handleModelChangedEvent(Model model, Object object, int index) {
        fireModelChanged();
    }

    /**
     * This method is called whenever a parameter is changed.
     * <p/>
     * It is strongly recommended that the model component sets a "dirty" flag and does no
     * further calculations. Recalculation is typically done when the model component is asked for
     * some information that requires them. This mechanism is 'lazy' so that this method
     * can be safely called multiple times with minimal computational cost.
     */
    protected final void handleParameterChangedEvent(Parameter parameter, int index, Parameter.ChangeType type) {
        fireModelChanged();
    }

    /**
     * Additional state information, outside of the sub-model is stored by this call.
     */
    protected void storeState() {
    }

    /**
     * After this call the model is guaranteed to have returned its extra state information to
     * the values coinciding with the last storeState call.
     * Sub-models are handled automatically and do not need to be considered in this method.
     */
    protected void restoreState() {
    }

    /**
     * This call specifies that the current state is accept. Most models will not need to do anything.
     * Sub-models are handled automatically and do not need to be considered in this method.
     */
    protected void acceptState() {
    }

    public abstract void getTipPartials(int nodeIndex, double[] tipPartials);

    protected int[][] states;
    protected boolean[] excluded;

    protected int patternCount = 0;
    protected int stateCount;

    protected TaxonList includeTaxa;
    protected TaxonList excludeTaxa;

    protected Tree tree;
}
