package dr.evomodel.treelikelihood;

import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id$
 */
public abstract class TipPartialsModel extends AbstractModel {

    /**
     * @param name Model Name
     */
    public TipPartialsModel(String name) {
        super(name);
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
    protected void handleParameterChangedEvent(Parameter parameter, int index) {
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

    public abstract double[] getTipPartials(int nodeIndex);
 
}
