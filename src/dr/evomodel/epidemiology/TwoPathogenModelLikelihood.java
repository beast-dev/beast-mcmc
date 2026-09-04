package dr.evomodel.epidemiology;

import dr.inference.model.AbstractModelLikelihood;
import dr.inference.model.Model;
import dr.inference.model.Variable;

public class TwoPathogenModelLikelihood extends AbstractModelLikelihood {

    private final TwoPathogenModel twoPathogenModel;
    private boolean likelihoodKnown = false;

    public TwoPathogenModelLikelihood(TwoPathogenModel twoPathogenModel) {

        super("TwoPathogenModelLikelihood");

        this.twoPathogenModel = twoPathogenModel;
        addModel(twoPathogenModel);
    }

    public double getLogLikelihood() {
        // always return 0.0 since trajectory prior and proposal cancel in MH ratio
        return 0.0;
    }

    public void makeDirty(){
        // do nothing
    }

    protected void handleModelChangedEvent(Model model, Object object, int index){
        fireModelChanged();
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Variable.ChangeType type) {
        fireModelChanged();
    }

    protected void storeState(){
        // do nothing
    }

    protected void restoreState(){
        // do nothing
    }

    protected void acceptState(){
        // do nothing
    }

    public Model getModel(){
        return this;
    }

}
