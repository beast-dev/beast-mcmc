package dr.inference.model;

import dr.inference.loggers.LogColumn;
import dr.inferencexml.model.DefaultModelParser;

/**
 * @author Marc Suchard
 */
public class DefaultModel extends AbstractModelLikelihood {

    public DefaultModel() {
        super(DefaultModelParser.DUMMY_MODEL);
    }

    public DefaultModel(Parameter parameter) {
        super(DefaultModelParser.DUMMY_MODEL);
        addVariable(parameter);
    }

    protected void handleModelChangedEvent(Model model, Object object, int index) {

    }

    protected final void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {

    }

    protected void storeState() {

    }

    protected void restoreState() {

    }

    protected void acceptState() {

    }

    public Model getModel() {
        return this;
    }

    public double getLogLikelihood() {
        return 0;
    }

    public void makeDirty() {

    }

    public LogColumn[] getColumns() {
        return new LogColumn[0];
    }

}


