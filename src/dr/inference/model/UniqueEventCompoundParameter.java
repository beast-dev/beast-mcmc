package dr.inference.model;

/**
 * Temporary HMC-focused compound parameter that propagates child change events
 * once per unique child parameter rather than once per dimension entry.
 */
public class UniqueEventCompoundParameter extends CompoundParameter {

    public UniqueEventCompoundParameter(String name) {
        super(name);
    }

    public UniqueEventCompoundParameter(String name, Parameter[] params) {
        super(name, params);
    }

    public UniqueEventCompoundParameter(CompoundParameter source) {
        super(source.getParameterName());
        if (source.getId() != null) {
            setId(source.getId());
        }
        for (int i = 0; i < source.getParameterCount(); ++i) {
            addParameter(source.getParameter(i));
        }
    }

    @Override
    public void fireParameterChangedEvent() {
        doNotPropagateChangeUp = true;
        for (Parameter p : uniqueParameters) {
            p.fireParameterChangedEvent();
        }
        doNotPropagateChangeUp = false;
        fireParameterChangedEvent(-1, ChangeType.ALL_VALUES_CHANGED);
    }
}
