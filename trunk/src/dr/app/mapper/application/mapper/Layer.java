package dr.app.mapper.application.mapper;

import dr.inference.trace.Trace;

import java.util.List;

/**
 * @author Andrew Rambaut
 * @version $Id$
 */
public class Layer {
    public Layer(String name, Trace trace, LayerType type, List<String> labels) {
        this.name = name;
        this.trace = trace;
        this.type = type;
        this.labels = labels;
    }

    public String getName() {
        return name;
    }

    public Trace getTrace() {
        return trace;
    }

    public LayerType getType() {
        return type;
    }

    public List<String> getLabels() {
        return labels;
    }

    private String name;
    private Trace trace;
    private LayerType type;

    private List<String> labels;
}
