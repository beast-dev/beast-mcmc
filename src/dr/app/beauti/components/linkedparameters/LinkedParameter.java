package dr.app.beauti.components.linkedparameters;

import dr.app.beauti.options.Parameter;
import dr.app.beauti.types.PriorType;

import java.util.List;

/**
 * @author Andrew Rambaut
 * @author Marc A. Suchard
 */
public class LinkedParameter {

    final private String name;
    final private List<Parameter> argumentParameterList;

    /**
     * A simple data class to store the definition of a linked parameter.
     * @param name the name of the parameter
     * @param argumentParameterList the list of 'control' parameters
     */
    public LinkedParameter(String name, List<Parameter> argumentParameterList) {
        this.name = name;
        this.argumentParameterList = argumentParameterList;
    }

    public List<Parameter> getArgumentParameterList() {
        return argumentParameterList;
    }

    public String getName() {
        return name;
    }

}
