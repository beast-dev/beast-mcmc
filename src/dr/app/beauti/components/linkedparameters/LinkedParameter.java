package dr.app.beauti.components.linkedparameters;

import dr.app.beauti.options.Parameter;
import dr.app.beauti.types.PriorType;

import java.util.List;

/**
 * @author Andrew Rambaut
 * @author Marc A. Suchard
 */
public class LinkedParameter {

    private String name;
    final private List<Parameter> argumentParameterList;
    final private LinkedParameterComponentOptions options;

    /**
     * A simple data class to store the definition of a linked parameter.
     * @param name the name of the parameter
     * @param argumentParameterList the list of 'control' parameters
     */
    public LinkedParameter(String name, List<Parameter> argumentParameterList, LinkedParameterComponentOptions options) {
        this.name = name;
        this.argumentParameterList = argumentParameterList;
        this.options = options;
    }

    public List<Parameter> getArgumentParameterList() {
        return argumentParameterList;
    }

    public List<Parameter> getDependentParameterList() {
        return options.getDependentParameters(this);
    }

    public void linkDependentParameters(List<Parameter> parameterList) {
        options.setDependentParameters(this, parameterList);
    }

    public String getName() {
        return name;
    }

    public void setName(String name) {
        this.name = name;
    }

}
