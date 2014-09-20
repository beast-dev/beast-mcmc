package dr.app.beauti.components.linkedparameters;

import dr.app.beauti.options.Operator;
import dr.app.beauti.options.Parameter;
import dr.app.beauti.types.PriorType;

import java.util.List;

/**
 * @author Andrew Rambaut
 * @author Marc A. Suchard
 */
public class LinkedParameter {

    private String name;
    final private Parameter argumentParameter;
    final private Operator argumentOperator;
    final private LinkedParameterComponentOptions options;

    /**
     * A simple data class to store the definition of a linked parameter.
     * @param name the name of the parameter
     * @param argumentParameter the control parameter
     * @param argumentOperator the control parameter's operator
     */
    public LinkedParameter(String name, Parameter argumentParameter, Operator argumentOperator, LinkedParameterComponentOptions options) {
        this.name = name;
        this.argumentParameter = argumentParameter;
        this.argumentOperator = argumentOperator;
        this.options = options;
    }

    public Parameter getArgumentParameter() {
        return argumentParameter;
    }

    public Operator getArgumentOperator() {
        return argumentOperator;
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
        argumentParameter.setName(name);
        argumentOperator.setName(name);
        this.name = name;

    }

}
