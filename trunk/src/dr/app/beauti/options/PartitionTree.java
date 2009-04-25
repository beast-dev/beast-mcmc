package dr.app.beauti.options;

import java.util.List;
import java.util.ArrayList;

/**
 * @author Andrew Rambaut
 * @version $Id$
 */
public class PartitionTree extends ModelOptions {

    public PartitionTree(BeautiOptions options, DataPartition partition) {
        this(options, partition.getName());
    }

    /**
     * A copy constructor
     *
     * @param options the beauti options
     * @param name    the name of the new model
     * @param source  the source model
     */
    public PartitionTree(BeautiOptions options, String name, PartitionTree source) {
        this(options, name);
    }

    public PartitionTree(BeautiOptions options, String name) {
        this.options = options;
        this.name = name;
  }

    public String getName() {
        return name;
    }

    public void setName(String name) {
        this.name = name;
    }

    public List<Operator> getOperators() {
        List<Operator> operators = new ArrayList<Operator>();

        return operators;
    }

    /**
     * @param includeRelativeRates true if relative rate parameters should be added
     * @return a list of parameters that are required
     */
    List<Parameter> getParameters(boolean includeRelativeRates) {

        List<Parameter> params = new ArrayList<Parameter>();


        return params;
    }

    public Parameter getParameter(String name) {

        if (name.startsWith(getName())) {
            name = name.substring(getName().length() + 1);
        }
        Parameter parameter = parameters.get(name);

        if (parameter == null) {
            throw new IllegalArgumentException("Parameter with name, " + name + ", is unknown");
        }

        parameter.setPrefix(getPrefix());

        return parameter;
    }

    public Operator getOperator(String name) {

        Operator operator = operators.get(name);

        if (operator == null) throw new IllegalArgumentException("Operator with name, " + name + ", is unknown");

        operator.setPrefix(getName());

        return operator;
    }

    public String toString() {
        return getName();
    }


    public String getPrefix() {
        String prefix = "";
        if (options.getActivePartitionTrees().size() > 1) {
            // There is more than one active partition model
            prefix += getName() + ".";
        }
        return prefix;
    }

// Instance variables

    private final BeautiOptions options;
    public String name;

}
