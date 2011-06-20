package dr.app.beauti.components.hpm;

import dr.app.beauti.options.Operator;
import dr.app.beauti.options.Parameter;
import dr.app.beauti.types.PriorType;

import java.util.List;

/**
 * @author Marc A. Suchard
 */
public class HierarchicalPhylogeneticModel {

    final private String name;
    final private List<Parameter> argumentParameterList;
    final private List<Parameter> conditionalParameterList;
    final private PriorType priorType;

    public HierarchicalPhylogeneticModel(String name, List<Parameter> argumentParameterList,
                                         List<Parameter> conditionalParameterList, PriorType priorType) {
        this.name = name;
        this.argumentParameterList = argumentParameterList;
        this.conditionalParameterList = conditionalParameterList;
        this.priorType = priorType;
    }

    public List<Parameter> getArgumentParameterList() {
        return argumentParameterList;
    }

    public List<Parameter> getConditionalParameterList() {
        return conditionalParameterList;
    }

    public String getName() {
        return name;
    }

    public PriorType getPriorType() {
        return priorType;
    }

    public boolean isEmpty() {
        return getArgumentParameterList().size() == 0;
    }

    public void selectOperators(final List<Operator> ops) {
        // TODO
    }

    public void removeParameter(Parameter parameter) {
        if (argumentParameterList.contains(parameter)) {
            argumentParameterList.remove(parameter);
        }
    } 
}
