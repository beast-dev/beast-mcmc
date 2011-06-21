package dr.app.beauti.components.hpm;

import dr.app.beauti.options.ModelOptions;
import dr.app.beauti.options.Operator;
import dr.app.beauti.options.Parameter;
import dr.app.beauti.types.OperatorType;
import dr.app.beauti.types.PriorType;
import dr.app.beauti.util.XMLWriter;

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

    public void selectOperators(ModelOptions modelOptions, final List<Operator> ops) {
        // Do nothing because Gibbs operator format do not fit into the current Operator implementation
    }

//    public void removeParameter(Parameter parameter) {
//        if (argumentParameterList.contains(parameter)) {
//            argumentParameterList.remove(parameter);
//        }
//    }
//
//    public void generateDistribution(final XMLWriter writer) {
//        System.err.println("Generating distribution for " + getName());
//    }
//
//    public void generatePriors(final XMLWriter writer) {
//        System.err.println("Generating priors for " + getName());
//    }
}
