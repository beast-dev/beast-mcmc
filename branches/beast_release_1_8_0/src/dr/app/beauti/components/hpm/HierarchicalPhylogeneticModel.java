package dr.app.beauti.components.hpm;

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

    public static final String TIP_TOOL = "<html>Select two or more parameters. " +
            "HPMs are described in<br>Suchard et al. (2003) Syst Biol 52, 649 - 664 " +
            "and <br>Edo-Matas et al. (2011) Mol Biol Evol 28, 1605 - 1616.</html>";

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

//    public void selectOperators(ModelOptions modelOptions, final List<Operator> ops) {
//        // Do nothing because Gibbs operator format do not fit into the current Operator implementation
//    }
//
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
