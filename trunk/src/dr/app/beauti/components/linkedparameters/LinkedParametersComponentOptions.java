package dr.app.beauti.components.linkedparameters;

import dr.app.beauti.components.hpm.HierarchicalModelComponentGenerator;
import dr.app.beauti.components.hpm.HierarchicalPhylogeneticModel;
import dr.app.beauti.options.*;
import dr.app.beauti.priorsPanel.PriorsPanel;
import dr.app.beauti.types.PriorScaleType;
import dr.app.beauti.types.PriorType;

import javax.swing.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * @author Andrew Rambaut
 * @author Marc A. Suchard
 * @version $Id$
 */
public class LinkedParametersComponentOptions implements ComponentOptions {

    public LinkedParametersComponentOptions(final BeautiOptions options) {
        this.options = options;
        linkedParametersMap = new HashMap<Parameter, LinkedParameters>();
    }

    public void createParameters(final ModelOptions modelOptions) {
        // Do nothing; this is only called at launch
    }

    public void selectParameters(final ModelOptions modelOptions, final List<Parameter> params) {
        for (LinkedParameters linkedParameters : getLinkedParametersList()) {
            params.addAll(linkedParameters.getArgumentParameterList());
        }
    }

    public void selectStatistics(final ModelOptions modelOptions, final List<Parameter> stats) {
        // No statistics
    }

    public void selectOperators(final ModelOptions modelOptions, final List<Operator> ops) {
        for (LinkedParameters linkedParameters : getLinkedParametersList()) {
            for (Parameter parameter : linkedParameters.getArgumentParameterList()) {
            }
        }
    }

    public boolean nameExits(String name) {
        boolean found = false;

        if (options.parameterExists(name)) {
            return true;
        }

        for (LinkedParameters linkedParameters : getLinkedParametersList()) {
            if (linkedParameters.getName().equalsIgnoreCase(name)) {
                found = true;
                break;
            }
        }
        return found;
    }

    public LinkedParameters addLinkedParameters(String name, List<Parameter> parameterList) {
        List<Parameter> argumentList = new ArrayList<Parameter>();

        Parameter sourceParameter = parameterList.get(0);

        //argumentList.add();

        LinkedParameters linkedParameters = new LinkedParameters(name, argumentList);
        for (Parameter parameter : parameterList) {
            linkedParametersMap.put(parameter, linkedParameters);
        }

        return linkedParameters;
    }

    public void removeParameters(PriorsPanel priorsPanel, List<Parameter> parametersToRemove, boolean caution) {
        for (Parameter parameter : parametersToRemove) {
            linkedParametersMap.remove(parameter);
        }
    }

    public List<LinkedParameters> getLinkedParametersList() {
        return new ArrayList<LinkedParameters>(linkedParametersMap.values());
    }

    public LinkedParameters getLinkedParametersForParameter(Parameter parameter) {
        return linkedParametersMap.get(parameter);
    }

    public List<Parameter> getLinkedParameters(LinkedParameters linkedParameters) {
        List<Parameter> parameters = new ArrayList<Parameter>();
        for (Parameter parameter : linkedParametersMap.keySet()) {
            if (linkedParametersMap.get(parameter) == linkedParameters) {
                parameters.add(parameter);
            }
        }
        return parameters;
    }



    final private BeautiOptions options;
    final private Map<Parameter, LinkedParameters> linkedParametersMap;


//    public void generatePriors(final XMLWriter writer) {
//        for (HierarchicalPhylogeneticModel hpm : hpmList) {
//            hpm.generatePriors(writer);
//        }
//    }

    public boolean isEmpty() {
        return linkedParametersMap.isEmpty();
    }
}