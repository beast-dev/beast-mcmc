package dr.app.beauti.components.hpm;

import dr.app.beauti.options.*;
import dr.app.beauti.priorsPanel.PriorsPanel;

import javax.swing.*;
import java.util.ArrayList;
import java.util.List;

/**
 * @author Marc A. Suchard
 * @version $Id$
 */
public class HierarchicalModelComponentOptions implements ComponentOptions {

    public HierarchicalModelComponentOptions(final BeautiOptions options) {
        this.options = options;
        hpmList = new ArrayList<HierarchicalPhylogeneticModel>();
    }

    public void createParameters(final ModelOptions modelOptions) {
        // Do nothing; this is only called at launch
    }

    public void selectParameters(final ModelOptions modelOptions, final List<Parameter> params) {
        for (HierarchicalPhylogeneticModel hpm : hpmList) {
            if (!hpm.isEmpty()) {
                List<Parameter> hpmParameterList = hpm.getConditionalParameterList();
                params.addAll(hpmParameterList);
            }
        }        
    }

    public void selectStatistics(final ModelOptions modelOptions, final List<Parameter> stats) {
        // No statistics
    }

    public void selectOperators(final ModelOptions modelOptions, final List<Operator> ops) {
        for (HierarchicalPhylogeneticModel hpm : hpmList) {
            if (!hpm.isEmpty()) {
                hpm.selectOperators(ops);
            }
        }
    }

    public boolean modelExists(String name) {
        boolean found = false;
        for (HierarchicalPhylogeneticModel hpm : hpmList) {
            if (hpm.getName().compareTo(name) == 0) {
                found = true;
                break;
            }
        }
        return found;        
    }

    public void addHPM(String text, List<Parameter> parameterList) {
        HierarchicalPhylogeneticModel hpm = new HierarchicalPhylogeneticModel(text, parameterList, null, null);
        hpmList.add(hpm);
    }

    public boolean isHierarchicalParameter(Parameter parameter) {
        boolean found = false;
        for (HierarchicalPhylogeneticModel hpm : hpmList) {
            if (hpm.getArgumentParameterList().contains(parameter)) {
                found = true;
                break;
            }
        }
        return found;
    }

    public int removeParameter(PriorsPanel priorsPanel, Parameter parameter) {
        HierarchicalPhylogeneticModel toRemove = null;
        for (HierarchicalPhylogeneticModel hpm : hpmList) {
            List<Parameter> parameterList = hpm.getArgumentParameterList();
            if (parameterList.contains(parameter)) {
                if (parameterList.size() == 2) {
                    String modelName = hpm.getName();
                    // Throw special warning
                    int option = JOptionPane.showConfirmDialog(priorsPanel,
                        "Removing this parameter from HPM '" + modelName + "' will result in one only\n" +
                        "parameter remaining the HPM.  Single parameter models are not recommended.\n" +
                        "Continue?",
                        "HPM warning",
                        JOptionPane.YES_NO_OPTION,
                        JOptionPane.WARNING_MESSAGE);
                    if (option == JOptionPane.NO_OPTION) {
                        return JOptionPane.NO_OPTION;
                    }
                }
                parameterList.remove(parameter);
            }           
            if (hpm.isEmpty())  {
                toRemove = hpm;
            }
        }
        if (toRemove != null) {
            hpmList.remove(toRemove);
        }
        return JOptionPane.YES_OPTION;
    }
        
    final private BeautiOptions options;
    final private List<HierarchicalPhylogeneticModel> hpmList;


}