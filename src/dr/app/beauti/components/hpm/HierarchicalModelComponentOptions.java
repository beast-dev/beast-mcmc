/*
 * HierarchicalModelComponentOptions.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
 *
 * This file is part of BEAST.
 * See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * BEAST is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 *  BEAST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package dr.app.beauti.components.hpm;

import dr.app.beauti.options.*;
import dr.app.beauti.priorspanel.PriorsPanel;
import dr.app.beauti.types.PriorScaleType;
import dr.app.beauti.types.PriorType;

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
        // Do nothing because priors are previously defined for idref in operator schedule
//        for (HierarchicalPhylogeneticModel hpm : hpmList) {
//            if (!hpm.isEmpty()) {
//                List<Parameter> hpmParameterList = hpm.getConditionalParameterList();
//                params.addAll(hpmParameterList);
//            }
//        }
    }

    public void selectStatistics(final ModelOptions modelOptions, final List<Parameter> stats) {
        // No statistics
    }

    public void selectOperators(final ModelOptions modelOptions, final List<Operator> ops) {
        // Do nothing because Gibbs operator format do not fit into the current Operator implementation
//        for (HierarchicalPhylogeneticModel hpm : hpmList) {
//            if (!hpm.isEmpty()) {
//                hpm.selectOperators(modelOptions, ops);
//            }
//        }
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

    public HierarchicalPhylogeneticModel addHPM(String text, List<Parameter> parameterList, PriorType priorType) {
        List<Parameter> argumentList = new ArrayList<Parameter>();

        // TODO May have to remove these constructors
        String meanName = text + HierarchicalModelComponentGenerator.MEAN_SUFFIX;
        Parameter mean = options.parameterExists(meanName) ?
                options.getParameter(meanName) :
                options.createParameterNormalPrior(meanName, "Unknown mean of HPM",
                PriorScaleType.NONE, 0.0, 0.0, 1.0, 0.0);
        argumentList.add(mean);

        String precisionName = text + HierarchicalModelComponentGenerator.PRECISION_SUFFIX;
        Parameter precision = options.parameterExists(precisionName) ?
                options.getParameter(precisionName) :
                options.createParameterGammaPrior(precisionName, "Unknown precision of HPM",
                PriorScaleType.NONE, 1.0, 0.001, 1000.0, true);
        argumentList.add(precision);

        HierarchicalPhylogeneticModel hpm = new HierarchicalPhylogeneticModel(text, parameterList, argumentList, priorType);
        hpmList.add(hpm);

        for (Parameter parameter : parameterList) {
            parameter.linkedName = hpm.getName();
        }
        return hpm;
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

    public int removeParameter(PriorsPanel priorsPanel, Parameter parameter, boolean caution) {
        HierarchicalPhylogeneticModel toRemove = null;
        for (HierarchicalPhylogeneticModel hpm : hpmList) {
            List<Parameter> parameterList = hpm.getArgumentParameterList();
            if (parameterList.contains(parameter)) {
                if (caution && parameterList.size() == 2 && priorsPanel != null) {
                    String modelName = hpm.getName();
                    // Throw special warning
                    int option = JOptionPane.showConfirmDialog(priorsPanel,
                        "Removing this parameter from HPM '" + modelName + "' will result in only one\n" +
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

    public List<HierarchicalPhylogeneticModel> getHPMList() {
        return hpmList;
    }

//    public void generateDistributions(final XMLWriter writer) {
//        for (HierarchicalPhylogeneticModel hpm : hpmList) {
//            hpm.generateDistribution(writer);
//        }
//    }

    final private BeautiOptions options;
    final private List<HierarchicalPhylogeneticModel> hpmList;


//    public void generatePriors(final XMLWriter writer) {
//        for (HierarchicalPhylogeneticModel hpm : hpmList) {
//            hpm.generatePriors(writer);
//        }
//    }

    public boolean isEmpty() {
        return hpmList.isEmpty();
    }
}