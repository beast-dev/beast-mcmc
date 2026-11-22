/*
 * HierarchicalModelComponentOptions.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
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
 *
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
 */
public class HierarchicalModelComponentOptions implements ComponentOptions {


    // explicit, so that files saved by this version still load after an unrelated change
    private static final long serialVersionUID = 8813437700244230279L;
    public HierarchicalModelComponentOptions(final BeautiOptions options) {
        this.options = options;
        hpmList = new ArrayList<HierarchicalPhylogeneticModel>();
    }

    public void createParameters(final ModelOptions modelOptions) {
        // Do nothing; this is only called at launch
    }

    public void selectParameters(final ModelOptions modelOptions, final List<Parameter> params) {
        // Do nothing because priors are previously defined for idref in operator schedule
    }

    public void selectStatistics(final ModelOptions modelOptions, final List<Parameter> stats) {
        // No statistics
    }

    public void selectOperators(final ModelOptions modelOptions, final List<Operator> ops) {
        // Do nothing because Gibbs operator format do not fit into the current Operator implementation
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

    /**
     * Drops from every hierarchical model the parameters that are no longer part of the analysis, and
     * discards any model left with fewer than two of them.
     *
     * A model holds its parameters directly, and changing the analysis - partitioning into codon positions,
     * say, or choosing a different substitution model - does not destroy the parameters it replaces, it only
     * stops them being used. The models were written out whatever had become of their parameters, so the
     * file could refer to parameters that it did not declare and BEAST would not read it.
     *
     * A parameter taken out of a model is given back an ordinary prior, because the density of a parameter
     * in a hierarchical model is written as part of that model and nothing would be written for it otherwise.
     *
     * @param activeParameters the parameters of the analysis as it now stands
     * @return what was dropped, for reporting to the user, or null if the models were already consistent
     */
    public String pruneToActiveParameters(List<Parameter> activeParameters) {

        List<String> messages = new ArrayList<String>();
        List<HierarchicalPhylogeneticModel> toRemove = new ArrayList<HierarchicalPhylogeneticModel>();
        List<Parameter> needingAPrior = new ArrayList<Parameter>();

        for (HierarchicalPhylogeneticModel hpm : hpmList) {
            List<Parameter> parameterList = hpm.getArgumentParameterList();

            List<Parameter> dropped = new ArrayList<Parameter>();
            for (Parameter parameter : parameterList) {
                if (!activeParameters.contains(parameter)) {
                    dropped.add(parameter);
                }
            }

            if (dropped.isEmpty()) {
                continue;
            }

            parameterList.removeAll(dropped);
            for (Parameter parameter : dropped) {
                detach(parameter);
            }

            if (parameterList.size() < 2) {
                // a model over a single parameter is not one, so the whole model goes
                for (Parameter parameter : new ArrayList<Parameter>(parameterList)) {
                    detach(parameter);
                    if (activeParameters.contains(parameter)) {
                        // this one is still in the analysis, so the user will want to give it a prior
                        needingAPrior.add(parameter);
                    }
                }
                parameterList.clear();
                toRemove.add(hpm);
                messages.add("the hierarchical model '" + hpm.getName() + "' has been removed");
            } else {
                messages.add(dropped.size() + (dropped.size() == 1 ? " parameter has" : " parameters have") +
                        " been removed from the hierarchical model '" + hpm.getName() + "'");
            }
        }

        hpmList.removeAll(toRemove);

        if (messages.isEmpty()) {
            return null;
        }

        StringBuilder message = new StringBuilder();
        for (String line : messages) {
            if (message.length() > 0) {
                message.append(", and ");
            }
            message.append(line);
        }
        message.append(", because ");
        message.append(messages.size() == 1 ? "a parameter it was over is" : "parameters they were over are");
        message.append(" no longer part of the analysis.");

        if (!needingAPrior.isEmpty()) {
            message.append("\n\n");
            message.append(needingAPrior.size() == 1 ? "The parameter " : "The parameters ");
            for (int i = 0; i < needingAPrior.size(); i++) {
                if (i > 0) {
                    message.append(i == needingAPrior.size() - 1 ? " and " : ", ");
                }
                message.append(needingAPrior.get(i).getName());
            }
            message.append(needingAPrior.size() == 1 ? " is still" : " are still");
            message.append(" in the analysis and now ");
            message.append(needingAPrior.size() == 1 ? "has a default prior" : "have default priors");
            message.append(" instead.");
        }

        return message.toString();
    }

    /**
     * Takes a parameter out of a hierarchical model, giving it back an ordinary prior: a uniform one where
     * the parameter is bounded, and the improper uniform otherwise. The prior is marked as not having been
     * set by the user, so that the Priors panel draws attention to it.
     */
    private void detach(Parameter parameter) {
        parameter.isLinked = false;
        parameter.linkedName = null;

        if (parameter.priorType == PriorType.NORMAL_HPM_PRIOR
                || parameter.priorType == PriorType.LOGNORMAL_HPM_PRIOR
                || parameter.priorType == PriorType.EXPONENTIAL_HPM_PRIOR) {
            boolean bounded = !Double.isInfinite(parameter.getLowerBound())
                    && !Double.isInfinite(parameter.getUpperBound());
            parameter.priorType = bounded ? PriorType.UNIFORM_PRIOR : PriorType.NONE_IMPROPER;
            parameter.setPriorEdited(false);
        }
    }


    final private BeautiOptions options;
    final private List<HierarchicalPhylogeneticModel> hpmList;



    public boolean isEmpty() {
        return hpmList.isEmpty();
    }
}