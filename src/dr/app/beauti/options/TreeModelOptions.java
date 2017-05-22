/*
 * TreeModelOptions.java
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

package dr.app.beauti.options;


import dr.app.beauti.types.PriorType;
import dr.evolution.util.Taxa;

import java.util.List;


/**
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @author Walter Xie
 * @version $Id$
 */
public class TreeModelOptions extends ModelOptions {
    private static final long serialVersionUID = 5328826852511460749L;

    // Instance variables
    private final BeautiOptions options;


    public TreeModelOptions(BeautiOptions options) {
        this.options = options;

        initModelParametersAndOpererators();
    }

    @Override
    public void initModelParametersAndOpererators() {

    }

    @Override
    public List<Parameter> selectParameters(List<Parameter> params) {
        return params;
    }

    @Override
    public List<Operator> selectOperators(List<Operator> ops) {
        return ops;
    }

    @Override
    public String getPrefix() {
        return null;
    }

    public double getExpectedAvgBranchLength(double rootHeight) {
        double sum = 0;
        int taxonCount = options.taxonList.getTaxonCount();

        for (int i = 2; i <= taxonCount; i++) {
            sum += (double) 1 / i;
        }

        return rootHeight * sum / (double) (2 * taxonCount - 2);
    }

    public int isNodeCalibrated(PartitionTreeModel treeModel) {
        if (treeModel == null) {
            return -1;
        }
        if (isNodeCalibrated(treeModel.getParameter("treeModel.rootHeight"))) {
            return 0; // root node
        } else if (options.getKeysFromValue(options.taxonSetsTreeModel, treeModel).size() > 0) {
            Taxa taxonSet = (Taxa) options.getKeysFromValue(options.taxonSetsTreeModel, treeModel).get(0);
            Parameter tmrca = options.statistics.get(taxonSet);
            if (tmrca != null && isNodeCalibrated(tmrca)) {
                return 1; // internal node (tmrca) with a proper prior
            }
            return -1;
        } else {
            return -1;
        }
    }

    public boolean isNodeCalibrated(Parameter parameter) {
        return (parameter.taxaId != null && hasProperPriorOn(parameter)) // param.taxa != null is TMRCA
                || (parameter.getBaseName().endsWith("treeModel.rootHeight") && hasProperPriorOn(parameter));
    }

    private boolean hasProperPriorOn(Parameter para) {
        return para.priorType == PriorType.EXPONENTIAL_PRIOR
                || (para.priorType == PriorType.UNIFORM_PRIOR && para.uniformLower > 0 && para.uniformUpper < Double.POSITIVE_INFINITY)
                || para.priorType == PriorType.LAPLACE_PRIOR
                || para.priorType == PriorType.NORMAL_PRIOR
                || para.priorType == PriorType.LOGNORMAL_PRIOR
                || para.priorType == PriorType.GAMMA_PRIOR
                || para.priorType == PriorType.INVERSE_GAMMA_PRIOR
                || para.priorType == PriorType.BETA_PRIOR
                || para.priorType == PriorType.CTMC_RATE_REFERENCE_PRIOR
                || para.priorType == PriorType.LOGNORMAL_HPM_PRIOR
                || para.priorType == PriorType.POISSON_PRIOR;
    }

}
