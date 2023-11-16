/*
 * DiscreteTraitBranchRateGradient.java
 *
 * Copyright (c) 2002-2022 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodel.treedatalikelihood.discrete;

import dr.evomodel.substmodel.GlmSubstitutionModel;
import dr.evomodel.treedatalikelihood.BeagleDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.inference.loggers.LogColumn;
import dr.inference.model.Parameter;
import dr.util.Citation;

import java.util.List;

/**
 * @author Filippo Monti
 * @author Marc A. Suchard
 */

public class GamGpSubstitutionModelGradient extends AbstractLogAdditiveSubstitutionModelGradient {

    public GamGpSubstitutionModelGradient(String traitName,
                                          TreeDataLikelihood treeDataLikelihood,
                                          BeagleDataLikelihoodDelegate likelihoodDelegate,
                                          GlmSubstitutionModel substitutionModel) {
        super(traitName, treeDataLikelihood, likelihoodDelegate, substitutionModel,
                ApproximationMode.FIRST_ORDER);
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    protected double preProcessNormalization(double[] differentials, double[] generator, boolean normalize) {
        return 0;
    }

    @Override
    double processSingleGradientDimension(int dim, double[] differentials, double[] generator, double[] pi,
                                          boolean normalize, double normalizationConstant) {
        return 0;
    }

    @Override
    public Parameter getParameter() {
        return null;
    }

    @Override
    public LogColumn[] getColumns() {
        return new LogColumn[0];
    }

    @Override
    public Citation.Category getCategory() {
        return null;
    }

    @Override
    public String getDescription() {
        return null;
    }

    @Override
    public List<Citation> getCitations() {
        return null;
    }
}
