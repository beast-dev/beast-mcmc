/*
 * FullyConjugateTreeTipsPotentialDerivative.java
 *
 * Copyright (c) 2002-2017 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodel.continuous.hmc;

import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.continuous.ContinuousDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.preorder.NewTipFullConditionalDistributionDelegate;
import dr.evomodel.treedatalikelihood.preorder.NormalSufficientStatistics;
import dr.inference.hmc.PrecisionVectorProductProvider;
import dr.inference.model.Parameter;
import dr.xml.Reportable;

import java.util.List;

/**
 * @author Max Tolkoff
 * @author Marc A. Suchard
 */
public class TreePrecisionTraitProductProvider implements PrecisionVectorProductProvider, Reportable {

    private final TreeTrait<List<NormalSufficientStatistics>> fullConditionalDensity;
    private final Tree tree;
    private final int nTaxa;

    public TreePrecisionTraitProductProvider(TreeDataLikelihood treeDataLikelihood,
                                             ContinuousDataLikelihoodDelegate likelihoodDelegate,
                                             String traitName) {

        String fcdName = NewTipFullConditionalDistributionDelegate.getName(traitName);
        if (treeDataLikelihood.getTreeTrait(fcdName) == null) {
            likelihoodDelegate.addNewFullConditionalDensityTrait(traitName);
        }

        this.fullConditionalDensity = treeDataLikelihood.getTreeTrait(fcdName);
        this.tree = treeDataLikelihood.getTree();
        this.nTaxa = tree.getExternalNodeCount();
    }

    @Override
    public double[] getMultiplicationResultant(Parameter vector) {

        // TODO ensure that vector == underlying data parameter for likelihoodDelegate

        double[] result = new double[vector.getDimension()];

        List<NormalSufficientStatistics> statistics = fullConditionalDensity.getTrait(tree, null);

        // TODO Fill in values

        return result;
    }

    @Override
    public String getReport() {
        return (new dr.math.matrixAlgebra.Vector(getMultiplicationResultant(null))).toString();
    }
}
