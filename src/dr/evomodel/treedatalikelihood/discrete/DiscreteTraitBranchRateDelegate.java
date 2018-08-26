/*
 * DiscreteTraitBranchRateDelegate.java
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

package dr.evomodel.treedatalikelihood.discrete;

import dr.evolution.tree.Tree;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.evomodel.treedatalikelihood.BeagleDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.preorder.AbstractDiscreteTraitDelegate;

/**
 * @author Xiang Ji
 * @author Marc A. Suchard
 */
public class DiscreteTraitBranchRateDelegate extends AbstractDiscreteTraitDelegate {
    public DiscreteTraitBranchRateDelegate(String name, Tree tree, BeagleDataLikelihoodDelegate likelihoodDelegate) {
        super(name, tree, likelihoodDelegate);
    }

    protected void cacheDifferentialMassMatrix(boolean cacheSquaredMatrix) {
        for (int i = 0; i < evolutionaryProcessDelegate.getSubstitutionModelCount(); i++) {
            double[] infinitesimalMatrix = new double[stateCount * stateCount];
            SubstitutionModel substitutionModel = evolutionaryProcessDelegate.getSubstitutionModel(i);
            substitutionModel.getInfinitesimalMatrix(infinitesimalMatrix);
            evolutionaryProcessDelegate.cacheFirstOrderDifferentialMatrix(beagle, i, infinitesimalMatrix);
            if (cacheSquaredMatrix) {
                double[] infinitesimalMatrixSquared = new double[stateCount * stateCount];
                for (int l = 0; l < stateCount; l++) {
                    for (int j = 0; j < stateCount; j++) {
                        double sumOverState = 0.0;
                        for (int k = 0; k < stateCount; k++) {
                            sumOverState += infinitesimalMatrix[l * stateCount + k] * infinitesimalMatrix[k * stateCount + j];
                        }
                        infinitesimalMatrixSquared[l * stateCount + j] = sumOverState;
                    }
                }
                evolutionaryProcessDelegate.cacheSecondOrderDifferentialMatrix(beagle, i, infinitesimalMatrixSquared);
            }
        }
    }

}
