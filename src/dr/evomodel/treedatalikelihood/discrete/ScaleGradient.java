/*
 * ScaleGradient.java
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

package dr.evomodel.treedatalikelihood.discrete;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.treedatalikelihood.BeagleDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.inference.model.Parameter;

/**
 * @author Marc A. Suchard
 */
public class ScaleGradient extends HyperParameterBranchRateGradient {

    //declaring String constants for use in BEAUti
    public static final String SCALE_GRADIENT = "scaleGradient";

    public ScaleGradient(String traitName, TreeDataLikelihood treeDataLikelihood,
                         BeagleDataLikelihoodDelegate likelihoodDelegate,
                         Parameter locationScaleParameter, boolean useHessian) {

        super(traitName, treeDataLikelihood, likelihoodDelegate, locationScaleParameter, useHessian);
    }

    @Override
    double[] getDifferential(Tree tree, NodeRef node) {
        double rate = branchRateModel.getBranchRate(tree, node);

        // TODO Can out out of this class (I think), if we provide both transform() and inverse() here.

        double tmp = (Math.log(rate / locationScaleTransform.getLocation(tree, node)) - locationScaleTransform.getTransformMu())
                /(locationScaleTransform.getTransformSigma() * locationScaleTransform.getTransformSigma()) - 1.0;

        return new double[] {tmp * rate * locationScaleTransform.getScale(tree, node) / (1.0 + locationScaleTransform.getScale(tree, node) * locationScaleTransform.getScale(tree, node))};
    }
}
