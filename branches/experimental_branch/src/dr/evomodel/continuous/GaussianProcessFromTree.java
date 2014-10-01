/*
 * GaussianProcessFromTree.java
 *
 * Copyright (c) 2002-2014 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodel.continuous;

import dr.evolution.tree.NodeRef;
import dr.math.MathUtils;
import dr.math.distributions.GaussianProcessRandomGenerator;

/**
 * @author Marc A. Suchard
 * @author Max R. Tolkoff
 */
public class GaussianProcessFromTree implements GaussianProcessRandomGenerator {

    private final FullyConjugateMultivariateTraitLikelihood traitModel;

    public GaussianProcessFromTree(FullyConjugateMultivariateTraitLikelihood traitModel) {
        this.traitModel = traitModel;
    }

    //    boolean firstTime=true;
    public double[] nextRandomFast(int i) {

        double[] random = new double[traitModel.getTreeModel().getExternalNodeCount()];
        NodeRef root = traitModel.getTreeModel().getRoot();
//        if(traitModel.getTreeModel().isExternal(root)) {
//            random[0] = traitModel.getTreeModel().getMultivariateNodeTrait(root, traitModel.getTraitName())[i];
//        }
//        else{


        double[][] var = MultivariateTraitUtils.computeTreeVariance(traitModel, true);
//        if(firstTime) {
//            for (int j = 0; j < var[0].length; j++) {
//                for (int k = 0; k < var[0].length; k++) {
//                    if(j!=k)
//                        var[j][k] = var[j][k] / Math.sqrt(var[k][k] * var[j][j]);
//                }
//            }
//
//
//
//            for (int j = 0; j < var[0].length; j++) {
//                String empty = "";
//                for (int k = 0; k < var[0].length; k++) {
//                    empty += Double.toString(var[j][k]) + "\t";
//                }
//                System.out.println(empty);
//            }
//            firstTime=false;
//        }
        nextRandomFast(0, root, random); // TODO 0 is not right;  draw root from its prior
//        }
        return random;
    }

    private void nextRandomFast(double currentValue, NodeRef currentNode, double[] random) {
        double rescaledLength;
        rescaledLength = traitModel.getRescaledBranchLengthForPrecision(currentNode);
        if (traitModel.getTreeModel().isExternal(currentNode)) {
            random[currentNode.getNumber()] = currentValue + MathUtils.nextGaussian() * Math.sqrt(rescaledLength);
        } else {
            int childCount = traitModel.getTreeModel().getChildCount(currentNode);
            double newValue = currentValue + Math.sqrt(rescaledLength) * MathUtils.nextGaussian();
            for (int i = 0; i < childCount; i++) {
                nextRandomFast(newValue, traitModel.getTreeModel().getChild(currentNode, i), random);
            }
        }
    }

    @Override
    public Object nextRandom() {
        return nextRandomFast(0);
    }

    @Override
    public double logPdf(Object x) {
        return 0;
    }
}
