/*
 * ModelExtensionProvider.java
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

package dr.evomodel.treedatalikelihood.preorder;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.evomodel.treedatalikelihood.continuous.ContinuousDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.continuous.ContinuousTraitPartialsProvider;
import dr.inference.model.MatrixParameterInterface;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;

import static dr.math.matrixAlgebra.missingData.MissingOps.safeInvert2;

public interface ModelExtensionProvider extends ContinuousTraitPartialsProvider {

    ContinuousExtensionDelegate getExtensionDelegate(ContinuousDataLikelihoodDelegate delegate,
                                                     TreeTrait treeTrait,
                                                     Tree tree);


    interface NormalExtensionProvider extends ModelExtensionProvider {

        boolean diagonalVariance();

        DenseMatrix64F getExtensionVariance();

        DenseMatrix64F getExtensionVariance(NodeRef node);

        MatrixParameterInterface getExtensionPrecision();

        void chainRuleWrtVariance(double[] gradient, NodeRef node);

        @Override
        default void updateTipDataGradient(DenseMatrix64F precision, DenseMatrix64F variance,
                                           NodeRef node, int offset, int dimGradient) {

            extendTipDataGradient(this, precision, variance, node, offset, dimGradient);
        }

        static void extendTipDataGradient(NormalExtensionProvider provider,
                                          DenseMatrix64F precision, DenseMatrix64F variance,
                                          NodeRef node, int offset, int dimGradient) {

            if (offset != 0 || dimGradient != provider.getTraitDimension()) {
                throw new RuntimeException("not implemented for subset of model.");
            }

            DenseMatrix64F samplingVariance = provider.getExtensionVariance(node);
            CommonOps.addEquals(samplingVariance, variance);
            safeInvert2(samplingVariance, precision, false);
        }

        @Override
        default boolean needToUpdateTipDataGradient(int offset, int dimGradient) {
            return true;
        }

    }
}


