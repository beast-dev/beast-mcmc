/*
 * AbstractValuesViaFullConditionalDelegate.java
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
import dr.evomodel.continuous.MultivariateDiffusionModel;
import dr.evomodel.treedatalikelihood.continuous.*;
import dr.math.matrixAlgebra.WrappedVector;
import org.ejml.data.DenseMatrix64F;

import static dr.math.matrixAlgebra.missingData.MissingOps.wrap;

/**
 * @author Marc A. Suchard
 */
public abstract class AbstractValuesViaFullConditionalDelegate extends TipFullConditionalDistributionDelegate {

    final private PartiallyMissingInformation missingInformation;

    protected boolean isLoggable() {
        return false;
    }

    public AbstractValuesViaFullConditionalDelegate(String name, Tree tree,
                                                    MultivariateDiffusionModel diffusionModel,
                                                    ContinuousTraitPartialsProvider dataModel,
                                                    ConjugateRootTraitPrior rootPrior,
                                                    ContinuousRateTransformation rateTransformation,
//                                                    BranchRateModel rateModel,
                                                    ContinuousDataLikelihoodDelegate likelihoodDelegate) {
        super(name, tree, diffusionModel, dataModel, rootPrior, rateTransformation, likelihoodDelegate);
        missingInformation = new PartiallyMissingInformation(tree, dataModel);
    }

    @Override
    protected double[] getTraitForNode(NodeRef node) {

        assert simulationProcess != null;
        assert node != null;

        final int nodeBuffer = likelihoodDelegate.getActiveNodeIndex(node.getNumber());

        if (node.getNumber() >= tree.getExternalNodeCount()) {   // Not external node
            return new double[0];
//                return new MeanAndVariance(new double[0]);
        }

        double[] conditionalNodeBuffer = null; //new double[dimPartial * numTraits];
        likelihoodDelegate.getPostOrderPartial(node.getNumber(), partialNodeBuffer);

        final double[] sample = new double[dimTrait * numTraits];

        int partialOffset = 0;
        int sampleOffset = 0;

        for (int trait = 0; trait < numTraits; ++trait) {
            if (missingInformation.isPartiallyMissing(node.getNumber(), trait)) {
                if (conditionalNodeBuffer == null) {
                    conditionalNodeBuffer = new double[dimPartial * numTraits];

                    simulationProcess.cacheSimulatedTraits(node);
                    cdi.getPreOrderPartial(nodeBuffer, conditionalNodeBuffer);
                }

                if (DEBUG)
                    System.err.println("Missing tip = " + node.getNumber() + " (" + nodeBuffer + "), trait = " + trait);


                final WrappedVector preMean = new WrappedVector.Raw(conditionalNodeBuffer, partialOffset, dimTrait);
                final DenseMatrix64F preVar = wrap(conditionalNodeBuffer, partialOffset + dimTrait + dimTrait * dimTrait, dimTrait, dimTrait);

                final WrappedVector postObs = new WrappedVector.Raw(partialNodeBuffer, partialOffset, dimTrait);

                if (DEBUG) {
                    System.err.println("post: " + postObs);
                    System.err.println("pre : " + preMean);
                    System.err.println("V: " + preVar);
                }


                if (!missingInformation.isCompletelyMissing(node.getNumber(), trait)) {

                    final PartiallyMissingInformation.HashedIntArray intArray =
                            missingInformation.getMissingIndices(node.getNumber(), trait);
                    final int[] missing = intArray.getArray();
                    final int[] observed = intArray.getComplement();

                    ConditionalVarianceAndTransform2 transform =
                            new ConditionalVarianceAndTransform2(
                                    preVar, missing, observed
                            );

                    final WrappedVector cM = transform.getConditionalMean(
                            partialNodeBuffer, partialOffset,      // Tip value
                            conditionalNodeBuffer, partialOffset); // Mean value

                    computeValueWithMissing(cM, // input mean
                            transform.getConditionalCholesky(), // input variance,
                            new WrappedVector.Indexed(sample, sampleOffset, missing), // output sample
                            transform.getTemporaryStorage());

                    if (DEBUG) {
                        System.err.println("cM: " + cM);
                        System.err.println("CV: " + transform.getConditionalVariance());
                        System.err.println("value: " + new WrappedVector.Raw(sample, sampleOffset, dimTrait));
                    }

                }

            } else {
                computeValueWithNoMissing(partialNodeBuffer, partialOffset, sample, sampleOffset, dimTrait);
            }

            partialOffset += dimPartial;
            sampleOffset += dimTrait;
        }

        return sample;
//            return new MeanAndVariance(sample);
    }

    abstract protected void computeValueWithNoMissing(final double[] mean, final int meanOffset,
                                                      final double[] output, final int outputOffset,
                                                      final int dim);

    abstract protected void computeValueWithMissing(final WrappedVector mean,
                                                    final double[][] cholesky,
                                                    final WrappedVector output,
                                                    final double[] buffer);

    private static Boolean DEBUG = false;
}
