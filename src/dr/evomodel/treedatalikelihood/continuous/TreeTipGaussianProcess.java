/*
 * TreeTipGradient.java
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

package dr.evomodel.treedatalikelihood.continuous;

import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.evomodel.treedatalikelihood.ProcessSimulationDelegate;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.continuous.cdi.PrecisionType;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.math.distributions.GaussianProcessRandomGenerator;

/**
 * @author Marc A. Suchard
 */
public class TreeTipGaussianProcess implements GaussianProcessRandomGenerator {

    private final String traitName;
    private final ContinuousDataLikelihoodDelegate likelihoodDelegate;
    private final TreeDataLikelihood treeDataLikelihood;
    private final TreeTrait<double[]> treeTraitProvider;
    private final Tree tree;
    private final Parameter traitParameter;

    private final int nTaxa;
    private final int nTraits;
    private final int dimTrait;
    private final int dimPartial;

    private final Parameter maskParameter;

//    private final PartiallyMissingInformation missingInformation;

//    private final boolean missingOnlyGradient;

    public TreeTipGaussianProcess(String traitName,
                                  TreeDataLikelihood treeDataLikelihood,
                                  ContinuousDataLikelihoodDelegate likelihoodDelegate,
                                  Parameter maskParameter) {
        this.traitName = traitName;
        this.treeDataLikelihood = treeDataLikelihood;
        this.likelihoodDelegate = likelihoodDelegate;
        this.tree = treeDataLikelihood.getTree();
//        this.traitParameter = traitParameter;
        this.maskParameter = maskParameter;

//        this.missingOnlyGradient = missingOnly;

        System.err.println(likelihoodDelegate.getIntegrator().getClass().getCanonicalName());
        System.exit(-1);


        assert(treeDataLikelihood != null);

        String name =
                ProcessSimulationDelegate.TipGradientViaFullConditionalDelegate.getTraitName(traitName);

        System.err.println("name: " + name);
        System.exit(-1);

        treeTraitProvider = treeDataLikelihood.getTreeTrait(name);
        nTaxa = treeDataLikelihood.getTree().getExternalNodeCount();
        nTraits = treeDataLikelihood.getDataLikelihoodDelegate().getTraitCount();
        dimTrait = treeDataLikelihood.getDataLikelihoodDelegate().getTraitDim();

//        missingInformation = new PartiallyMissingInformation(treeDataLikelihood.getTree(),
//                likelihoodDelegate.getDataModel(), likelihoodDelegate);

        PrecisionType precisionType = likelihoodDelegate.getPrecisionType();
        dimPartial = precisionType.getMatrixLength(dimTrait);

        if (precisionType != PrecisionType.SCALAR) {
            throw new RuntimeException("Not yet implemented for full precision");
        }

        if (nTraits != 1) {
            throw new RuntimeException("Not yet implemented for >1 traits");
        }


        this.traitParameter = likelihoodDelegate.getDataModel().getParameter();

        if (maskParameter != null &&
                (maskParameter.getDimension() != traitParameter.getDimension())) {
            throw new RuntimeException("Trait and mask parameters must be the same size");
        }
    }

    @Override
    public Likelihood getLikelihood() {
        return treeDataLikelihood;
    }

//    @Override
//    public Parameter getParameter() {
//        return traitParameter;
//    }

    @Override
    public int getDimension() {
        return 0;
    }

    @Override
    public double[][] getPrecisionMatrix() {
        return new double[0][];
    }

//    @Override
//    public double[] getGradientLogDensity() {
//
//        double[] gradient = new double[nTaxa  * dimTrait * nTraits];
//
//        int offsetOutput = 0;
//        for (int taxon = 0; taxon < nTaxa; ++taxon) {
//            double[] taxonGradient = treeTraitProvider.getTrait(tree, tree.getExternalNode(taxon));
//            System.arraycopy(taxonGradient, 0, gradient, offsetOutput, taxonGradient.length);
//            offsetOutput += taxonGradient.length;
//        }
//
//        if (maskParameter != null) {
//            for (int i = 0; i < maskParameter.getDimension(); ++i) {
//                if (maskParameter.getParameterValue(i) == 0.0) {
//                    gradient[i] = 0.0;
//                }
//            }
//        }
//
//        return gradient;
//    }

    @Override
    public Object nextRandom() {
        return null;
    }

    @Override
    public double logPdf(Object x) {
        return 0;
    }
}
