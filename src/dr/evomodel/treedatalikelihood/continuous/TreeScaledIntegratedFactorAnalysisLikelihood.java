/*
 * TreeScaledIntegratedFactorAnalysisLikelihood.java
 *
 * Copyright (c) 2002-2019 Alexei Drummond, Andrew Rambaut and Marc Suchard
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
import dr.evomodel.continuous.hmc.TaxonTaskPool;
import dr.inference.model.CompoundParameter;
import dr.inference.model.MatrixParameterInterface;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;

import java.util.List;

/**
 * @author Marc A. Suchard
 * @author Paul Bastide
 */

public class TreeScaledIntegratedFactorAnalysisLikelihood extends IntegratedFactorAnalysisLikelihood {

    public TreeScaledIntegratedFactorAnalysisLikelihood(String name,
                                                        CompoundParameter traitParameter,
                                                        List<Integer> missingIndices,
                                                        MatrixParameterInterface loadings,
                                                        MatrixParameterInterface traitPrecision,
                                                        double nuggetPrecision,
                                                        TaxonTaskPool taxonTaskPool) {
        super(name, traitParameter, missingIndices, loadings, traitPrecision, nuggetPrecision, taxonTaskPool);
    }


//    @Override
//    protected void storeState() {
//        super.storeState();
//        System.arraycopy(tipHeights, 0, storedTipHeights, 0, tipHeights.length);
//
//    }

//    @Override
//    protected void restoreState() {
//        super.restoreState();
//
//        double[] tmp3 = tipHeights;
//        tipHeights = storedTipHeights;
//        storedTipHeights = tmp3;
//    }


    // Private class functions

//    @Override
//    protected void setupStatistics() {
//
//        if (tipHeights == null) {
//            tipHeights = new double[numTaxa];
//            storedTipHeights = new double[numTaxa];
//        }
//
//        super.setupStatistics();
//    }


    @Override
    protected DenseMatrix64F getTraitVariance(final int taxon) {

        DenseMatrix64F V = super.getTraitVariance(taxon);
        if (delegate != null) CommonOps.scale(getTipHeight(taxon), V);
        return V;
    }

    private double getTipHeight(int taxon) {
        Tree tree = delegate.getCallbackLikelihood().getTree();
        double time = tree.getNodeHeight(tree.getRoot()) - tree.getNodeHeight(tree.getExternalNode(taxon));
        return time * delegate.getRateTransformation().getNormalization();
    }


//    private double[] tipHeights;
//    private double[] storedTipHeights;

    @Override
    public void setLikelihoodDelegate(ContinuousDataLikelihoodDelegate delegate) {
        super.setLikelihoodDelegate(delegate);
        makeDirty();
        fireModelChanged();
    }

    @Override
    protected void fillTreeScales(double[] tmp, double[][] treeSharedLengths) {
        for (int i = 0; i < tmp.length; i++) {
            tmp[i] = treeSharedLengths[i][i];
        }
    }

}
