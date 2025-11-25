/*
 * TreeTipGradient.java
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

package dr.evomodel.treedatalikelihood.continuous;

import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.inference.model.MaskedParameter;
import dr.inference.model.Parameter;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Marc A. Suchard
 */
public class SubsetTreeTipGradient extends TreeTipGradient {

    private final MaskedParameter maskedParameter;
    private final List<MaskMap> maskMapList;

    public SubsetTreeTipGradient(String traitName,
                                 MaskedParameter maskedParameter,
                                 TreeDataLikelihood treeDataLikelihood,
                                 ContinuousDataLikelihoodDelegate likelihoodDelegate) {

        super(traitName, maskedParameter.getUnmaskedParameter(), treeDataLikelihood, likelihoodDelegate, null);
        this.maskedParameter = maskedParameter;

        this.maskMapList = makeMaskMapList();
    }

    class MaskMap {
        int nodeIndex;
        int[] dimIndices;

        MaskMap(int nodeIndex, int[] dimIndices) {
            this.nodeIndex = nodeIndex;
            this.dimIndices = dimIndices;
        }
    }

    private List<MaskMap> makeMaskMapList() {

        List<MaskMap> maskMapList = new ArrayList<>();
        int offset = 0;
        for (int taxon = 0; taxon < nTaxa; ++taxon) {
            List<Integer> indices = new ArrayList<>();
            for (int i = 0; i < dimTrait * nTraits; ++i) {
                if (maskedParameter.getParameterMaskValue(offset) == 1.0) {
                    indices.add(i);
                }
                ++offset;
            }
            if (indices.size() > 0) {
                maskMapList.add(new MaskMap(taxon, indices.stream()
                        .mapToInt(Integer::intValue)
                        .toArray()));
            }
        }

        return maskMapList;
    }

    @Override
    public Parameter getParameter() {
        return maskedParameter;
    }

    @Override
    public double[] getGradientLogDensity() {

        double[] gradient = new double[getDimension()];

        int offsetOutput = 0;
        for (MaskMap map : maskMapList) {
            double[] taxonGradient = (double[]) treeTraitProvider.getTrait(
                    tree, tree.getExternalNode(map.nodeIndex));
            for (int i : map.dimIndices) {
                gradient[offsetOutput] = taxonGradient[i];
                ++offsetOutput;
            }
        }


//
//        for (int taxon = 0; taxon < nTaxa; ++taxon) {
//            double[] taxonGradient = (double[]) treeTraitProvider.getTrait(tree, tree.getExternalNode(taxon));
//            System.arraycopy(taxonGradient, 0, gradient, offsetOutput, dimTrait);
//            offsetOutput += dimTrait;
//        }

        return gradient;
    }

    @Override
    double getTolerance() {
        return 1E-2;
    }
}
