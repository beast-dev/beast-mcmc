/*
 * IntegratedOrnsteinUhlenbeckDataModel.java
 *
 * Copyright (c) 2002-2018 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

import dr.evomodel.treedatalikelihood.continuous.cdi.PrecisionType;
import dr.inference.model.CompoundParameter;

import java.util.List;

import static dr.math.matrixAlgebra.missingData.MissingOps.blockUnwrap;
import static dr.math.matrixAlgebra.missingData.MissingOps.wrap;

/**
 * @author Marc A. Suchard
 * @author Paul Bastide
 */

public class IntegratedProcessTraitDataModel extends
        ContinuousTraitDataModel implements ContinuousTraitPartialsProvider {


    public IntegratedProcessTraitDataModel(String name,
                                           CompoundParameter parameter,
                                           boolean[] missingIndicators,
                                           boolean useMissingIndices,
                                           final int dimTrait) {
        super(name, parameter, missingIndicators, useMissingIndices, dimTrait, PrecisionType.FULL);
    }

    public IntegratedProcessTraitDataModel(String name,
                                           CompoundParameter parameter,
                                           boolean[] missingIndicators,
                                           boolean useMissingIndices,
                                           final int dimTrait,
                                           PrecisionType precisionType) {
        this(name, parameter, missingIndicators, useMissingIndices, dimTrait);
        assert precisionType == PrecisionType.FULL : "Integrated Process is only implemented for full precision type.";
    }

    @Override
    public double[] getTipPartial(int taxonIndex, boolean fullyObserved) {

        assert numTraits == 1;
        assert precisionType == PrecisionType.FULL;

        double[] partial = super.getTipPartial(taxonIndex, fullyObserved);

        int dimTraitDouble = 2 * dimTrait;
        int dimPartialDouble = precisionType.getPartialsDimension(dimTraitDouble);
        double[] partialDouble = new double[dimPartialDouble];

        // Traits [0, traitsPosition]
        System.arraycopy(partial, 0, partialDouble, dimTrait, dimTrait);
        // Precision [0, 0; 0, precisionPosistion]
        blockUnwrap(wrap(partial, precisionType.getPrecisionOffset(dimTrait), dimTrait, dimTrait),
                partialDouble, precisionType.getPrecisionOffset(dimTraitDouble),
                dimTrait, dimTrait, dimTraitDouble);
        // Variance [Inf, 0; 0, variancePosisition]
        int offsetVar = precisionType.getVarianceOffset(dimTraitDouble);
        for (int i = 0; i < dimTrait; i++) {
            partialDouble[offsetVar + i * dimTraitDouble + i] = Double.POSITIVE_INFINITY;
        }
        blockUnwrap(wrap(partial, precisionType.getVarianceOffset(dimTrait), dimTrait, dimTrait),
                partialDouble, offsetVar,
                dimTrait, dimTrait, dimTraitDouble);

        int effDimOffset = precisionType.getEffectiveDimensionOffset(dimTrait);
        int effDim = (int) Math.round(partial[effDimOffset]);
        precisionType.fillEffDimInPartials(partialDouble, 0, effDim, dimTraitDouble);

        int detDimOffset = precisionType.getDeterminantOffset(dimTrait);
        double det = partial[detDimOffset];
        precisionType.fillDeterminantInPartials(partialDouble, 0, det, dimTraitDouble);


        return partialDouble;
    }

    @Override
    public int getTraitDimension() {
        return 2 * dimTrait;
    }

}