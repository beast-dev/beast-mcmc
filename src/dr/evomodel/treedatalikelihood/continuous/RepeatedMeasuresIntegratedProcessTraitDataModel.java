/*
 * RepeatedMeasuresTraitDataModel.java
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
import dr.inference.model.MatrixParameterInterface;
import dr.math.matrixAlgebra.Matrix;

import static dr.math.matrixAlgebra.missingData.MissingOps.blockUnwrap;
import static dr.math.matrixAlgebra.missingData.MissingOps.wrap;

/**
 * @author Marc A. Suchard
 * @author Gabriel Hassler
 * @author Paul Bastide
 */
public class RepeatedMeasuresIntegratedProcessTraitDataModel extends RepeatedMeasuresTraitDataModel {

    private IntegratedProcessTraitDataModel integratedProcessDataModel = null;
    private int dimProcess;

    public RepeatedMeasuresIntegratedProcessTraitDataModel(String name,
                                                           ContinuousTraitPartialsProvider childModel,
                                                           CompoundParameter parameter,
                                                           boolean[] missindIndicators,
                                                           boolean useMissingIndices,
                                                           final int dimTrait,
                                                           final int numTraits,
                                                           MatrixParameterInterface samplingPrecision) {

        super(name, childModel, parameter, missindIndicators, useMissingIndices, dimTrait, numTraits, samplingPrecision, PrecisionType.FULL);

        integratedProcessDataModel = new IntegratedProcessTraitDataModel(name, parameter, missindIndicators, useMissingIndices, dimTrait);
        dimProcess = 2 * dimTrait;
    }

    @Override
    public double[] getTipPartial(int taxonIndex, boolean fullyObserved) {

        assert (numTraits == 1);
        assert (samplingPrecision.rows() == dimProcess && samplingPrecision.columns() == dimProcess);

        recomputeVariance();

        if (fullyObserved) {
            throw new RuntimeException("Incompatible with this model.");
        }

        double[] partial = integratedProcessDataModel.getTipPartial(taxonIndex, fullyObserved);

        scalePartialwithSamplingPrecision(partial, taxonIndex, dimProcess);

        return partial;
    }

    @Override
    public boolean[] getDataMissingIndicators() {
        return integratedProcessDataModel.getDataMissingIndicators();
    }

    @Override
    public int getTraitDimension() {
        return dimProcess;
    }

    @Override
    public boolean isIntegratedProcess() {
        return true;
    }

    @Override
    protected int getParameterPartialDimension() {
        return 2 * getParameter().getDimension();
    }

    @Override
    protected void calculatePrecisionInfo() {
        int dimProcess = 2 * dimTrait;
        double[] precisionBuffer = new double[dimProcess * dimProcess];
        // Precision [Id, 0; 0, precisionPosistion]
        blockUnwrap(wrap(samplingPrecisionParameter.getParameterValues(), 0, dimTrait, dimTrait),
                precisionBuffer, 0,
                dimTrait, dimTrait, dimProcess);
        for (int i = 0; i < dimTrait; i++) {
            precisionBuffer[i * dimProcess + i] = 1.0;
        }
        samplingPrecision = new Matrix(precisionBuffer, dimProcess, dimProcess);
    }

}

