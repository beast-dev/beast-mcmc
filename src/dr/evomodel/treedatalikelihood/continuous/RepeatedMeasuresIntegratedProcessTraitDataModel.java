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

import dr.evolution.tree.MutableTreeModel;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.treedatalikelihood.continuous.cdi.PrecisionType;
import dr.evomodel.treedatalikelihood.preorder.ContinuousExtensionDelegate;
import dr.evomodelxml.treelikelihood.TreeTraitParserUtilities;
import dr.inference.model.CompoundParameter;
import dr.inference.model.MatrixParameterInterface;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.math.matrixAlgebra.CholeskyDecomposition;
import dr.math.matrixAlgebra.IllegalDimension;
import dr.math.matrixAlgebra.Matrix;
import dr.math.matrixAlgebra.WrappedVector;
import dr.math.matrixAlgebra.missingData.MissingOps;
import dr.xml.*;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;

import java.util.Arrays;

/**
 * @author Marc A. Suchard
 * @author Gabriel Hassler
 * @author Paul Bastide
 */
public class RepeatedMeasuresIntegratedProcessTraitDataModel extends RepeatedMeasuresTraitDataModel {

    private IntegratedProcessTraitDataModel integratedProcessDataModel = null;
    private int dimProcess;

    public RepeatedMeasuresIntegratedProcessTraitDataModel(String name,
                                                           CompoundParameter parameter,
                                                           boolean[] missindIndicators,
                                                           boolean useMissingIndices,
                                                           final int dimTrait,
                                                           MatrixParameterInterface samplingPrecision) {

        super(name, parameter, missindIndicators, useMissingIndices, dimTrait, samplingPrecision, PrecisionType.FULL);

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

        scalePartialwithSamplingPrecision(partial, dimProcess);

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
    protected int getParameterPartialDimension() {
        return 2 * getParameter().getDimension();
    }

}

