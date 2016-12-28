/*
 * ContinuousTraitData.java
 *
 * Copyright (c) 2002-2016 Alexei Drummond, Andrew Rambaut and Marc Suchard
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
import dr.inference.model.*;

import java.util.*;

/**
 * @author Marc A. Suchard
 */
public class ContinuousTraitDataModel extends AbstractModel {

    private final CompoundParameter parameter;
    private final List<Integer> missingIndices;

    private final int numTraits;
    private final int dimTrait;
    private final PrecisionType precisionType;

    public ContinuousTraitDataModel(String name,
                                    CompoundParameter parameter,
                                    List<Integer> missingIndices,
                                    final int dimTrait) {
        super(name);
        this.parameter = parameter;
        this.missingIndices = missingIndices;
        addVariable(parameter);

        this.dimTrait = dimTrait;
        this.numTraits = getParameter().getParameter(0).getDimension() / dimTrait;
        this.precisionType = PrecisionType.SCALAR;
    }

    public boolean bufferTips() { return true; }

    public int getTraitCount() {  return numTraits; }

    public int getTraitDimension() { return dimTrait; }

    public PrecisionType getPrecisionType() {
        return precisionType;
    }

    public String getName() { return super.getModelName(); }

    public CompoundParameter getParameter() { return parameter; }

    public List<Integer> getMissingIndices() { return missingIndices; }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {

    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        if (variable == parameter) {
            if (type == Parameter.ChangeType.VALUE_CHANGED) {
                fireModelChanged(this, getTaxonIndex(index));
            } else if (type == Parameter.ChangeType.ALL_VALUES_CHANGED){
                fireModelChanged(this);
            } else {
                throw new RuntimeException("Unhandled parameter change type");
            }
        }
    }

    private int getTaxonIndex(int parameterIndex) {
        return parameterIndex / (dimTrait * numTraits);
    }

    @Override
    protected void storeState() { }

    @Override
    protected void restoreState() { }

    @Override
    protected void acceptState() { }

//    public double[] getTipMean(int index) {
//        return parameter.getParameter(index).getParameterValues();
//    }
//
//    public double[] getTipPrecision(int index) {
//        if (numTraits == 1) {
//            return NON_MISSING;
//        } else {
//            double[] missing = new double[numTraits];
//            Arrays.fill(missing, Double.POSITIVE_INFINITY);
//            return missing;
//        }
//    }

//    public boolean getAnyPartiallyMissing(int taxonIndex) {
//        boolean missing = false;
//
//        return missing;
//    }
//
//    private void buildMissing
//
//    public boolean[] getPartiallyMissing(int taxonIndex) {
//
//        boolean[] missing = new boolean[numTraits * dimTrait];
//        if (missingIndices != null) {
//            for (int i = 0; i < numTraits; ++i) {
//                for (int j = 0; j < dimTrait; ++j) {
//                    final int index = i * dimTrait + j;
//                    final int missingIndex = index + dimTrait * numTraits * taxonIndex;
//                    missing[index] = missingIndices.contains(missingIndex);
//                }
//            }
//        }
//        return missing;
//    }

    public double[] getTipPartial(int taxonIndex) {
        double[] partial = new double[numTraits * (dimTrait + 1)];
        final Parameter p = parameter.getParameter(taxonIndex);
        int offset = 0;
        for (int i = 0; i < numTraits; ++i) {
            boolean missing = false;
            for (int j = 0; j < dimTrait; ++j) {
                final int index = i * dimTrait + j;
                final int missingIndex = index + dimTrait * numTraits * taxonIndex;
                partial[offset + j] = p.getParameterValue(index);
                if (missingIndices != null && missingIndices.contains(missingIndex)) {
                    missing = true;
                }
            }
            partial[offset + dimTrait] = missing ? 0.0 : Double.POSITIVE_INFINITY;
            offset += dimTrait + 1;
        }
        return partial;
    }

    private static double[] NON_MISSING = new double[] { Double.POSITIVE_INFINITY };

    private Map<Integer, boolean[]> missingCache;
    private boolean[] hasAnyMissing;

    /**
     * For partially observed tips: (y_1, y_2)^t \sim N(\mu, \Sigma) where
     *
     *      \mu = (\mu_1, \mu_2)^t
     *      \Sigma = ((\Sigma_{11}, \Sigma_{12}), (\Sigma_{21}, \Sigma_{22})^t
     *
     * then  y_1 | y_2 \sim N (\bar{\mu}, \bar{\Sigma}), where
     *
     *      \bar{\mu} = \mu_1 + \Sigma_{12}\Sigma_{22}^{-1}(y_2 - \mu_2), and
     *      \bar{\Sigma} = \Sigma_{11} - \Sigma_{12}\Sigma_{22}^1\Sigma{21}
     *
     */
}
