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

import java.util.List;

/**
 * @author Marc A. Suchard
 */
public class EmptyTraitDataModel implements ContinuousTraitPartialsProvider {

    private final String name;
    private final int dimTrait;
    private final PrecisionType precisionType;

    public EmptyTraitDataModel(String name, final int dimTrait, PrecisionType precisionType) {
        this.name = name;
        this.dimTrait = dimTrait;
        this.precisionType = precisionType;
    }

    @Override
    public boolean bufferTips() { return true; } // TODO maybe should be false

    @Override
    public int getTraitCount() {  return 1; }

    @Override
    public int getTraitDimension() { return dimTrait; }

    @Override
    public PrecisionType getPrecisionType() {
        return precisionType;
    }

    @Override
    public CompoundParameter getParameter() { return null; }

    @Override
    public String getModelName() {
        return name;
    }

    @Override
    public List<Integer> getMissingIndices() { return null; }

    @Override
    public double[] getTipPartial(int taxonIndex, boolean fullyObserved) {
        return new double[dimTrait + precisionType.getMatrixLength(dimTrait)];
    }
}
