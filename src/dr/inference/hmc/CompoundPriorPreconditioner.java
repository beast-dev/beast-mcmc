/*
 * CompoundPriorPreconditioner.java
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

package dr.inference.hmc;

import dr.inference.model.PriorPreconditioningProvider;
import dr.math.matrixAlgebra.Vector;
import dr.xml.Reportable;

import java.util.Arrays;
import java.util.List;

/**
 * @author Alexander Fisher
 * @author Andy Magee
 */

public class CompoundPriorPreconditioner implements PriorPreconditioningProvider, Reportable {
// todo: look at refactoring into MassPreconditioner.CompoundPreconditioning
    private final int smallDim;
    private final int totalDim;
    private List<PriorPreconditioningProvider> priorPreconditionerList;
    private int[] iMap;
    private int[] jMap;

    public CompoundPriorPreconditioner(List<PriorPreconditioningProvider> priorPreconditionerList) {
        this.smallDim = priorPreconditionerList.get(0).getDimension();
        this.priorPreconditionerList = priorPreconditionerList;

        int tempDim = 0;
        for (int i = 0; i < priorPreconditionerList.size(); i++) {
            tempDim += priorPreconditionerList.get(i).getDimension();
        }
        this.totalDim = tempDim;

        int[] mapI = new int[totalDim];
        int[] mapJ = new int[totalDim];
        int k = 0;
        for (int i = 0; i < priorPreconditionerList.size(); i++) {
            for (int j = 0; j < priorPreconditionerList.get(i).getDimension(); j++) {
                mapI[k] = i;
                mapJ[k] = j;
                k++;
            }
        }
        this.iMap = mapI;
        this.jMap = mapJ;
    }

    @Override
    public double getStandardDeviation(int index) {
//        int derivativeIndex = (int) Math.floor(index / smallDim);
//        int standardDeviationIndex = index % smallDim;
//        return priorPreconditionerList.get(derivativeIndex).getStandardDeviation(standardDeviationIndex);
        int i = iMap[index];
        int j = jMap[index];
        return priorPreconditionerList.get(i).getStandardDeviation(j);
    }

    @Override
    public int getDimension() {
        return totalDim;
    }

    public String getReport() {

        StringBuilder sb = new StringBuilder("compoundPriorPreconditioner Report\n\n");

        sb.append("totalDim: " + totalDim + "\n\n");

        sb.append("priorPreconditionerList size: " + priorPreconditionerList.size() + "\n\n");

        int n = 0;
        for (int i = 0; i < priorPreconditionerList.size(); i++) {
            n += priorPreconditionerList.get(i).getDimension();
        }

        double[] values = new double[n];
        for (int i = 0; i < n; i++) {
            values[i] = getStandardDeviation(i);
        }
        sb.append("Prior SDs: ");
        sb.append(new Vector(values));
        sb.append("\n\n");

        return sb.toString();
    }
}
