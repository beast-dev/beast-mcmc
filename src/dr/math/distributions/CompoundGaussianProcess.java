/*
 * CompoundGaussianProcess.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.math.distributions;

import dr.inference.model.CompoundLikelihood;
import dr.inference.model.Likelihood;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Marc A. Suchard
 */

public class CompoundGaussianProcess implements GaussianProcessRandomGenerator {

    private final List<GaussianProcessRandomGenerator> gpList;
//    private final List<Likelihood> likelihoodList;
    private final CompoundLikelihood compoundLikelihood;

    public CompoundGaussianProcess(List<GaussianProcessRandomGenerator> gpList, List<Likelihood> likelihoodList) {
        this.gpList = gpList;
//        this.likelihoodList = likelihoodList;
        compoundLikelihood = new CompoundLikelihood(likelihoodList);
    }

    @Override
    public Likelihood getLikelihood() { return compoundLikelihood; }

    @Override
    public Object nextRandom() {

        int size = 0;
        List<double[]> randomList = new ArrayList<double[]>();
        for (GaussianProcessRandomGenerator gp : gpList) {
            double[] vector = (double[]) gp.nextRandom();
            randomList.add(vector);
            size += vector.length;
//            System.err.println("Drew len = " + vector.length);
        }

        double[] result = new double[size];
        int offset = 0;
        for (double[] vector : randomList) {
            System.arraycopy(vector, 0, result, offset, vector.length);
            offset += vector.length;
        }

//        System.err.println("Done with draw\n");

        return result;
    }

    @Override
    public double logPdf(Object x) {
        throw new RuntimeException("Not yet implemented");
    }
}
