/*
 * SubstitutionRateMatrixMixture.java
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

package dr.evomodel.substmodel;

import dr.evomodel.substmodel.aminoacid.EmpiricalAminoAcidModel;
import dr.inference.model.DesignMatrix;
import dr.inference.model.Parameter;
import dr.math.matrixAlgebra.Vector;
import dr.math.matrixAlgebra.Matrix;

import java.util.List;

public class SubstitutionRateMatrixMixture extends DesignMatrix {

    public SubstitutionRateMatrixMixture(List<SubstitutionModel> rateMatrixList) {
        super(getName(rateMatrixList), getParameters(rateMatrixList), false);
        if (!modelsAreSameSize(rateMatrixList)) {
            throw new RuntimeException("SubstitutionRateMatrixMixture requires all component models have same dimension.");
        }
    }

    private static String getName(List<SubstitutionModel> rateMatrixList) {
        StringBuilder sb = new StringBuilder("designMatrix");
        for (SubstitutionModel model : rateMatrixList) {
            String id = model.getId();
            if (id != null && id.length() > 0) {
                sb.append(".").append(id);
            }
        }
        return sb.toString();
    }

    private static Parameter[] getParameters(List<SubstitutionModel> rateMatrixList) {
        Parameter[] p = new Parameter[rateMatrixList.size()];

        int alphabetSize = rateMatrixList.get(0).getFrequencyModel().getFrequencyCount();
        int fullSize = alphabetSize * alphabetSize;
        int size = alphabetSize * (alphabetSize - 1);
        int halfSize = size/2;

        for (int m = 0; m < p.length; ++m) {
            // Get Q-matrix
            SubstitutionModel model = rateMatrixList.get(m);
            // Full matrix stored row-wise
            double[] qMatrix = new double[fullSize];
            model.getInfinitesimalMatrix(qMatrix);

            // get base freqs
            double[] freqs = model.getFrequencyModel().getFrequencies();
            for (int i = 0; i < alphabetSize; i++) {
                freqs[i] = Math.log(freqs[i]);
            }

            double[] rates = new double[size];
            int idx = 0;
            for (int i = 0; i < (alphabetSize - 1); i++) {
                for (int j = (i + 1); j < alphabetSize; j++) {
                    rates[idx] = Math.log(qMatrix[i*alphabetSize + j]) - freqs[j];
                    rates[idx + halfSize] = Math.log(qMatrix[j*alphabetSize + i]) - freqs[i];
                    idx++;
                }
            }

            p[m] = new Parameter.Default(rates);
        }
        return p;
    }

//   // TODO if this class is to subsume aminoAcidMixtureModel this will need to be made to work
//    private boolean anyRatesAreZero(List<SubstitutionModel> rateMatrixList) {
//        boolean anyZeros = false;
//        int s = rateMatrixList.get(0).getFrequencyModel().getFrequencyCount();
//        int s2 = s * s;
//        for (int i = 0; i < rateMatrixList.size(); ++i) {
//            SubstitutionModel model = rateMatrixList.get(i);
//            double[] qMatrix = new double[s2];
//            model.getInfinitesimalMatrix(qMatrix);
//            int x = 0;
//            for ( int j = 0; j < qMatrix.length; j++ ) {
//                if ((!(j % s == 0)) && qMatrix[j] < Double.MIN_VALUE) {
//                    anyZeros = true;
//                }
//            }
//        }
//        return anyZeros;
//    }

    private boolean modelsAreSameSize(List<SubstitutionModel> rateMatrixList) {
        boolean sameSize = true;

        int size0 = rateMatrixList.get(0).getFrequencyModel().getFrequencyCount();

        for (int i = 0; i < rateMatrixList.size(); ++i) {
            int size = rateMatrixList.get(i).getFrequencyModel().getFrequencyCount();
            if (size != size0) {
                sameSize = false;
            }
        }
        return sameSize;
    }

}
