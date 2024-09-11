/*
 * AminoAcidMixture.java
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

import java.util.List;

public class AminoAcidMixture extends DesignMatrix {

    public AminoAcidMixture(List<EmpiricalAminoAcidModel> rateMatrix) {
        super(getName(rateMatrix), getParameters(rateMatrix), false);
        if ( anyRatesAreZero(rateMatrix)) { throw new RuntimeException("AminoAcidMixtureModel cannot be used for rate matrices which have entries which are 0.0.");}
    }

    private static String getName(List<EmpiricalAminoAcidModel> rateMatrix) {
        StringBuilder sb = new StringBuilder("designMatrix");
        for (EmpiricalAminoAcidModel matrix : rateMatrix) {
            String id = matrix.getId();
            if (id != null && id.length() > 0) {
                sb.append(".").append(id);
            }
        }
        return sb.toString();
    }

    private static Parameter[] getParameters(List<EmpiricalAminoAcidModel> rateMatrix) {
        Parameter[] p = new Parameter[rateMatrix.size()];

        for (int i = 0; i < p.length; ++i) {
            EmpiricalAminoAcidModel model = rateMatrix.get(i);
            // this is a flat representation of the lower triangular order
            double[] rates = model.getEmpiricalRateMatrix().getEmpiricalRates();

            // We need upper and lower triangular rates, both logged
            double[] ratesUpperLower = new double[2*rates.length];
            for ( int j = 0; j < rates.length; j++ ) {
                ratesUpperLower[j] = ratesUpperLower[rates.length + j] = Math.log(rates[j]);
            }
            p[i] = new Parameter.Default(ratesUpperLower);
        }
        return p;
    }

    private boolean anyRatesAreZero(List<EmpiricalAminoAcidModel> rateMatrix) {
        boolean anyZeros = false;

        for (int i = 0; i < rateMatrix.size(); ++i) {
            EmpiricalAminoAcidModel model = rateMatrix.get(i);
            // this is a flat representation of the lower triangular order
            double[] rates = model.getEmpiricalRateMatrix().getEmpiricalRates();
            for ( int j = 0; j < rates.length; j++ ) {
                if ( rates[j] < Double.MIN_VALUE) {
                    anyZeros = true;
                }
            }
        }
        return anyZeros;
    }
}
