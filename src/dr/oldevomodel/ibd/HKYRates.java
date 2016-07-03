/*
 * HKYRates.java
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

package dr.oldevomodel.ibd;

import dr.oldevomodel.substmodel.FrequencyModel;
import dr.oldevomodel.substmodel.HKY;
import dr.inference.model.Parameter;

/**
 * Package: dr.evomodel.ibd
 * Description:
 * <p/>
 * <p/>
 * Created by
 * avaleks (alexander.alekseyenko@gmail.com)
 * Date: 05-Aug-2008
 * Time: 09:32:35
 */
public class HKYRates extends HKY {
    /**
     * Constructor
     */
    public HKYRates(Parameter kappaParameter, FrequencyModel freqModel) {
        super(kappaParameter, freqModel);
    }

    public double[] getRelativeRates(double[] rateMatrix) {
        double kappa = getKappa();
        double[] freq = getFrequencyModel().getFrequencies();
        // A - C - G - T
        rateMatrix[0] = -(freq[1] + freq[3]) - freq[2] * kappa;
        rateMatrix[1] = freq[1];
        rateMatrix[2] = freq[2] * kappa;
        rateMatrix[3] = freq[3];

        rateMatrix[4] = freq[0];
        rateMatrix[5] = -(freq[0] + freq[2]) - freq[3] * kappa;
        rateMatrix[6] = freq[2];
        rateMatrix[7] = freq[3] * kappa;

        rateMatrix[8] = freq[0] * kappa;
        rateMatrix[9] = freq[1];
        rateMatrix[10] = -(freq[1] + freq[3]) - freq[0] * kappa;
        rateMatrix[11] = freq[3];

        rateMatrix[12] = freq[0];
        rateMatrix[13] = freq[1] * kappa;
        rateMatrix[14] = freq[2];
        rateMatrix[15] = -(freq[0] + freq[1]) - freq[1] * kappa;


        return rateMatrix;
    }
}
