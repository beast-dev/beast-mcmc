/*
 * MultiEpochExponentialTest.java
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

package test.dr.evomodel.coalescent;

import dr.evolution.coalescent.ExponentialExponential;
import dr.evolution.coalescent.ExponentialGrowth;
import dr.evolution.coalescent.MultiEpochExponential;
import dr.evolution.util.Units;
import dr.evomodel.coalescent.smooth.SmoothSkygridLikelihood;
import junit.framework.TestCase;

/**
 * @author Marc A. Suchard
 */
public class SmoothSkygridTest extends TestCase {

    public void testSmoothSkyrideIntegrals() {


        double x = 0.0;


//        x = SmoothSkygridLikelihood.getLogPopSizeInInterval(0.1, 0, 1, 2.0);

        x = SmoothSkygridLikelihood.getPopSizeInInterval(0.6,
                0, 1,
                Math.log(5), Math.log(10), 10.0);

        System.err.println(x);

        x = SmoothSkygridLikelihood.getIntensityInInterval(0.6,
                0, 1,
                Math.log(5), Math.log(10), 10.0);

        System.err.println(x);

    }

    private final static double tolerance1 = 1E-10;
    private final static double tolerance2 = 1E-5;
}
