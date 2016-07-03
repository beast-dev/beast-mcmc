/*
 * TN93Test.java
 *
 * Copyright (c) 2002-2012 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package test.dr.app.beagle;

import dr.evomodel.substmodel.EigenDecomposition;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.nucleotide.TN93;
import dr.evolution.datatype.Nucleotides;
import dr.inference.model.Parameter;
import dr.math.matrixAlgebra.Vector;
import test.dr.math.MathTestCase;

/**
 * @author Marc Suchard
 */
public class TN93Test extends MathTestCase {

    public void testTN93() {

        Parameter kappa1 = new Parameter.Default(5.0);
        Parameter kappa2 = new Parameter.Default(2.0);
        double[] pi = new double[]{0.40, 0.20, 0.30, 0.10};
        double time = 0.1;

        FrequencyModel freqModel = new FrequencyModel(Nucleotides.INSTANCE, pi);
        TN93 tn = new TN93(kappa1, kappa2, freqModel);

        EigenDecomposition decomp = tn.getEigenDecomposition();

        Vector eval = new Vector(decomp.getEigenValues());
        System.out.println("Eval = " + eval);

        double[] probs = new double[16];
        tn.getTransitionProbabilities(time, probs);
        System.out.println("new probs = " + new Vector(probs));

        // check against old implementation
        dr.oldevomodel.substmodel.FrequencyModel oldFreq = new dr.oldevomodel.substmodel.FrequencyModel(Nucleotides.INSTANCE, pi);
        dr.oldevomodel.substmodel.TN93 oldTN = new dr.oldevomodel.substmodel.TN93(kappa1, kappa2, oldFreq);

        double[] oldProbs = new double[16];
        oldTN.getTransitionProbabilities(time, oldProbs);
        System.out.println("old probs = " + new Vector(oldProbs));
        assertEquals(probs, oldProbs, 10E-6);
    }

}
