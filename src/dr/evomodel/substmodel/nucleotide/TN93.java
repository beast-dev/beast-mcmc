/*
 * TN93.java
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

package dr.evomodel.substmodel.nucleotide;

import dr.evomodel.substmodel.BaseSubstitutionModel;
import dr.evomodel.substmodel.EigenDecomposition;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evolution.datatype.Nucleotides;
import dr.inference.model.Parameter;
import dr.util.Author;
import dr.util.Citable;
import dr.util.Citation;

import java.util.Collections;
import java.util.List;

/**
 * Tamura-Nei model of nucleotide evolution
 *
 * @author Marc A. Suchard
 */
public class TN93 extends BaseSubstitutionModel implements Citable {

    private Parameter kappaParameter1 = null;
    private Parameter kappaParameter2 = null;

    public TN93(Parameter kappaParameter1, Parameter kappaParameter2, FrequencyModel freqModel) {

        super("TN93", Nucleotides.INSTANCE, freqModel);

        this.kappaParameter1 = kappaParameter1;
        addVariable(kappaParameter1);
        kappaParameter1.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));

        this.kappaParameter2 = kappaParameter2;
        addVariable(kappaParameter2);
        kappaParameter2.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));
    }

    public final double getKappa1() {
        return kappaParameter1.getParameterValue(0);
    }

    public final double getKappa2() {
        return kappaParameter2.getParameterValue(0);
    }

    protected void frequenciesChanged() {
    }

    protected void ratesChanged() {
    }

    protected void setupRelativeRates(double[] rates) {
        double kappa1 = getKappa1();
        double kappa2 = getKappa2();
        rates[0] = 1.0;
        rates[1] = kappa1;
        rates[2] = 1.0;
        rates[3] = 1.0;
        rates[4] = kappa2;
        rates[5] = 1.0;
    }

    public synchronized EigenDecomposition getEigenDecomposition() {

        if (eigenDecomposition == null) {
            double[] evec = new double[stateCount * stateCount];
            double[] ivec = new double[stateCount * stateCount];
            double[] eval = new double[stateCount];
            eigenDecomposition = new EigenDecomposition(evec, ivec, eval);

            ivec[2 * stateCount + 1] = 1; // left eigenvectors 3 = (0,1,0,-1); 4 = (1,0,-1,0)
            ivec[2 * stateCount + 3] = -1;

            ivec[3 * stateCount + 0] = 1;
            ivec[3 * stateCount + 2] = -1;

            evec[0 * stateCount + 0] = 1; // right eigenvector 1 = (1,1,1,1)'
            evec[1 * stateCount + 0] = 1;
            evec[2 * stateCount + 0] = 1;
            evec[3 * stateCount + 0] = 1;

        }

        if (updateMatrix) {

            double[] evec = eigenDecomposition.getEigenVectors();
            double[] ivec = eigenDecomposition.getInverseEigenVectors();
            double[] pi = freqModel.getFrequencies();
            double piR = pi[0] + pi[2];
            double piY = pi[1] + pi[3];

            // left eigenvector #1
            ivec[0 * stateCount + 0] = pi[0]; // or, evec[0] = pi;
            ivec[0 * stateCount + 1] = pi[1];
            ivec[0 * stateCount + 2] = pi[2];
            ivec[0 * stateCount + 3] = pi[3];

            // left eigenvector #2
            ivec[1 * stateCount + 0] = pi[0] * piY;
            ivec[1 * stateCount + 1] = -pi[1] * piR;
            ivec[1 * stateCount + 2] = pi[2] * piY;
            ivec[1 * stateCount + 3] = -pi[3] * piR;

            // right eigenvector #2
            evec[0 * stateCount + 1] = 1.0 / piR;
            evec[1 * stateCount + 1] = -1.0 / piY;
            evec[2 * stateCount + 1] = 1.0 / piR;
            evec[3 * stateCount + 1] = -1.0 / piY;

            // right eigenvector #3
            evec[1 * stateCount + 2] = pi[3] / piY;
            evec[3 * stateCount + 2] = -pi[1] / piY;

            // right eigenvector #4
            evec[0 * stateCount + 3] = pi[2] / piR;
            evec[2 * stateCount + 3] = -pi[0] / piR;

            // eigenvectors
            double[] eval = eigenDecomposition.getEigenValues();

            final double kappa1 = getKappa1();
            final double kappa2 = getKappa2();
            final double beta = -1.0 / (2.0 * (piR * piY + kappa1 * pi[0] * pi[2] + kappa2 * pi[1] * pi[3]));
            final double A_R = 1.0 + piR * (kappa1 - 1);
            final double A_Y = 1.0 + piY * (kappa2 - 1);

            eval[1] = beta;
            eval[2] = beta * A_Y;
            eval[3] = beta * A_R;

            updateMatrix = false;
        }

        return eigenDecomposition;
    }

    @Override
    public Citation.Category getCategory() {
        return Citation.Category.SUBSTITUTION_MODELS;
    }

    @Override
    public String getDescription() {
        return "Tamura-Nei nucleotide substitution model";
    }

    @Override
    public List<Citation> getCitations() {
        return Collections.singletonList(CITATION);
    }

    public static Citation CITATION = new Citation(
            new Author[]{
                    new Author("K", "Tamura"),
                    new Author("M", "Nei")
            },
            "Estimation of the number of nucleotide substitutions in the control region of mitochondrial DNA in humans and chimpanzees",
            1993,
            "Mol Biol Evol",
            10, 512, 526
    );
}