/*
 * DependentProductChainSubstitutionModel.java
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

package dr.evomodel.substmodel;

import dr.evomodel.siteratemodel.SiteRateModel;
import dr.inference.model.Parameter;
import dr.math.KroneckerOperation;

import java.util.List;

/**
 * A derived class from ProductChainSubstitutionModel that introduces dependence between sites.
 *
 * @author Marc A. Suchard
 * @author Yu-Nong Gong
 * @version 08/01/2010
 */
public class DependentProductChainSubstitutionModel extends ProductChainSubstitutionModel {

    public DependentProductChainSubstitutionModel(String name, List<SubstitutionModel> baseModels,
                                                  Parameter dependenceParameter) {
        this(name, baseModels, null, dependenceParameter);
    }

    public DependentProductChainSubstitutionModel(String name, List<SubstitutionModel> baseModels,
                                                  List<SiteRateModel> rateModels, Parameter dependenceParameter) {
        super(name, baseModels, rateModels);
        this.dependenceParameter = dependenceParameter;

        eigenSystemCES = new ColtEigenSystem(stateCount);

    }

    public EigenDecomposition getEigenDecomposition() {
        synchronized (this) {
            if (updateMatrix) {
                computeKroneckerSumsAndProducts();
            }
        }
        return eigenDecomposition;
    }

    public void getInfinitesimalMatrix(double[] out) {
        getEigenDecomposition(); // Updates rate matrix if necessary
        System.arraycopy(rateMatrix, 0, out, 0, stateCount * stateCount);
    }

    private void computeKroneckerSumsAndProducts() {

        int currentStateSize = stateSizes[0];
        double[] currentRate = new double[currentStateSize * currentStateSize];
        baseModels.get(0).getInfinitesimalMatrix(currentRate);
        currentRate = scaleForProductChain(currentRate, 0);

        for (int i = 1; i < numBaseModel; i++) {
            SubstitutionModel nextModel = baseModels.get(i);
            int nextStateSize = stateSizes[i];
            double[] nextRate = new double[nextStateSize * nextStateSize];
            nextModel.getInfinitesimalMatrix(nextRate);
            nextRate = scaleForProductChain(nextRate, i);
            currentRate = KroneckerOperation.sum(currentRate, currentStateSize, nextRate, nextStateSize);
            currentStateSize *= nextStateSize;
        }

        double dependence = dependenceParameter.getParameterValue(0);

        /* put dependence parameter into anti-diagonal elements*/
        for (int j = 0; j < stateCount; j++) {
            currentRate[(j + 1) * stateCount - (j + 1)] = dependence;
        }

        /* put dependence parameter into co-varying sites */

            /* take nucleotide for example, A, T, C, and G as 0, 1, 2, and 3 */
        String[] letA = new String[(int) Math.sqrt(stateCount)];
        String[] letB = new String[(int) Math.sqrt(stateCount)];
        String[] letAB = new String[stateCount];
        int count = 0;

        for (int letI = 0; letI < letA.length; letI++) {
            letA[letI] = String.valueOf(count);
            letB[letI] = String.valueOf(count);
            count++;
        }

            /* generate two-state array {(00), (01), (02), (03), (10), ..., (33)} */
        int letABi = 0;
        for (int letAi = 0; letAi < letA.length; letAi++) {
            for (int letBi = 0; letBi < letB.length; letBi++) {
                letAB[letABi] = letA[letAi] + " " + letB[letBi];
                letABi++;
            }
        }

            /* assign dependence parameter to co-varying sites */
        for (int j = 0; j < letAB.length; j++) {
            for (int k = 0; k < letAB.length; k++) {
                String[] Spl1 = letAB[j].split(" ");
                String[] Spl2 = letAB[k].split(" ");

                if (!Spl1[0].equals(Spl2[0]) && !Spl1[1].equals(Spl2[1])) {
                                                  /* Is co-varying site */
                    currentRate[j * stateCount + k] = dependence;
                }

            }
        }

        /* Rescale diagonal elements in currentRate */
        double diagonalSum = 0;
        for (int k = 0; k < stateCount; k++) {
            for (int l = k * stateCount; l < ((k + 1) * stateCount); l++) {
                if (l != (k * stateCount + k)) {
                    diagonalSum += currentRate[l];
                }
            }
                /* assign value to diagonal element */
            currentRate[k * stateCount + k] = -(diagonalSum);
            diagonalSum = 0;
        }

        rateMatrix = currentRate;

        // Call standard eigendecomposition on currentRate
        decompose();

    }

    private void decompose() {

        double[][] q = new double[stateCount][stateCount];
        for (int i = 0; i < stateCount; i++) {
            System.arraycopy(rateMatrix, i * stateCount, q[i], 0, stateCount);
        }

        eigenDecomposition = eigenSystemCES.decomposeMatrix(q);
        //System.err.println("eigenSystem null? " + (eigenSystem == null ? "yes" : "no"));
        //eigenDecomposition = eigenSystem.decomposeMatrix(q);

        updateMatrix = false;
    }

    private final Parameter dependenceParameter;

    private final EigenSystem eigenSystemCES; //= new ColtEigenSystem();

}
