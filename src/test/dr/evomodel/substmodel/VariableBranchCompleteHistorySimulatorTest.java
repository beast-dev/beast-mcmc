/*
 * VariableBranchCompleteHistorySimulatorTest.java
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

package test.dr.evomodel.substmodel;

import dr.evolution.tree.TreeUtils;
import dr.evomodel.siteratemodel.GammaSiteRateModel;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.nucleotide.HKY;
import dr.evolution.datatype.Nucleotides;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.branchratemodel.DefaultBranchRateModel;
import dr.inference.markovjumps.MarkovJumpsCore;
import dr.inference.model.Parameter;

/**
 * @author Marc A. Suchard
 */
public class VariableBranchCompleteHistorySimulatorTest extends CompleteHistorySimulatorTest {

    public int N = 5000;

    public void setUp() throws Exception {
        super.setUp();
    }

    public void testHKYVariableSimulation() {

        System.out.println("Starting HKY variable branch simulation");
        Parameter kappa = new Parameter.Default(1, 2.0);
        double[] pi = {0.45, 0.05, 0.25, 0.25};
        Parameter freqs = new Parameter.Default(pi);
        FrequencyModel f = new FrequencyModel(Nucleotides.INSTANCE, freqs);
        HKY hky = new HKY(kappa, f);
        int stateCount = hky.getDataType().getStateCount();

        Parameter mu = new Parameter.Default(1, 0.5);
        Parameter alpha = new Parameter.Default(1, 0.5);
        GammaSiteRateModel siteModel = new GammaSiteRateModel("gammaModel", mu, 1.0, alpha, 4, null);
        siteModel.setSubstitutionModel(hky);
        BranchRateModel branchRateModel = new DefaultBranchRateModel();

        double analyticResult = TreeUtils.getTreeLength(tree, tree.getRoot()) * mu.getParameterValue(0);
        int nSites = 200;

        double[] register1 = new double[stateCount * stateCount];
        double[] register2 = new double[stateCount * stateCount];

        MarkovJumpsCore.fillRegistrationMatrix(register1, stateCount); // Count all jumps

        // Move some jumps from 1 to 2
        register1[1 * stateCount + 2] = 0;
        register2[1 * stateCount + 2] = 1;

        register1[1 * stateCount + 3] = 0;
        register2[1 * stateCount + 3] = 1;

        register1[2 * stateCount + 3] = 0;
        register2[2 * stateCount + 3] = 1;

        double[] branchValues = { 10.0, 10.0, 10.0, 10.0, 10.0 };

        Parameter branchValuesParam = new Parameter.Default(branchValues);

        runSimulation(N, tree, siteModel, branchRateModel, nSites,
                new double[][] {register1, register2}, analyticResult, kappa, branchValuesParam);
    }

    public void testCodonSimulation() {
        // Do nothing
    }
}
