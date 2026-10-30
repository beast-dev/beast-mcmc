/*
 * BastaMigrationRateGradientTest.java
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

package test.dr.evomodel.coalescent;

import dr.evolution.alignment.SimpleSiteList;
import dr.evolution.datatype.GeneralDataType;
import dr.evolution.io.NewickImporter;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxa;
import dr.evolution.util.Taxon;
import dr.evomodel.branchratemodel.StrictClockBranchRates;
import dr.evomodel.coalescent.basta.AdjointGenericBastaLikelihoodDelegate;
import dr.evomodel.coalescent.basta.BastaLikelihood;
import dr.evomodel.coalescent.basta.BastaLikelihoodDelegate;
import dr.evomodel.coalescent.basta.ConstantPopulationSizeModel;
import dr.evomodel.coalescent.basta.StructuredCoalescentLikelihoodGradient;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.SVSComplexSubstitutionModel;
import dr.evomodel.tree.DefaultTreeModel;
import dr.inference.model.Parameter;
import junit.framework.TestCase;

import java.util.Arrays;

/**
 * Migration-rate gradients of BASTA under BSSVS, checked against finite differences.
 */
public class BastaMigrationRateGradientTest extends TestCase {

    private static final String TREE = "((((A:1.0,B:0.6):0.8,C:2.1):0.5,(D:0.7,E:1.2):2.0):0.6,F:3.5);";
    private static final String[] TAXA = {"A", "B", "C", "D", "E", "F"};
    private static final int[] DEMES = {0, 1, 2, 0, 1, 2};

    /**
     * A rate whose BSSVS indicator is 0 does not enter the likelihood, so its gradient must be exactly 0. Index 3
     * (1 to 0) is where the rates-parameter order and the row-major migration order differ.
     */
    public void testBssvsGradientsMatchFiniteDifferences() throws Exception {
        for (boolean normalize : new boolean[]{false, true}) {
            String label = normalize ? "normalized" : "not normalized";

            Parameter rates = new Parameter.Default("location.rates", new double[]{0.8, 1.3, 0.5, 1.7, 0.9, 1.1});
            Parameter indicators = new Parameter.Default("location.indicators",
                    new double[]{1.0, 0.0, 1.0, 0.0, 1.0, 1.0});
            BastaLikelihood likelihood = basta(rates, indicators, normalize);

            double[] gradient = new StructuredCoalescentLikelihoodGradient(likelihood,
                    likelihood.getSubstitutionModel(),
                    StructuredCoalescentLikelihoodGradient.WrtParameter.MIGRATION_RATE).getGradientLogDensity();

            assertEquals(label, 0.0, gradient[1], 0.0);
            assertEquals(label, 0.0, gradient[3], 0.0);

            double[] numeric = finiteDifferences(likelihood, rates);
            for (int i = 0; i < gradient.length; ++i) {
                assertEquals(label + " rate " + i + ": " + Arrays.toString(gradient) + " vs " +
                        Arrays.toString(numeric), numeric[i], gradient[i], 1.0e-5 * Math.max(1.0, Math.abs(numeric[i])));
            }
        }
    }

    private static BastaLikelihood basta(Parameter rates, Parameter indicators, boolean normalize) throws Exception {
        GeneralDataType dataType = new GeneralDataType(new String[]{"d0", "d1", "d2"});

        Taxa taxa = new Taxa();
        for (String id : TAXA) {
            taxa.addTaxon(new Taxon(id));
        }
        SimpleSiteList patterns = new SimpleSiteList(dataType, taxa);
        patterns.addPattern(DEMES);

        FrequencyModel frequencies = new FrequencyModel(dataType, new double[]{0.2, 0.3, 0.5});
        SVSComplexSubstitutionModel migrationModel =
                new SVSComplexSubstitutionModel("location.model", dataType, frequencies, rates, indicators);
        migrationModel.setNormalization(normalize);

        Tree tree = new DefaultTreeModel(new NewickImporter(TREE).importTree(null));
        BastaLikelihoodDelegate delegate =
                new AdjointGenericBastaLikelihoodDelegate("basta", tree, dataType.getStateCount(), true, false);

        Parameter popSizes = new Parameter.Default("popSizes", new double[]{1.5, 0.7, 2.2});

        return new BastaLikelihood("basta", tree, patterns, migrationModel,
                new StrictClockBranchRates(new Parameter.Default(1.0)),
                new ConstantPopulationSizeModel("popSizeModel", popSizes, dataType.getStateCount(), -1), 2,
                delegate, true, dataType, "states", false);
    }

    private static double[] finiteDifferences(BastaLikelihood likelihood, Parameter rates) {
        double[] gradient = new double[rates.getDimension()];
        for (int i = 0; i < gradient.length; ++i) {
            double x = rates.getParameterValue(i);
            double h = 1.0e-6 * Math.max(1.0, Math.abs(x));
            rates.setParameterValue(i, x + h);
            double plus = likelihood.getLogLikelihood();
            rates.setParameterValue(i, x - h);
            double minus = likelihood.getLogLikelihood();
            rates.setParameterValue(i, x);
            gradient[i] = (plus - minus) / (2.0 * h);
        }
        likelihood.getLogLikelihood();
        return gradient;
    }
}
