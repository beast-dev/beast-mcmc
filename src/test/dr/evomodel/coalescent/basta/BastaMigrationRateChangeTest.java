/*
 * BastaMigrationRateChangeTest.java
 *
 * Copyright (c) 2002-2026 Alexei Drummond, Andrew Rambaut and Marc Suchard
 */

package test.dr.evomodel.coalescent.basta;

import dr.evolution.alignment.Patterns;
import dr.evolution.datatype.GeneralDataType;
import dr.evolution.io.NewickImporter;
import dr.evomodel.branchratemodel.StrictClockBranchRates;
import dr.evomodel.coalescent.basta.BastaLikelihood;
import dr.evomodel.coalescent.basta.ConstantPopulationSizeModel;
import dr.evomodel.coalescent.basta.GenericBastaLikelihoodDelegate;
import dr.evomodel.substmodel.ComplexSubstitutionModel;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.SVSComplexSubstitutionModel;
import dr.evomodel.tree.DefaultTreeModel;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Parameter;
import junit.framework.TestCase;

/**
 * Changing a migration rate (or a BSSVS indicator) must invalidate BASTA's cached log likelihood, and
 * restoreModelState must bring back both the rates and the matching likelihood. "Full" below is the
 * likelihood after makeDirty(), i.e. recomputed from scratch.
 */
public class BastaMigrationRateChangeTest extends TestCase {

    private static final double TOL = 1e-9;

    private Parameter rates;
    private Parameter indicators;
    private Parameter popSizes;
    private BastaLikelihood likelihood;

    private void build(boolean svs) throws Exception {
        GeneralDataType dataType = new GeneralDataType(new String[]{"0", "1", "2"});
        TreeModel tree = new DefaultTreeModel("tree", new NewickImporter(
                "(((A:1.0,B:1.0):1.2,(C:1.5,D:1.5):0.7):0.8,(E:2.0,F:2.0):1.0);").importTree(null));

        Patterns patterns = new Patterns(dataType, tree);
        int[] tipStates = {0, 1, 2, 0, 1, 2};
        int[] pattern = new int[tree.getTaxonCount()];
        for (int i = 0; i < pattern.length; i++) {
            pattern[i] = tipStates[tree.getTaxonId(i).charAt(0) - 'A'];
        }
        patterns.addPattern(pattern);

        rates = new Parameter.Default("rates", new double[]{0.11, 0.23, 0.17, 0.31, 0.29, 0.37});
        FrequencyModel frequencies = new FrequencyModel(dataType, new double[]{1 / 3.0, 1 / 3.0, 1 / 3.0});
        ComplexSubstitutionModel substitutionModel;
        if (svs) {
            indicators = new Parameter.Default("indicators", new double[]{1, 1, 1, 1, 1, 1});
            substitutionModel = new SVSComplexSubstitutionModel("migration", dataType, frequencies, rates, indicators);
        } else {
            substitutionModel = new ComplexSubstitutionModel("migration", dataType, frequencies, rates);
        }
        substitutionModel.setNormalization(false);
        substitutionModel.setScaleRatesByFrequencies(false);

        popSizes = new Parameter.Default("popSizes", new double[]{0.18, 0.47, 0.34});
        likelihood = new BastaLikelihood("basta", tree, patterns, substitutionModel,
                new StrictClockBranchRates(new Parameter.Default("clock", 1.0)),
                new ConstantPopulationSizeModel("popSizeModel", popSizes, 3, -1), 1,
                new GenericBastaLikelihoodDelegate("basta", tree, 3, true), true, dataType, "states", false);
    }

    private double full() {
        likelihood.makeDirty();
        return likelihood.getLogLikelihood();
    }

    public void testEachRateChangeInvalidatesLikelihood() throws Exception {
        build(false);
        for (int i = 0; i < rates.getDimension(); i++) {
            double before = likelihood.getLogLikelihood();
            likelihood.storeModelState();
            rates.setParameterValue(i, rates.getParameterValue(i) * 1.5);

            double incremental = likelihood.getLogLikelihood();
            assertTrue("rate " + i + " left the likelihood unchanged", Math.abs(incremental - before) > TOL);
            assertEquals("rate " + i, full(), incremental, TOL);

            likelihood.acceptModelState();
        }
    }

    public void testStoreAccept() throws Exception {
        build(false);
        likelihood.getLogLikelihood();
        likelihood.storeModelState();
        rates.setParameterValue(2, 0.6);
        double proposed = likelihood.getLogLikelihood();
        likelihood.acceptModelState();

        assertEquals(0.6, rates.getParameterValue(2), 0.0);
        assertEquals(proposed, likelihood.getLogLikelihood(), 0.0);
        assertEquals(proposed, full(), TOL);
    }

    public void testStoreRestore() throws Exception {
        build(false);
        double original = likelihood.getLogLikelihood();
        likelihood.storeModelState();
        rates.setParameterValue(2, 0.6);
        likelihood.getLogLikelihood();
        likelihood.restoreModelState();

        assertEquals(0.17, rates.getParameterValue(2), 0.0);
        assertEquals(original, likelihood.getLogLikelihood(), 0.0);

        // A later move must see transition matrices that match the restored rates.
        likelihood.storeModelState();
        popSizes.setParameterValue(1, 0.6);
        double incremental = likelihood.getLogLikelihood();
        assertEquals(full(), incremental, TOL);
        likelihood.acceptModelState();

        popSizes.setParameterValue(1, 0.47);
        assertEquals(original, full(), TOL);
    }

    public void testIndicatorToggleInvalidatesLikelihood() throws Exception {
        build(true);
        double original = likelihood.getLogLikelihood();
        likelihood.storeModelState();
        indicators.setParameterValue(3, 0.0);

        double incremental = likelihood.getLogLikelihood();
        assertTrue(Math.abs(incremental - original) > TOL);
        assertEquals(full(), incremental, TOL);

        likelihood.restoreModelState();
        assertEquals(1.0, indicators.getParameterValue(3), 0.0);
        assertEquals(original, likelihood.getLogLikelihood(), 0.0);
        assertEquals(original, full(), TOL);
    }
}
