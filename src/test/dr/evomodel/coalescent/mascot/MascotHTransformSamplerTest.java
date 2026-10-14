/*
 * MascotHTransformSamplerTest.java
 *
 * Copyright (c) 2002-2026 Alexei Drummond, Andrew Rambaut and Marc Suchard
 */

package test.dr.evomodel.coalescent.mascot;

import dr.evomodel.coalescent.mascot.MascotLikelihood;
import dr.evomodel.coalescent.mascot.MascotStochasticMap;
import dr.evomodel.coalescent.mascot.StochasticMapSampler;
import dr.math.MathUtils;
import junit.framework.TestCase;

/**
 * The h-transform sampler ({@link StochasticMapSampler#H_TRANSFORM}) targets the continuous frozen-environment
 * surrogate, so its draws are compared with Pi_h on a four-times finer grid (an O(h^2) proxy for the continuum).
 */
public class MascotHTransformSamplerTest extends TestCase {

    private static final int FINER = 4;

    public void testEveryHistoryIsCompleteAndConsistent() throws Exception {
        MathUtils.setSeed(201);
        for (MascotMappingFixtures.Spec spec : MascotMappingFixtures.oracleCases()) {
            MascotStochasticMapTest.checkStructure(hTransform(spec, 0.01), 200);
            MascotStochasticMapTest.checkStructure(hTransform(spec, 0.25), 200);
        }
        MascotStochasticMapTest.checkStructure(hTransform(MascotMappingFixtures.singleMobile(), 0.05), 200);
        MascotStochasticMapTest.checkStructure(hTransform(MascotMappingFixtures.highPressure(), 0.01), 50);
    }

    public void testC1MatchesContinuum() throws Exception {
        check(MascotMappingFixtures.c1(), 0.01, 1);
    }

    public void testC2MatchesContinuum() throws Exception {
        check(MascotMappingFixtures.c2(), 0.01, 2);
    }

    public void testC3AmbiguousTipMatchesContinuum() throws Exception {
        check(MascotMappingFixtures.c3(), 0.02, 3);
    }

    public void testC4EpochsMatchContinuum() throws Exception {
        check(MascotMappingFixtures.c4(), 0.01, 4);
    }

    public void testC5StrictClockMatchesContinuum() throws Exception {
        check(MascotMappingFixtures.c5(), 0.01, 5);
    }

    public void testC6BranchRatesMatchContinuum() throws Exception {
        check(MascotMappingFixtures.c6(), 0.02, 6);
    }

    public void testSwitchingSamplerRedraws() throws Exception {
        MathUtils.setSeed(7);
        MascotMappingFixtures.Built built = MascotMappingFixtures.c3().build(0.02);
        MascotLikelihood likelihood = built.likelihood;
        assertEquals(StochasticMapSampler.GRID, likelihood.getStochasticMapSampler());
        MascotStochasticMap grid = likelihood.getStochasticMap();
        assertEquals(StochasticMapSampler.GRID, grid.getSampler());
        likelihood.setStochasticMapSampler(StochasticMapSampler.H_TRANSFORM);
        MascotStochasticMap hTransform = likelihood.getStochasticMap();
        assertEquals(StochasticMapSampler.H_TRANSFORM, hTransform.getSampler());
        assertNotSame(grid, hTransform);
        assertSame(hTransform, likelihood.getStochasticMap());
        assertEquals(likelihood.getLogLikelihood(), hTransform.getLogLikelihood(), 0.0);
    }

    /** On stiff steps the h-transform is biased, so such states are drawn, exactly, by the grid sampler. */
    public void testStiffStepsFallBackToGrid() throws Exception {
        MascotMappingFixtures.Built stiff = hTransform(MascotMappingFixtures.c2(), 0.3);
        MascotStochasticMap map = stiff.likelihood.getStochasticMap();
        assertTrue(map.getMaxRateTimesLength() > 1.0);
        assertEquals(StochasticMapSampler.GRID, map.getSampler());
        assertEquals(StochasticMapSampler.H_TRANSFORM, stiff.likelihood.getStochasticMapSampler());
        MascotStochasticMapTest.checkAgainstOracle(stiff, 20000, 9, stiff.oracle());

        MascotMappingFixtures.Built mild = hTransform(MascotMappingFixtures.c2(), 0.01);
        assertEquals(StochasticMapSampler.H_TRANSFORM, mild.likelihood.getStochasticMap().getSampler());
    }

    private static void check(MascotMappingFixtures.Spec spec, double maxStep, long seed) throws Exception {
        MascotMappingFixtures.Built built = hTransform(spec, maxStep);
        MascotStochasticMapTest.checkAgainstOracle(built, 20000, seed, built.oracle(maxStep / FINER));
    }

    private static MascotMappingFixtures.Built hTransform(MascotMappingFixtures.Spec spec, double maxStep)
            throws Exception {
        MascotMappingFixtures.Built built = spec.build(maxStep);
        built.likelihood.setStochasticMapSampler(StochasticMapSampler.H_TRANSFORM);
        return built;
    }
}
