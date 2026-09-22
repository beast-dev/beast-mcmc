/*
 * AgeDependentBirthDeathSimulatorTest.java
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

package test.dr.evomodel.speciation;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeUtils;
import dr.evomodel.speciation.agedependent.agehazard.AgeHazard;
import dr.evomodel.speciation.agedependent.agehazard.ExpAgeHazard;
import dr.evomodel.speciation.agedependent.agehazard.LinExpAgeHazard;
import dr.evomodel.speciation.agedependent.agehazard.OffsetLinExpAgeHazard;
import dr.evomodel.speciation.agedependent.simulation.AgeDependentBirthDeathSimulator;
import dr.inference.model.Parameter;
import dr.math.MathUtils;
import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

/**
 * Tests for the age-dependent birth-death tree simulator.
 *
 * @author Frederik M. Andersen
 */
public class AgeDependentBirthDeathSimulatorTest extends TestCase {

    private static final double ORIGIN = 8.0;
    private static final int MAX_LINEAGES = 100000;
    private static final int MAX_ATTEMPTS = 2000;
    private static final int MIN_TIPS = 10;

    public AgeDependentBirthDeathSimulatorTest(String name) {
        super(name);
    }

    private static Parameter p(double... values) {
        return new Parameter.Default(values);
    }

    /**
     * One epoch, no serial sampling, complete extant sampling: the ultrametric process.
     */
    private static Tree simulateUltrametric(long seed, AgeHazard birth, AgeHazard death) {
        MathUtils.setSeed(seed);
        return new AgeDependentBirthDeathSimulator(
                new double[]{1.5}, new double[]{0.8}, new double[]{0.0}, 1.0,
                birth, death, new double[]{}, ORIGIN, true, MAX_LINEAGES
        ).simulate(MIN_TIPS, MAX_ATTEMPTS);
    }

    /**
     * Two epochs, serial sampling at psi and extant sampling at rho.
     */
    private static Tree simulateSerial(long seed, AgeHazard birth, AgeHazard death,
                                       double psi, double rho) {
        MathUtils.setSeed(seed);
        return new AgeDependentBirthDeathSimulator(
                new double[]{1.2, 1.8}, new double[]{0.6, 1.0}, new double[]{psi}, rho,
                birth, death, new double[]{3.0}, ORIGIN, true, MAX_LINEAGES
        ).simulate(MIN_TIPS, MAX_ATTEMPTS);
    }

    private static LinExpAgeHazard linExp(double r, double gamma) {
        return new LinExpAgeHazard(p(r), p(gamma));
    }

    /**
     * Regression test: the simulator assigns node heights only, so the tree it returns must
     * still report branch lengths derived from those heights. Constructing the FlexibleTree
     * as "lengths already known" silently left every branch at the unset default of 0, which
     * no height-based check would notice.
     */
    public void testBranchLengthsAreAssigned() {
        Tree tree = simulateUltrametric(42, linExp(0.4, 0.25), linExp(0.3, 0.15));

        double sumReported = 0.0;
        double sumImplied = 0.0;
        int zeroLength = 0;
        int branches = 0;

        for (int i = 0; i < tree.getNodeCount(); i++) {
            NodeRef node = tree.getNode(i);
            if (tree.isRoot(node)) continue;
            branches++;
            double reported = tree.getBranchLength(node);
            double implied = tree.getNodeHeight(tree.getParent(node)) - tree.getNodeHeight(node);
            if (reported == 0.0) zeroLength++;
            sumReported += reported;
            sumImplied += implied;
        }

        assertTrue("simulated tree should have branches", branches > 0);
        assertEquals("every branch length should be assigned", 0, zeroLength);
        assertTrue("total tree length should be positive", sumReported > 0.0);
        assertEquals("branch lengths must agree with node heights",
                sumImplied, sumReported, 1e-10);
    }

    /**
     * With psi = 0 and rho = 1 the process is ultrametric: every tip sits at the present.
     */
    public void testUltrametricModeHasContemporaneousTips() {
        Tree tree = simulateUltrametric(42, linExp(0.4, 0.25), linExp(0.3, 0.15));

        for (int i = 0; i < tree.getExternalNodeCount(); i++) {
            assertEquals("tip should be at the present",
                    0.0, tree.getNodeHeight(tree.getExternalNode(i)), 1e-12);
        }
    }

    /**
     * With psi > 0 the tree carries sampled-through-time tips above the present.
     */
    public void testSerialSamplingProducesTipsAboveThePresent() {
        Tree tree = simulateSerial(42, linExp(0.4, 0.25), linExp(0.3, 0.15), 0.3, 0.5);

        int above = 0;
        for (int i = 0; i < tree.getExternalNodeCount(); i++) {
            if (tree.getNodeHeight(tree.getExternalNode(i)) > 1e-12) above++;
        }
        assertTrue("serial sampling should yield tips above the present", above > 0);
        assertTrue("serial tips cannot outnumber all tips",
                above <= tree.getExternalNodeCount());
    }

    /**
     * All three hazard shapes reduce to h(a) = 1 for these parameters, so driven by the same
     * seed they must produce the very same tree. Catches a hazard wired in inconsistently.
     */
    public void testFlatHazardsAgree() {
        Tree fromLinExp = simulateUltrametric(99, linExp(0.0, 0.0), linExp(0.0, 0.0));
        Tree fromExp = simulateUltrametric(99,
                new ExpAgeHazard(p(0.0)), new ExpAgeHazard(p(0.0)));
        Tree fromOffset = simulateUltrametric(99,
                new OffsetLinExpAgeHazard(p(0.0), p(0.0), p(0.0)),
                new OffsetLinExpAgeHazard(p(0.0), p(0.0), p(0.0)));

        String linExpNewick = TreeUtils.newick(fromLinExp);
        assertEquals("expAgeHazard(0) should match linExpAgeHazard(0, 0)",
                linExpNewick, TreeUtils.newick(fromExp));
        assertEquals("offsetLinExpAgeHazard(0, 0, 0) should match linExpAgeHazard(0, 0)",
                linExpNewick, TreeUtils.newick(fromOffset));
    }

    /**
     * A hazard given the same values in every epoch must behave exactly like the shared-shape
     * hazard covering one epoch.
     */
    public void testPerEpochHazardMatchesSharedHazard() {
        Tree shared = simulateSerial(7, linExp(0.4, 0.25), linExp(0.3, 0.15), 0.3, 0.5);
        Tree perEpoch = simulateSerial(7,
                new LinExpAgeHazard(p(0.4, 0.4), p(0.25, 0.25)),
                new LinExpAgeHazard(p(0.3, 0.3), p(0.15, 0.15)), 0.3, 0.5);

        assertEquals("per-epoch hazard with identical values should match the shared hazard",
                TreeUtils.newick(shared), TreeUtils.newick(perEpoch));
    }

    /**
     * A hazard covering neither 1 nor numEpochs epochs is a specification error.
     */
    public void testEpochCountMismatchIsRejected() {
        try {
            simulateSerial(1,
                    new LinExpAgeHazard(p(0.4, 0.5, 0.6), p(0.25, 0.25, 0.25)),
                    linExp(0.3, 0.15), 0.3, 0.5);
            fail("a hazard covering 3 epochs should be rejected when there are 2");
        } catch (IllegalArgumentException expected) {
            assertTrue("message should name the offending hazard",
                    expected.getMessage().contains("birthHazard"));
        }
    }

    /**
     * The non-linExp hazards must be simulable at all -- they could not be before the
     * simulator took an AgeHazard.
     */
    public void testNonLinExpHazardsSimulate() {
        Tree fromExp = simulateUltrametric(3, new ExpAgeHazard(p(0.25)), new ExpAgeHazard(p(0.15)));
        assertTrue("expAgeHazard should produce a tree",
                fromExp.getExternalNodeCount() >= MIN_TIPS);

        Tree mixed = simulateUltrametric(3, linExp(0.4, 0.25), new ExpAgeHazard(p(0.15)));
        assertTrue("mixed hazard families should produce a tree",
                mixed.getExternalNodeCount() >= MIN_TIPS);
    }

    public static Test suite() {
        return new TestSuite(AgeDependentBirthDeathSimulatorTest.class);
    }
}
