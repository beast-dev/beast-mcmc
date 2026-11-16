/*
 * ReflectiveHMCGeneralizedTest.java
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

package test.dr.inference.operators.hmc;

import dr.evolution.io.NewickImporter;
import dr.evomodel.tree.DefaultTreeModel;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.treedatalikelihood.discrete.NodeHeightBounds;
import dr.evomodel.treedatalikelihood.discrete.NodeHeightProxyParameter;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.GraphicalParameterBound;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.inference.model.PriorPreconditioningProvider;
import dr.inference.operators.AdaptationMode;
import dr.inference.operators.hmc.HamiltonianMonteCarloOperator;
import dr.inference.operators.hmc.MassPreconditioner;
import dr.inference.operators.hmc.MassPreconditioningOptions;
import dr.inference.operators.hmc.ReflectiveHamiltonianMonteCarloOperator;
import dr.inference.operators.hmc.ReflectiveHamiltonianMonteCarloOperator.ReflectionRule;
import dr.math.matrixAlgebra.WrappedVector;
import junit.framework.TestCase;

import java.lang.reflect.Field;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.util.Arrays;
import java.util.Random;

/**
 * Tests of reflective HMC's wall handling (plan 06, tests U1-U3, U5 and U8): Gamma on all walls hit at once against a
 * dense projection, exact ties, reversibility, identical interior drifts under both rules, tied starts, the zero-time
 * cap and maxEventsPerStep. Drifts are run through the operator's leap-frog engine with a zero gradient.
 */
public class ReflectiveHMCGeneralizedTest extends TestCase {

    // ------------------------------------------------------------------ set-up helpers

    static final class NoGrad implements GradientWrtParameterProvider {
        final Parameter p;

        NoGrad(Parameter p) { this.p = p; }

        public Likelihood getLikelihood() { return null; }

        public Parameter getParameter() { return p; }

        public int getDimension() { return p.getDimension(); }

        public double[] getGradientLogDensity() { return new double[p.getDimension()]; }
    }

    static MassPreconditioner diagonal(final double[] inverseMass) {
        return new MassPreconditioner.PriorPreconditioner(new PriorPreconditioningProvider() {
            public double getStandardDeviation(int index) { return Math.sqrt(inverseMass[index]); }

            public int getDimension() { return inverseMass.length; }
        }, null);
    }

    static ReflectiveHamiltonianMonteCarloOperator operator(Parameter p, GraphicalParameterBound bounds,
                                                            double[] inverseMass, ReflectionRule rule, long maxEvents) {
        MassPreconditioningOptions mpo = new MassPreconditioningOptions.Default(0, 0, 2, 0, false, null, null);
        HamiltonianMonteCarloOperator.Options options = new HamiltonianMonteCarloOperator.Options(0.1, 1, 0.0, mpo,
                0, 0.0, 10, 0.1, 0.8, HamiltonianMonteCarloOperator.InstabilityHandler.REJECT);
        return new ReflectiveHamiltonianMonteCarloOperator(AdaptationMode.ADAPTATION_OFF, 1.0, new NoGrad(p), p, null,
                null, options, diagonal(inverseMass), bounds, rule, maxEvents);
    }

    static ReflectiveHamiltonianMonteCarloOperator operator(Parameter p, GraphicalParameterBound bounds,
                                                            double[] inverseMass, ReflectionRule rule) {
        return operator(p, bounds, inverseMass, rule, ReflectiveHamiltonianMonteCarloOperator.DEFAULT_MAX_EVENTS_PER_STEP);
    }

    /** one drift of length step through the operator's engine; position and momentum are updated in place */
    static void drift(ReflectiveHamiltonianMonteCarloOperator op, double[] x, double[] p, double step)
            throws HamiltonianMonteCarloOperator.NumericInstabilityException {
        try {
            Field f = HamiltonianMonteCarloOperator.class.getDeclaredField("leapFrogEngine");
            f.setAccessible(true);
            Object engine = f.get(op);
            Method update = null;
            for (Class<?> c = engine.getClass(); c != null && update == null; c = c.getSuperclass()) {
                try {
                    update = c.getDeclaredMethod("updatePosition", double[].class, WrappedVector.class, double.class);
                } catch (NoSuchMethodException e) {
                    // look in the superclass
                }
            }
            update.setAccessible(true);
            update.invoke(engine, x, new WrappedVector.Raw(p), step);
        } catch (InvocationTargetException e) {
            if (e.getCause() instanceof HamiltonianMonteCarloOperator.NumericInstabilityException) {
                throw (HamiltonianMonteCarloOperator.NumericInstabilityException) e.getCause();
            }
            throw new RuntimeException(e.getCause());
        } catch (NoSuchFieldException | IllegalAccessException e) {
            throw new RuntimeException(e);
        }
    }

    /** a tree with all tips at 0 and the given newick branch lengths; heights through a rootNode="true" proxy */
    static final class Tree {
        final TreeModel model;
        final Parameter heights;
        final NodeHeightBounds bounds;

        Tree(String newick) throws Exception {
            model = new DefaultTreeModel("tree", new NewickImporter(newick).importTree(null), true, false);
            heights = new NodeHeightProxyParameter("heights", model, true);
            bounds = new NodeHeightBounds(heights, model);
        }
    }

    static String caterpillar(int k, double branch) {
        String s = "t0:" + branch;
        for (int i = 1; i <= k; i++) {
            s = "(" + s + ",t" + i + ":" + (i == 1 ? branch : branch * i) + "):" + branch;
        }
        return s.substring(0, s.lastIndexOf(':')) + ";";
    }

    static String randomTree(int tips, Random rng) {
        java.util.List<String> nodes = new java.util.ArrayList<>();
        java.util.List<Double> depth = new java.util.ArrayList<>();
        for (int i = 0; i < tips; i++) {
            nodes.add("t" + i);
            depth.add(0.0);
        }
        double h = 0.0;
        while (nodes.size() > 1) {
            int a = rng.nextInt(nodes.size());
            String na = nodes.remove(a);
            double da = depth.remove(a);
            int b = rng.nextInt(nodes.size());
            String nb = nodes.remove(b);
            double db = depth.remove(b);
            h += 0.1 + rng.nextDouble();
            nodes.add("(" + na + ":" + (h - da) + "," + nb + ":" + (h - db) + ")");
            depth.add(h);
        }
        return nodes.get(0) + ";";
    }

    static double[] inverseMasses(int dim, double logSd, Random rng) {
        double[] a = new double[dim];
        for (int i = 0; i < dim; i++) {
            a[i] = Math.exp(logSd * rng.nextGaussian());
        }
        return a;
    }

    static double[] momenta(double[] inverseMass, Random rng) {
        double[] p = new double[inverseMass.length];
        for (int i = 0; i < p.length; i++) {
            p[i] = rng.nextGaussian() / Math.sqrt(inverseMass[i]);
        }
        return p;
    }

    static double kinetic(double[] p, double[] inverseMass) {
        double k = 0.0;
        for (int i = 0; i < p.length; i++) {
            k += 0.5 * p[i] * p[i] * inverseMass[i];
        }
        return k;
    }

    // ------------------------------------------------------------------ a dense reference for Gamma

    /** walls of a tree: normals (rows) of the edge walls (child minus parent) and the floor walls (-e_i) */
    static final class Walls {
        final java.util.List<double[]> normals = new java.util.ArrayList<>();
        final java.util.List<double[]> offsets = new java.util.ArrayList<>();   // {b} with slack b - a.x

        Walls(Tree tree, double[] x) {
            final int dim = x.length;
            for (int i = 0; i < dim; i++) {
                int[] connected = tree.bounds.getConnectedParameterIndices(i);
                if (connected != null) {
                    for (int j : connected) {
                        if (j > i) {
                            double[] a = new double[dim];
                            int child = tree.bounds.isBelow(i, j) ? i : j;
                            int parent = child == i ? j : i;
                            a[child] = 1.0;
                            a[parent] = -1.0;
                            normals.add(a);
                            offsets.add(new double[]{0.0});
                        }
                    }
                }
                double lower = tree.bounds.getFixedLowerBound(i);
                if (lower > Double.NEGATIVE_INFINITY) {
                    double[] a = new double[dim];
                    a[i] = -1.0;
                    normals.add(a);
                    offsets.add(new double[]{-lower});
                }
            }
        }
    }

    /**
     * Resolves every wall at zero slack being hit, Gamma on all of them at once per pass, by projection in the
     * M^-1 metric (Gram-Schmidt on the sqrt(a)-scaled normals), until none is hit. Returns the number of passes.
     */
    static int denseTimeZeroResolution(Walls walls, double[] x, double[] p, double[] inverseMass, int maxPasses) {
        final int dim = x.length;
        for (int pass = 0; pass < maxPasses; pass++) {
            java.util.List<double[]> hit = new java.util.ArrayList<>();
            for (int w = 0; w < walls.normals.size(); w++) {
                double[] a = walls.normals.get(w);
                double slack = walls.offsets.get(w)[0];
                double rate = 0.0;
                for (int i = 0; i < dim; i++) {
                    slack -= a[i] * x[i];
                    rate += a[i] * inverseMass[i] * p[i];
                }
                if (slack <= 0 && rate > 0) {
                    hit.add(a);
                }
            }
            if (hit.isEmpty()) {
                return pass;
            }
            // u = M^{-1/2} p; reflect u across span{M^{-1/2} a_w}, then p = M^{1/2} u
            java.util.List<double[]> basis = new java.util.ArrayList<>();
            for (double[] a : hit) {
                double[] q = new double[dim];
                for (int i = 0; i < dim; i++) {
                    q[i] = a[i] * Math.sqrt(inverseMass[i]);
                }
                double norm0 = Math.sqrt(dot(q, q));
                for (int repeat = 0; repeat < 2; repeat++) {
                    for (double[] e : basis) {
                        double c = dot(q, e);
                        for (int i = 0; i < dim; i++) {
                            q[i] -= c * e[i];
                        }
                    }
                }
                double norm = Math.sqrt(dot(q, q));
                if (norm > 1e-10 * norm0) {
                    for (int i = 0; i < dim; i++) {
                        q[i] /= norm;
                    }
                    basis.add(q);
                }
            }
            double[] u = new double[dim];
            for (int i = 0; i < dim; i++) {
                u[i] = p[i] * Math.sqrt(inverseMass[i]);
            }
            double[] proj = new double[dim];
            for (double[] e : basis) {
                double c = dot(u, e);
                for (int i = 0; i < dim; i++) {
                    proj[i] += c * e[i];
                }
            }
            for (int i = 0; i < dim; i++) {
                p[i] = (u[i] - 2.0 * proj[i]) / Math.sqrt(inverseMass[i]);
            }
        }
        return -1;
    }

    static double dot(double[] a, double[] b) {
        double s = 0.0;
        for (int i = 0; i < a.length; i++) {
            s += a[i] * b[i];
        }
        return s;
    }

    // ------------------------------------------------------------------ tests

    /**
     * U1 and U2: from fully tied trees (every node and tip at 0), a drift too short to move resolves only the time-0
     * walls; the engine's momenta under the generalized rule match the dense reference, pass for pass.
     */
    public void testGammaMatchesDenseReferenceAtTiedStarts() throws Exception {
        Random rng = new Random(20261115);
        int checked = 0;
        for (int draw = 0; draw < 300; draw++) {
            int tips = 3 + rng.nextInt(18);
            String newick = draw % 3 == 0 ? caterpillar(tips - 1, 0.0) : randomTree(tips, rng)
                    .replaceAll(":[0-9.Ee-]+", ":0.0");
            Tree tree = new Tree(newick);
            double[] x = tree.heights.getParameterValues();
            final int dim = x.length;
            double[] a = inverseMasses(dim, draw % 2 == 0 ? 0.0 : 1.0, rng);
            double[] p = momenta(a, rng);

            double[] pDense = p.clone();
            int passes = denseTimeZeroResolution(new Walls(tree, x), x.clone(), pDense, a, 100000);
            assertTrue("dense reference ends", passes >= 0);

            ReflectiveHamiltonianMonteCarloOperator op = operator(tree.heights, tree.bounds, a,
                    ReflectionRule.GENERALIZED);
            double[] xe = x.clone();
            double[] pe = p.clone();
            drift(op, xe, pe, 1e-300);
            for (int i = 0; i < dim; i++) {
                assertEquals("momentum " + i + " in draw " + draw, pDense[i], pe[i],
                        1e-9 * (1.0 + Math.abs(pDense[i])));
            }
            assertEquals("energy in draw " + draw, kinetic(p, a), kinetic(pe, a), 1e-9 * (1.0 + kinetic(p, a)));
            assertEquals("events in draw " + draw, passes, op.getEventCount());
            assertTrue(tree.bounds.isWithinBounds(xe));
            ++checked;
        }
        assertEquals(300, checked);
    }

    /** U5: away from exact ties the two rules give bit-identical drifts. */
    public void testInteriorDriftsIdenticalUnderBothRules() throws Exception {
        Random rng = new Random(11);
        long events = 0;
        for (int draw = 0; draw < 400; draw++) {
            Tree tree = new Tree(randomTree(4 + rng.nextInt(30), rng));
            double[] x = tree.heights.getParameterValues();
            double[] a = inverseMasses(x.length, 0.5 * (draw % 4), rng);
            double[] p = momenta(a, rng);
            double step = 0.5 + 3.0 * rng.nextDouble();

            double[] xg = x.clone(), pg = p.clone(), xs = x.clone(), ps = p.clone();
            ReflectiveHamiltonianMonteCarloOperator g = operator(tree.heights, tree.bounds, a, ReflectionRule.GENERALIZED);
            drift(g, xg, pg, step);
            ReflectiveHamiltonianMonteCarloOperator s = operator(tree.heights, tree.bounds, a, ReflectionRule.SEQUENTIAL);
            drift(s, xs, ps, step);
            assertTrue("positions in draw " + draw, Arrays.equals(xg, xs));
            assertTrue("momenta in draw " + draw, Arrays.equals(pg, ps));
            assertEquals(g.getEventCount(), s.getEventCount());
            assertEquals(0, g.getMultiWallEventCount());
            events += g.getEventCount();
        }
        assertTrue("drifts with events", events > 400);
    }

    /** U3: forward, flip the momentum, forward again returns to the start, from interior and from tied starts. */
    public void testReversibility() throws Exception {
        Random rng = new Random(7);
        for (int draw = 0; draw < 200; draw++) {
            boolean tied = draw % 2 == 1;
            String newick = randomTree(4 + rng.nextInt(16), rng);
            if (tied) {
                newick = newick.replaceAll(":[0-9.Ee-]+", ":0.0");
            }
            Tree tree = new Tree(newick);
            double[] x0 = tree.heights.getParameterValues();
            double[] a = inverseMasses(x0.length, 0.5, rng);
            double[] p0 = momenta(a, rng);
            double step = tied ? 0.3 : 1.0 + 2.0 * rng.nextDouble();
            ReflectiveHamiltonianMonteCarloOperator op = operator(tree.heights, tree.bounds, a,
                    ReflectionRule.GENERALIZED);
            double[] x = x0.clone(), p = p0.clone();
            drift(op, x, p, step);
            for (int i = 0; i < p.length; i++) {
                p[i] = -p[i];
            }
            drift(op, x, p, step);
            for (int i = 0; i < x.length; i++) {
                assertEquals("position " + i + " in draw " + draw, x0[i], x[i], 1e-9 * (1.0 + Math.abs(x0[i])));
                if (!tied) {
                    assertEquals("momentum " + i + " in draw " + draw, -p0[i], p[i], 1e-9 * (1.0 + Math.abs(p0[i])));
                }
            }
        }
    }

    /** U2 for boxes: a start in a corner, moving out of both walls, flips both momenta in one generalized event. */
    public void testBoxCornerTie() throws Exception {
        for (ReflectionRule rule : ReflectionRule.values()) {
            Parameter q = new Parameter.Default(new double[]{0.0, 1.0});
            q.addBounds(new Parameter.DefaultBounds(1.0, 0.0, 2));
            GraphicalParameterBound box = new GraphicalParameterBound.FixedBound(q);
            ReflectiveHamiltonianMonteCarloOperator op = operator(q, box, new double[]{1.0, 1.0}, rule);
            double[] x = {0.0, 1.0};
            double[] p = {-1.0, 1.0};
            drift(op, x, p, 0.25);
            assertEquals(1.0, p[0]);
            assertEquals(-1.0, p[1]);
            assertEquals(0.25, x[0], 1e-15);
            assertEquals(0.75, x[1], 1e-15);
            assertEquals(rule == ReflectionRule.GENERALIZED ? 1 : 2, op.getEventCount());
            assertEquals(rule == ReflectionRule.GENERALIZED ? 1 : 0, op.getMultiWallEventCount());
        }
    }

    /** U8: a fully tied caterpillar of 300 nodes with equal masses leaves its walls under both rules. */
    public void testTiedCaterpillarLeavesUnderBothRules() throws Exception {
        Random rng = new Random(3);
        for (ReflectionRule rule : ReflectionRule.values()) {
            for (int draw = 0; draw < 3; draw++) {
                Tree tree = new Tree(caterpillar(300, 0.0));
                double[] x = tree.heights.getParameterValues();
                double[] a = new double[x.length];
                Arrays.fill(a, 1.0);
                double[] p = momenta(a, rng);
                ReflectiveHamiltonianMonteCarloOperator op = operator(tree.heights, tree.bounds, a, rule);
                drift(op, x, p, 0.1);
                assertTrue(tree.bounds.isWithinBounds(x));
                assertEquals(0, op.getZeroTimeCapCount());
                if (rule == ReflectionRule.GENERALIZED) {
                    assertTrue("about 2k events: " + op.getEventCount(), op.getEventCount() < 2000);
                } else {
                    assertTrue("about k^2/2 events: " + op.getEventCount(), op.getEventCount() > 30000);
                }
            }
        }
    }

    /** An endless sequence of time-0 events (a degenerate box) stops after exactly 10^6 + 1 of them. */
    public void testZeroTimeCap() throws Exception {
        for (ReflectionRule rule : ReflectionRule.values()) {
            Parameter q = new Parameter.Default(new double[]{1.0});
            q.addBounds(new Parameter.DefaultBounds(1.0, 1.0, 1));
            ReflectiveHamiltonianMonteCarloOperator op = operator(q, new GraphicalParameterBound.FixedBound(q),
                    new double[]{1.0}, rule);
            try {
                drift(op, new double[]{1.0}, new double[]{1.0}, 0.1);
                fail("an endless sequence must be stopped");
            } catch (HamiltonianMonteCarloOperator.NumericInstabilityException e) {
                assertEquals(1, op.getZeroTimeCapCount());
                assertEquals(1000001L, op.getEventCount());
            }
        }
    }

    /** maxEventsPerStep rejects a drift with more events than it allows. */
    public void testMaxEventsPerStep() throws Exception {
        Tree tree = new Tree(caterpillar(30, 0.0));
        double[] x = tree.heights.getParameterValues();
        double[] a = new double[x.length];
        Arrays.fill(a, 1.0);
        ReflectiveHamiltonianMonteCarloOperator op = operator(tree.heights, tree.bounds, a, ReflectionRule.SEQUENTIAL,
                10);
        try {
            drift(op, x, momenta(a, new Random(5)), 0.1);
            fail("about 450 events cannot pass a limit of 10");
        } catch (HamiltonianMonteCarloOperator.NumericInstabilityException e) {
            assertEquals(1, op.getMaxEventsRejectionCount());
            assertEquals(11, op.getEventCount());
        }
    }

    public void testReflectionRuleNames() {
        assertEquals(ReflectionRule.GENERALIZED, ReflectionRule.parse("generalized"));
        assertEquals(ReflectionRule.SEQUENTIAL, ReflectionRule.parse("Sequential"));
        try {
            ReflectionRule.parse("corner");
            fail();
        } catch (IllegalArgumentException e) {
            // expected
        }
    }
}
