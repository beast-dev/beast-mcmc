/*
 * MultiTypeBirthDeathSerialSamplingModel.java
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

package dr.evomodel.speciation;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxon;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Set;

/**
 * A PROTOTYPE episodic multi-type birth-death-sampling model, generalizing
 * {@link NewBirthDeathSerialSamplingModel} (Stadler 2010; Gavryushkina et al. 2014) from
 * K=1 type to K&gt;=1 types, following the multi-type birth-death (MTBD) parameterization
 * of Stadler &amp; Bonhoeffer (2013) and Kuehnert et al. (2016, "BDMM").
 * <p>
 * Unlike {@code NewBirthDeathSerialSamplingModel}, which for K=1 exploits a closed-form
 * solution of the (scalar) Riccati equation for the extinction probability and plugs
 * into the highly-optimized interval-processing machinery of
 * {@link EfficientSpeciationLikelihood} / {@link SpeciationModelGradientProvider} (which
 * tracks only the *number* of contemporaneous lineages, not their individual states), a
 * genuine multi-type model must track a per-lineage, per-type probability vector, which
 * is fundamentally incompatible with that lineage-count abstraction. This prototype
 * therefore implements the classic, direct {@link SpeciationModel#calculateTreeLogLikelihood}
 * entry point instead (paired with the plain {@link SpeciationLikelihood} wrapper), and
 * solves the underlying nonlinear (extinction probability) and linear (edge/lineage
 * probability) ODE systems with a simple fixed-step 4th order Runge-Kutta integrator.
 * <p>
 * This is deliberately a numerically-exact-given-fine-enough-steps baseline: it is the
 * natural reference against which the fast MASCOT-inspired approximations discussed
 * separately (matrix-exponential epoch caching, mean-field/linearized extinction-probability
 * surrogates, weak-migration perturbation expansion around the exact K=1 solution, etc.)
 * should be validated, and against which it should be checked that this class exactly
 * reproduces {@code NewBirthDeathSerialSamplingModel} in the K=1 special case.
 * <p>
 * Modeling conventions and known limitations of this prototype (see also the design note
 * "multitype-BDS-fast-approximations.md" in notes/):
 * <ul>
 *     <li>Birth events follow the common MTBD convention: a birth of a type-<i>i</i> parent
 *     produces two lineages of types <i>i</i> and <i>j</i> at rate birthRate<sub>ij</sub>
 *     (i=j is an ordinary within-type birth). There is no support for births producing two
 *     daughters of two different, non-parental types.</li>
 *     <li>Trees must be strictly bifurcating; sampled-ancestor (degree-2) nodes are not
 *     supported (an exception is thrown if encountered).</li>
 *     <li>Tip types are read from a {@link Taxon} attribute (name configurable via
 *     {@code typeAttributeName}, default {@code "state"}), matched (as a string) against
 *     the supplied {@code typeLabels}, or given directly as an integer index.</li>
 *     <li>{@code conditionOnSurvival} divides by &Sigma;<sub>i</sub>(1-E<sub>i</sub>(origin)),
 *     i.e. it uses the same implicit "uniform sum over unknown origin type" convention as is
 *     used for the (equally unnormalized) numerator &Sigma;<sub>i</sub> p<sub>i</sub>(origin);
 *     no separate origin-type prior is currently exposed. A consequence worth noting: since
 *     every declared type is given equal weight-1 in this sum regardless of whether it can
 *     actually generate the observed tree, declaring an extra type that the data never visits
 *     still dilutes the conditioning denominator (it has its own nonzero survival probability)
 *     and so a K-type model with K-1 "unreachable" types is *not* numerically identical to the
 *     true K=1 model with conditionOnSurvival=true, even though their unconditioned likelihoods
 *     agree exactly -- see {@code MultiTypeBirthDeathLikelihoodTest} for a worked example. This
 *     is arguably the least satisfying modeling choice in this prototype and a good candidate
 *     for revisiting (e.g. exposing an explicit origin-type prior) before this goes beyond
 *     prototype status.</li>
 *     <li>No analytic gradient is provided yet (see the design note, section 3.9, for the
 *     planned adjoint-based approach); this class is a likelihood-only prototype.</li>
 * </ul>
 *
 * @author Marc A. Suchard
 * @author Claude
 */
public class MultiTypeBirthDeathSerialSamplingModel extends SpeciationModel {

    public static final String DEFAULT_TYPE_ATTRIBUTE = "state";

    private final int numTypes;
    private final String[] typeLabels;
    private final String typeAttributeName;

    // birthRate_{ij}: parent type i gives birth to lineages of type i and type j.
    // Size numTypes*numTypes (constant across epochs) or numTypes*numTypes*numIntervals
    // (row-major within each epoch block: index = (model*K + i)*K + j).
    private final Parameter birthRate;

    // migrationRate_{ij} (i != j): rate at which a lineage of type i becomes type j.
    // Same size convention as birthRate; diagonal entries are ignored.
    private final Parameter migrationRate;

    // deathRate_i, serialSamplingRate_i, treatmentProbability_i, samplingProbability_i:
    // size numTypes (constant across epochs) or numTypes*numIntervals.
    private final Parameter deathRate;
    private final Parameter serialSamplingRate;
    private final Parameter treatmentProbability;
    private final Parameter samplingProbability;

    private final Parameter originTime;

    private final boolean conditionOnSurvival;

    private final int numIntervals;
    private final double gridEnd;
    private final double[] modelStartTimes;

    // RK4 steps per unit time; controls the fineness of both the E(t) grid and the
    // per-branch integration of the linear lineage-probability ODE.
    private final int stepsPerUnitTime;

    private static final double TIME_TOL = 1e-8;

    // Cache: dense grid of the (tree-independent) extinction probability trajectory E(t),
    // solved once per parameter update and re-used (via linear interpolation) for every
    // lineage on the tree -- see design note section 3.1/3.3.
    private double[] eGridTimes;
    private double[][] eGridValues;
    private boolean eGridKnown = false;

    public MultiTypeBirthDeathSerialSamplingModel(
            String modelName,
            int numTypes,
            String[] typeLabels,
            Parameter birthRate,
            Parameter migrationRate,
            Parameter deathRate,
            Parameter serialSamplingRate,
            Parameter treatmentProbability,
            Parameter samplingProbability,
            Parameter originTime,
            boolean conditionOnSurvival,
            int numIntervals,
            double gridEnd,
            String typeAttributeName,
            int stepsPerUnitTime,
            Type units) {

        super(modelName, units);

        if (numTypes < 1) {
            throw new IllegalArgumentException("numTypes must be >= 1");
        }
        if (typeLabels.length != numTypes) {
            throw new IllegalArgumentException("typeLabels must have length numTypes");
        }

        this.numTypes = numTypes;
        this.typeLabels = typeLabels;
        this.typeAttributeName = (typeAttributeName != null) ? typeAttributeName : DEFAULT_TYPE_ATTRIBUTE;
        this.numIntervals = numIntervals;
        this.gridEnd = gridEnd;
        this.stepsPerUnitTime = stepsPerUnitTime;
        this.conditionOnSurvival = conditionOnSurvival;

        final int K = numTypes;
        final int KK = K * K;

        checkSize("birthRate", birthRate, KK, KK * numIntervals);
        checkSize("migrationRate", migrationRate, KK, KK * numIntervals);
        checkSize("deathRate", deathRate, K, K * numIntervals);
        checkSize("serialSamplingRate", serialSamplingRate, K, K * numIntervals);
        checkSize("treatmentProbability", treatmentProbability, K, K * numIntervals);
        checkSize("samplingProbability", samplingProbability, K, K * numIntervals);
        if (originTime.getDimension() != 1) {
            throw new IllegalArgumentException("originTime must have dimension 1");
        }

        this.birthRate = birthRate;
        addVariable(birthRate);
        this.migrationRate = migrationRate;
        addVariable(migrationRate);
        this.deathRate = deathRate;
        addVariable(deathRate);
        this.serialSamplingRate = serialSamplingRate;
        addVariable(serialSamplingRate);
        this.treatmentProbability = treatmentProbability;
        addVariable(treatmentProbability);
        this.samplingProbability = samplingProbability;
        addVariable(samplingProbability);
        this.originTime = originTime;
        addVariable(originTime);

        this.modelStartTimes = new double[numIntervals];
        for (int i = 1; i < numIntervals; ++i) {
            modelStartTimes[i] = i * (gridEnd / numIntervals);
        }
    }

    private static void checkSize(String name, Parameter p, int sizeConstant, int sizeEpisodic) {
        int size = p.getDimension();
        if (size != sizeConstant && size != sizeEpisodic) {
            throw new IllegalArgumentException("Parameter '" + name + "' has dimension " + size +
                    ", expected " + sizeConstant + " (constant across epochs) or " + sizeEpisodic +
                    " (one block of " + sizeConstant + " per epoch)");
        }
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        eGridKnown = false;
    }

    // ----------------------------------------------------------------------------------
    // Rate look-ups (handle both "constant across epochs" and "one block per epoch" sizes)
    // ----------------------------------------------------------------------------------

    private double rateK(Parameter p, int model, int i) {
        return (p.getDimension() == numTypes) ? p.getParameterValue(i) : p.getParameterValue(model * numTypes + i);
    }

    private double rateKK(Parameter p, int model, int i, int j) {
        int KK = numTypes * numTypes;
        return (p.getDimension() == KK) ? p.getParameterValue(i * numTypes + j)
                : p.getParameterValue(model * KK + i * numTypes + j);
    }

    private int epochIndexForTime(double t) {
        int m = 0;
        while (m < numIntervals - 1 && t >= modelStartTimes[m + 1]) {
            ++m;
        }
        return m;
    }

    // ----------------------------------------------------------------------------------
    // ODE right-hand sides
    // ----------------------------------------------------------------------------------

    /** Nonlinear (matrix Riccati) ODE for the tree-independent extinction probabilities E(t). */
    private double[] derivativeE(double[] E, int model) {
        int K = numTypes;
        double[] out = new double[K];
        for (int i = 0; i < K; ++i) {
            double mu = rateK(deathRate, model, i);
            double psi = rateK(serialSamplingRate, model, i);

            double birthSum = 0.0;
            double quadSum = 0.0; // sum_j birthRate_ij * E_j
            for (int j = 0; j < K; ++j) {
                double b = rateKK(birthRate, model, i, j);
                birthSum += b;
                quadSum += b * E[j];
            }

            double migOutSum = 0.0;
            double migInSum = 0.0;
            for (int j = 0; j < K; ++j) {
                if (j == i) continue;
                double m = rateKK(migrationRate, model, i, j);
                migOutSum += m;
                migInSum += m * E[j];
            }

            double totalOutRate = mu + psi + birthSum + migOutSum;
            out[i] = mu - totalOutRate * E[i] + migInSum + E[i] * quadSum;
        }
        return out;
    }

    /** Linear ODE (given E(t)) for a single lineage's per-type probability vector p(t). */
    private double[] derivativeP(double[] P, double[] E, int model) {
        int K = numTypes;
        double[] out = new double[K];
        for (int i = 0; i < K; ++i) {
            double mu = rateK(deathRate, model, i);
            double psi = rateK(serialSamplingRate, model, i);

            double birthSum = 0.0;
            double quadSum = 0.0;   // sum_j birthRate_ij * E_j
            double birthPSum = 0.0; // sum_j birthRate_ij * P_j
            for (int j = 0; j < K; ++j) {
                double b = rateKK(birthRate, model, i, j);
                birthSum += b;
                quadSum += b * E[j];
                birthPSum += b * P[j];
            }

            double migOutSum = 0.0;
            double migInSum = 0.0;
            for (int j = 0; j < K; ++j) {
                if (j == i) continue;
                double m = rateKK(migrationRate, model, i, j);
                migOutSum += m;
                migInSum += m * P[j];
            }

            double totalOutRate = mu + psi + birthSum + migOutSum;
            out[i] = -totalOutRate * P[i] + migInSum + P[i] * quadSum + E[i] * birthPSum;
        }
        return out;
    }

    private static double[] addScaled(double[] a, double[] b, double s) {
        double[] out = new double[a.length];
        for (int i = 0; i < a.length; ++i) {
            out[i] = a[i] + s * b[i];
        }
        return out;
    }

    private double[] rk4StepE(double[] y, double t, double dt, int model) {
        double[] k1 = derivativeE(y, model);
        double[] k2 = derivativeE(addScaled(y, k1, dt / 2), model);
        double[] k3 = derivativeE(addScaled(y, k2, dt / 2), model);
        double[] k4 = derivativeE(addScaled(y, k3, dt), model);
        double[] out = new double[y.length];
        for (int i = 0; i < y.length; ++i) {
            out[i] = y[i] + (dt / 6.0) * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]);
        }
        return out;
    }

    private double[] rk4StepP(double[] p, double t, double dt, int model) {
        double[] E1 = getE(t);
        double[] k1 = derivativeP(p, E1, model);
        double[] Emid = getE(t + dt / 2);
        double[] k2 = derivativeP(addScaled(p, k1, dt / 2), Emid, model);
        double[] k3 = derivativeP(addScaled(p, k2, dt / 2), Emid, model);
        double[] E4 = getE(t + dt);
        double[] k4 = derivativeP(addScaled(p, k3, dt), E4, model);
        double[] out = new double[p.length];
        for (int i = 0; i < p.length; ++i) {
            out[i] = p[i] + (dt / 6.0) * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]);
        }
        return out;
    }

    // ----------------------------------------------------------------------------------
    // E(t): computed once on a dense grid from t=0 (present) to the origin, then
    // interpolated. Handles "rho pulse" (intensive/present-day-style sampling) jump
    // conditions at interior epoch boundaries by storing duplicate (pre-/post-jump)
    // entries at the same time stamp.
    // ----------------------------------------------------------------------------------

    private void computeEGrid() {
        int K = numTypes;
        double origin = originTime.getParameterValue(0);

        List<Double> times = new ArrayList<Double>();
        List<double[]> values = new ArrayList<double[]>();

        double[] y = new double[K];
        for (int k = 0; k < K; ++k) {
            y[k] = 1.0 - rateK(samplingProbability, 0, k);
        }
        times.add(0.0);
        values.add(y.clone());

        double t = 0.0;
        int model = 0;
        while (t < origin - TIME_TOL) {
            double segEnd = (model < numIntervals - 1) ? Math.min(modelStartTimes[model + 1], origin) : origin;
            int steps = Math.max(1, (int) Math.ceil((segEnd - t) * stepsPerUnitTime));
            double dt = (segEnd - t) / steps;
            for (int s = 0; s < steps; ++s) {
                y = rk4StepE(y, t, dt, model);
                t += dt;
                times.add(t);
                values.add(y.clone());
            }
            if (model < numIntervals - 1 && Math.abs(t - modelStartTimes[model + 1]) < TIME_TOL) {
                int nextModel = model + 1;
                double[] jumped = new double[K];
                for (int k = 0; k < K; ++k) {
                    double rho = rateK(samplingProbability, nextModel, k);
                    jumped[k] = (1.0 - rho) * y[k];
                }
                y = jumped;
                times.add(t);
                values.add(y.clone());
                model = nextModel;
            }
        }

        eGridTimes = new double[times.size()];
        eGridValues = new double[values.size()][];
        for (int i = 0; i < times.size(); ++i) {
            eGridTimes[i] = times.get(i);
            eGridValues[i] = values.get(i);
        }
        eGridKnown = true;
    }

    private double[] getE(double t) {
        if (!eGridKnown) {
            computeEGrid();
        }
        int hi = eGridTimes.length - 1;
        if (t <= eGridTimes[0]) {
            return eGridValues[0].clone();
        }
        if (t >= eGridTimes[hi]) {
            return eGridValues[hi].clone();
        }
        int lo = 0;
        while (lo < hi) {
            int mid = (lo + hi + 1) / 2;
            if (eGridTimes[mid] <= t) {
                lo = mid;
            } else {
                hi = mid - 1;
            }
        }
        int idx = lo;
        double t0 = eGridTimes[idx];
        double t1 = eGridTimes[idx + 1];
        if (t1 - t0 < 1e-12) {
            return eGridValues[idx + 1].clone();
        }
        double frac = (t - t0) / (t1 - t0);
        double[] out = new double[numTypes];
        for (int k = 0; k < numTypes; ++k) {
            out[k] = eGridValues[idx][k] + frac * (eGridValues[idx + 1][k] - eGridValues[idx][k]);
        }
        return out;
    }

    // ----------------------------------------------------------------------------------
    // Advance a single lineage's probability vector along the linear ODE from t0 to t1
    // (t0 <= t1), stepping across epoch boundaries and applying rho-pulse jumps.
    // ----------------------------------------------------------------------------------

    private double[] advanceP(double[] y, double t0, double t1) {
        if (t1 <= t0 + TIME_TOL) {
            return y;
        }
        double t = t0;
        int model = epochIndexForTime(t0);
        double[] cur = y;
        while (t < t1 - TIME_TOL) {
            double segEnd = (model < numIntervals - 1) ? Math.min(modelStartTimes[model + 1], t1) : t1;
            int steps = Math.max(1, (int) Math.ceil((segEnd - t) * stepsPerUnitTime));
            double dt = (segEnd - t) / steps;
            for (int s = 0; s < steps; ++s) {
                cur = rk4StepP(cur, t, dt, model);
                t += dt;
            }
            if (model < numIntervals - 1
                    && t1 > modelStartTimes[model + 1] + TIME_TOL
                    && Math.abs(t - modelStartTimes[model + 1]) < TIME_TOL) {
                int nextModel = model + 1;
                double[] jumped = new double[numTypes];
                for (int k = 0; k < numTypes; ++k) {
                    double rho = rateK(samplingProbability, nextModel, k);
                    jumped[k] = (1.0 - rho) * cur[k];
                }
                cur = jumped;
                model = nextModel;
            }
        }
        return cur;
    }

    // ----------------------------------------------------------------------------------
    // Tree traversal
    // ----------------------------------------------------------------------------------

    private int resolveType(Taxon taxon) {
        Object attr = taxon.getAttribute(typeAttributeName);
        if (attr == null) {
            throw new IllegalArgumentException("Taxon '" + taxon.getId() + "' has no '" +
                    typeAttributeName + "' attribute giving its type");
        }
        if (attr instanceof Number) {
            return ((Number) attr).intValue();
        }
        String s = attr.toString().trim();
        for (int k = 0; k < numTypes; ++k) {
            if (typeLabels[k].equals(s)) {
                return k;
            }
        }
        throw new IllegalArgumentException("Unknown type label '" + s + "' for taxon '" + taxon.getId() +
                "'; expected one of " + Arrays.toString(typeLabels));
    }

    private double[] initializeTipVector(Tree tree, NodeRef node, double height) {
        Taxon taxon = tree.getNodeTaxon(node);
        int type = resolveType(taxon);
        int model = epochIndexForTime(height);
        double[] vec = new double[numTypes];
        if (height < TIME_TOL) {
            double rho = rateK(samplingProbability, model, type);
            if (!(rho > 0.0)) {
                throw new IllegalArgumentException("Tip '" + taxon.getId() +
                        "' is sampled at time 0 but samplingProbability for its type is 0");
            }
            vec[type] = rho;
        } else {
            double psi = rateK(serialSamplingRate, model, type);
            double r = rateK(treatmentProbability, model, type);
            double[] Et = getE(height);
            vec[type] = psi * (r + (1.0 - r) * Et[type]);
        }
        return vec;
    }

    private double[] combineBirth(double[] v1, double[] v2, double height) {
        int model = epochIndexForTime(height);
        double[] out = new double[numTypes];
        for (int i = 0; i < numTypes; ++i) {
            double sum = rateKK(birthRate, model, i, i) * v1[i] * v2[i];
            for (int j = 0; j < numTypes; ++j) {
                if (j == i) continue;
                double b = rateKK(birthRate, model, i, j);
                sum += b * (v1[i] * v2[j] + v1[j] * v2[i]);
            }
            out[i] = sum;
        }
        return out;
    }

    private double[] computeVectorAt(Tree tree, NodeRef node) {
        double height = tree.getNodeHeight(node);
        if (tree.isExternal(node)) {
            return initializeTipVector(tree, node, height);
        }
        int childCount = tree.getChildCount(node);
        if (childCount != 2) {
            throw new UnsupportedOperationException(
                    "MultiTypeBirthDeathSerialSamplingModel currently requires strictly bifurcating trees " +
                            "(encountered a node with " + childCount + " children; sampled-ancestor / direct-ancestor " +
                            "nodes are not yet supported in this prototype)");
        }
        NodeRef c1 = tree.getChild(node, 0);
        NodeRef c2 = tree.getChild(node, 1);
        double[] v1 = advanceP(computeVectorAt(tree, c1), tree.getNodeHeight(c1), height);
        double[] v2 = advanceP(computeVectorAt(tree, c2), tree.getNodeHeight(c2), height);
        return combineBirth(v1, v2, height);
    }

    @Override
    public double calculateTreeLogLikelihood(Tree tree) {
        if (!eGridKnown) {
            computeEGrid();
        }
        NodeRef root = tree.getRoot();
        double rootHeight = tree.getNodeHeight(root);
        double origin = originTime.getParameterValue(0);
        if (origin < rootHeight) {
            return Double.NEGATIVE_INFINITY;
        }

        double[] rootVec = computeVectorAt(tree, root);
        double[] atOrigin = advanceP(rootVec, rootHeight, origin);

        double total = 0.0;
        for (int k = 0; k < numTypes; ++k) {
            total += atOrigin[k];
        }
        if (!(total > 0.0) || Double.isNaN(total)) {
            return Double.NEGATIVE_INFINITY;
        }
        double logL = Math.log(total);

        if (conditionOnSurvival) {
            double[] eOrigin = getE(origin);
            double norm = 0.0;
            for (int k = 0; k < numTypes; ++k) {
                norm += (1.0 - eOrigin[k]);
            }
            if (!(norm > 0.0)) {
                return Double.NEGATIVE_INFINITY;
            }
            logL -= Math.log(norm);
        }

        return logL;
    }

    @Override
    public double calculateTreeLogLikelihood(Tree tree, Set<Taxon> exclude) {
        if (exclude != null && exclude.size() > 0) {
            throw new UnsupportedOperationException(
                    "Excluding taxa is not yet supported by MultiTypeBirthDeathSerialSamplingModel");
        }
        return calculateTreeLogLikelihood(tree);
    }

    public int getNumTypes() {
        return numTypes;
    }

    public String[] getTypeLabels() {
        return typeLabels;
    }

    public String getDescription() {
        return "A prototype episodic multi-type birth-death-sampling model (likelihood only, no gradient yet).";
    }
}
