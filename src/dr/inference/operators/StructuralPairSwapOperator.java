/*
 * StructuralPairSwapOperator.java
 *
 * Copyright © 2002-2026 the BEAST Development Team
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
 * BEAST is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package dr.inference.operators;

import dr.inference.model.Parameter;
import dr.math.MathUtils;

import java.io.PrintWriter;

/**
 * List-free swap of an asymmetric CTMC rate theta_a (q = theta^2, rate i->j) and its Bayesian-bridge local scale
 * with those of a structurally related rate b: its reverse j->i (class R), a rate k->j into the same target (T)
 * or a rate i->k out of the same source (S).
 * <p>
 * The anchor a is drawn with probability proportional to min(q_a, anchorCap)^2 and the class independently of the
 * state. For T and S, the third location k is drawn, with probability linkProbability, uniformly from C_a, the
 * locations linked to the changing endpoint e (i for T, j for S) by an on rate (q > onThreshold) in either
 * direction, and otherwise uniformly from all K - 2 locations other than i and j.
 * <p>
 * The swap keeps the multiset of values, so the anchor terms cancel, and the link sets C_a and C_b never involve
 * rates a or b, so they are the same before and after the swap. The log Hastings ratio is therefore
 * log rho(|C_b|, L) - log rho(|C_a|, L), with L = [k in C_a] and
 * rho(n, L) = linkProbability * L / n + (1 - linkProbability) / (K - 2) for n > 0, 1 / (K - 2) otherwise.
 *
 * @author Marc A. Suchard
 */
public class StructuralPairSwapOperator extends SimpleMCMCOperator {

    public static final String OPERATOR_NAME = "structuralPairSwapOperator";

    static final int REVERSE = 0;
    static final int SAME_TARGET = 1;
    static final int SAME_SOURCE = 2;
    private static final String[] CLASS_NAMES = {"R", "T", "S"};

    /**
     * @param stateCount K, with rates and localScale of dimension K(K - 1) in the order of
     *                   AbstractLogAdditiveSubstitutionModelGradient.makeAsymmetricMap
     * @param trace      one line per proposal, or null
     * The same-source class has probability 1 - reverseProbability - sameTargetProbability.
     */
    public StructuralPairSwapOperator(Parameter rates, Parameter localScale, int stateCount,
                                      double onThreshold, double anchorCap,
                                      double reverseProbability, double sameTargetProbability,
                                      double linkProbability, PrintWriter trace, double weight) {
        super(weight);
        this.rates = rates;
        this.localScale = localScale;
        this.stateCount = stateCount;
        this.halfDim = stateCount * (stateCount - 1) / 2;
        this.onThreshold = onThreshold;
        this.anchorCap = anchorCap;
        this.reverseProbability = reverseProbability;
        this.sameTargetProbability = sameTargetProbability;
        this.linkProbability = linkProbability;
        this.trace = trace;

        int dim = 2 * halfDim;
        src = new int[dim];
        dst = new int[dim];
        int k = 0;
        for (int i = 0; i < stateCount; ++i) {
            for (int j = i + 1; j < stateCount; ++j) {
                src[k] = i;
                dst[k++] = j;
            }
        }
        for (int j = 0; j < stateCount; ++j) {
            for (int i = j + 1; i < stateCount; ++i) {
                src[k] = i;
                dst[k++] = j;
            }
        }
        anchorWeights = new double[dim];
        linked = new int[stateCount];

        if (trace != null) {
            trace.println("a\tb\tclass\tq_a\tq_b\th\taccepted");
        }
    }

    @Override
    public double doOperation() {
        for (int k = 0; k < anchorWeights.length; ++k) {
            double v = Math.min(q(k), anchorCap);
            anchorWeights[k] = v * v;
        }
        final int a = MathUtils.randomChoicePDF(anchorWeights);
        final int i = src[a];
        final int j = dst[a];

        final double u = MathUtils.nextDouble();
        final int cls = u < reverseProbability ? REVERSE :
                (u < reverseProbability + sameTargetProbability ? SAME_TARGET : SAME_SOURCE);

        int k = -1;
        if (cls != REVERSE) {
            int nA = fillLinked(cls == SAME_TARGET ? i : j, cls == SAME_TARGET ? j : i);
            if (nA > 0 && MathUtils.nextDouble() < linkProbability) {
                k = linked[MathUtils.nextInt(nA)];
            } else {
                k = MathUtils.nextInt(stateCount - 2);
                if (k >= Math.min(i, j)) ++k;
                if (k >= Math.max(i, j)) ++k;
            }
        }
        return swap(a, cls, k);
    }

    /**
     * Swaps anchor a with its class-cls partner, whose third location is k (unused for REVERSE),
     * and returns the log Hastings ratio.
     */
    double swap(int a, int cls, int k) {
        final int i = src[a];
        final int j = dst[a];
        final int b;
        double logRatio = 0.0;
        if (cls == REVERSE) {
            b = index(j, i);
        } else {
            final int e = (cls == SAME_TARGET) ? i : j; // endpoint of a that b does not share
            final int f = (cls == SAME_TARGET) ? j : i; // endpoint shared by a and b
            b = (cls == SAME_TARGET) ? index(k, j) : index(i, k);
            final boolean link = isLinked(e, k);
            logRatio = Math.log(rho(fillLinked(k, f), link)) - Math.log(rho(fillLinked(e, f), link));
        }

        if (trace != null) {
            lastMove = (a + 1) + "\t" + (b + 1) + "\t" + CLASS_NAMES[cls] + "\t" + q(a) + "\t" + q(b) + "\t" + logRatio;
        }

        swapValues(rates, a, b);
        swapValues(localScale, a, b);
        return logRatio;
    }

    @Override
    public void accept(double deviation) {
        super.accept(deviation);
        if (trace != null) {
            trace.println(lastMove + "\t1");
        }
    }

    @Override
    public void reject() {
        super.reject();
        if (trace != null) {
            trace.println(lastMove + "\t0");
        }
    }

    @Override
    protected String getBaseOperatorName() {
        return OPERATOR_NAME + "(" + rates.getParameterName() + ")";
    }

    private static void swapValues(Parameter parameter, int a, int b) {
        final double va = parameter.getParameterValue(a);
        parameter.setParameterValueQuietly(a, parameter.getParameterValue(b));
        parameter.setParameterValueQuietly(b, va);
        parameter.fireParameterChangedEvent();
    }

    private double q(int k) {
        final double theta = rates.getParameterValue(k);
        return theta * theta;
    }

    private int index(int u, int v) {
        return u < v ? upper(u, v) : halfDim + upper(v, u);
    }

    private int upper(int u, int v) {
        return u * stateCount - u * (u + 1) / 2 + (v - u - 1);
    }

    private boolean isLinked(int e, int m) {
        return q(index(e, m)) > onThreshold || q(index(m, e)) > onThreshold;
    }

    /**
     * Fills linked[] with the locations m not in {e, f} linked to e and returns their number.
     */
    private int fillLinked(int e, int f) {
        int n = 0;
        for (int m = 0; m < stateCount; ++m) {
            if (m != e && m != f && isLinked(e, m)) {
                linked[n++] = m;
            }
        }
        return n;
    }

    private double rho(int n, boolean link) {
        final double uniform = 1.0 / (stateCount - 2);
        return n > 0 ? linkProbability * (link ? 1.0 : 0.0) / n + (1.0 - linkProbability) * uniform : uniform;
    }

    private final Parameter rates;
    private final Parameter localScale;
    private final int stateCount;
    private final int halfDim;
    private final double onThreshold;
    private final double anchorCap;
    private final double reverseProbability;
    private final double sameTargetProbability;
    private final double linkProbability;
    private final PrintWriter trace;

    private final int[] src;
    private final int[] dst;
    private final double[] anchorWeights;
    private final int[] linked;
    private String lastMove;
}
