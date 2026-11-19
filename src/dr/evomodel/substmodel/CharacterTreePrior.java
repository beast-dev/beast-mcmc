/*
 * CharacterTreePrior.java
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

package dr.evomodel.substmodel;

import dr.inference.model.Likelihood;

/**
 * The prior on the topology of a sampled character tree: log pi_T(q(M)) - log R(M), with q the number of non-root
 * internal nodes of the contracted tree M and R(M) the number of binary trees that resolve it. Dividing by R(M)
 * corrects for the skeleton: exactly R(M) pairs (binary skeleton, collapse state) give M, so the marginal prior on M is
 * pi_T. Types, each a function of q:
 * uniform: every multifurcating tree has probability 1 / N(S) (Schroeder's fourth problem);
 * polytomy: weight c^(S - 2 - q), the polytomy prior of Lewis, Holder and Holsinger (2005), c > 1 favours polytomies;
 * resolutionClass: probability 1 / (S - 1) for each q, uniform within a class.
 * Tree counts by q come from T(n + 1, k) = k T(n, k) + (n + k - 1) T(n, k - 1), T(2, 1) = 1, in log space.
 *
 * @author Xiang Ji
 * @author Marc A. Suchard
 */
public class CharacterTreePrior extends Likelihood.Abstract {

    public enum Type {
        UNIFORM, POLYTOMY, RESOLUTION_CLASS;

        public static Type parse(String value) {
            if ("uniform".equalsIgnoreCase(value)) return UNIFORM;
            if ("polytomy".equalsIgnoreCase(value)) return POLYTOMY;
            if ("resolutionClass".equalsIgnoreCase(value)) return RESOLUTION_CLASS;
            throw new IllegalArgumentException("type must be uniform, polytomy or resolutionClass; found '" + value + "'");
        }
    }

    private final CharacterTreeModel characterTree;
    private final Type type;
    private final int stateCount;
    private final double logC;
    private final double[] logTreesByQ;
    private final double logTreeCount;
    private final double logPolytomyNormaliser;

    public CharacterTreePrior(CharacterTreeModel characterTree, Type type, double c) {
        super(characterTree);
        if (type == Type.POLYTOMY && !(c > 0.0)) {
            throw new IllegalArgumentException("the polytomy prior needs c > 0");
        }
        this.characterTree = characterTree;
        this.type = type;
        this.stateCount = characterTree.getExternalNodeCount();
        this.logC = Math.log(c);
        this.logTreesByQ = logTreesByQ(stateCount);
        this.logTreeCount = logSumExp(logTreesByQ, 0.0);
        final double[] weighted = new double[logTreesByQ.length];
        for (int q = 0; q < weighted.length; q++) {
            weighted[q] = logTreesByQ[q] + (stateCount - 2 - q) * logC;
        }
        this.logPolytomyNormaliser = logSumExp(weighted, 0.0);
    }

    /**
     * log T(S, q + 1) for q = 0 .. S - 2: the log number of rooted trees on S labelled tips with q + 1 internal
     * nodes, each with at least two children.
     */
    public static double[] logTreesByQ(int stateCount) {
        double[] row = new double[]{Double.NEGATIVE_INFINITY, 0.0};
        for (int n = 2; n < stateCount; n++) {
            final double[] next = new double[n + 1];
            next[0] = Double.NEGATIVE_INFINITY;
            for (int k = 1; k <= n; k++) {
                final double a = k < row.length ? Math.log(k) + row[k] : Double.NEGATIVE_INFINITY;
                final double b = k - 1 < row.length ? Math.log(n + k - 1) + row[k - 1] : Double.NEGATIVE_INFINITY;
                next[k] = logAdd(a, b);
            }
            row = next;
        }
        final double[] out = new double[Math.max(1, stateCount - 1)];
        for (int q = 0; q < out.length; q++) {
            out[q] = q + 1 < row.length ? row[q + 1] : Double.NEGATIVE_INFINITY;
        }
        return out;
    }

    private static double logAdd(double a, double b) {
        if (a == Double.NEGATIVE_INFINITY) return b;
        if (b == Double.NEGATIVE_INFINITY) return a;
        final double max = Math.max(a, b);
        return max + Math.log(Math.exp(a - max) + Math.exp(b - max));
    }

    private static double logSumExp(double[] values, double start) {
        double sum = Double.NEGATIVE_INFINITY;
        for (double v : values) {
            sum = logAdd(sum, v);
        }
        return sum;
    }

    public double getLogTreeCount() {
        return logTreeCount;
    }

    /**
     * The log prior probability of one multifurcating tree with q non-root internal nodes.
     */
    public double logTreeProbability(int q) {
        switch (type) {
            case UNIFORM:
                return -logTreeCount;
            case POLYTOMY:
                return (stateCount - 2 - q) * logC - logPolytomyNormaliser;
            default:
                return -Math.log(stateCount - 1) - logTreesByQ[q];
        }
    }

    @Override
    protected double calculateLogLikelihood() {
        final int q = characterTree.getContractedInternalNodeCount() - 1;
        return logTreeProbability(q) - characterTree.getLogResolutionCount();
    }

    @Override
    protected boolean getLikelihoodKnown() {
        return false;
    }

    @Override
    public String prettyName() {
        return "characterTreePrior(" + characterTree.getId() + ", " + type.name().toLowerCase() + ")";
    }
}
