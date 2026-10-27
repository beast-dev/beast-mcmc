/*
 * OperatorPerturbation.java
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

package dr.evomodel.tree.search;

import dr.evomodel.operators.ExchangeOperator;
import dr.evomodel.operators.NNI;
import dr.evomodel.operators.SubtreeLeapOperator;
import dr.evomodel.operators.WilsonBalding;
import dr.evomodel.tree.TreeModel;
import dr.inference.operators.AdaptationMode;
import dr.inference.operators.SimpleMCMCOperator;
import dr.math.MathUtils;

import java.util.ArrayList;
import java.util.List;

/**
 * Perturbation of the search (plan 7 §5.9, §7.2) by BEAST's own topology operators, called through their doOperation as
 * MCMC would and with the Hastings ratio ignored: rooted NNI with a uniform height (NNI), wide exchange
 * (ExchangeOperator), subtree leap (SubtreeLeapOperator) and Wilson-Balding (WilsonBalding). Each draw is kept only if
 * the operator does not reject it, the constraints hold and the objective stays finite; otherwise it is undone and
 * another drawn. A draw may give the root another node (subtree leap above the root): snapshots taken before it
 * restore by clade. A perturbed tree is
 * always hill-climbed and compared with the best, never accepted on its own.
 *
 * @author Xiang Ji
 * @author Marc A. Suchard
 */
public final class OperatorPerturbation {

    private static final int MAX_DRAWS_PER_MOVE = 20;

    private final TreeModel tree;
    private final TrialEvaluator evaluator;
    private final TopologyConstraints constraints;
    private final List<SimpleMCMCOperator> operators = new ArrayList<SimpleMCMCOperator>();
    private final List<String> names = new ArrayList<String>();
    private final SubtreeLeapOperator leap;
    private final double leapSize;

    private long draws = 0;
    private long kept = 0;

    /**
     * @param operatorNames comma-separated: nni, narrow, wide, leap, wb
     * @param leapSize      SubtreeLeap's distance scale relative to the root height
     */
    public OperatorPerturbation(TreeModel tree, TrialEvaluator evaluator, TopologyConstraints constraints,
                                String operatorNames, double leapSize) {
        this.tree = tree;
        this.evaluator = evaluator;
        this.constraints = constraints;
        this.leapSize = leapSize;
        SubtreeLeapOperator leap = null;
        for (String name : operatorNames.replace(" ", "").split(",")) {
            if (name.isEmpty()) {
                continue;
            }
            final SimpleMCMCOperator operator;
            if (name.equals("nni")) {
                operator = new NNI(tree, 1.0);
            } else if (name.equals("narrow")) {
                operator = new ExchangeOperator(ExchangeOperator.NARROW, tree, 1.0);
            } else if (name.equals("wide")) {
                operator = new ExchangeOperator(ExchangeOperator.WIDE, tree, 1.0);
            } else if (name.equals("leap")) {
                leap = new SubtreeLeapOperator(tree, 1.0, 1.0, SubtreeLeapOperator.DistanceKernelType.NORMAL,
                        AdaptationMode.ADAPTATION_OFF, 0.234);
                operator = leap;
            } else if (name.equals("wb")) {
                operator = new WilsonBalding(tree, 1.0);
            } else {
                throw new IllegalArgumentException("Unknown perturbation operator " + name +
                        " (nni, narrow, wide, leap or wb)");
            }
            operators.add(operator);
            names.add(name);
        }
        if (operators.isEmpty()) {
            throw new IllegalArgumentException("No perturbation operator");
        }
        this.leap = leap;
    }

    /**
     * Applies up to moves operator draws.
     *
     * @return the number of draws kept
     */
    public int perturb(int moves) {
        if (leap != null) {
            leap.setSize(leapSize * tree.getNodeHeight(tree.getRoot()));
        }
        int done = 0;
        for (int m = 0; m < moves; ++m) {
            for (int attempt = 0; attempt < MAX_DRAWS_PER_MOVE; ++attempt) {
                final SimpleMCMCOperator operator = operators.get(MathUtils.nextInt(operators.size()));
                draws++;
                final boolean ok = evaluator.tryEdit(new TrialEvaluator.Edit() {
                    public boolean apply() {
                        final double logHastings = operator.doOperation();
                        return logHastings > Double.NEGATIVE_INFINITY && constraints.areSatisfied(tree);
                    }
                });
                if (ok) {
                    kept++;
                    done++;
                    break;
                }
            }
        }
        return done;
    }

    public List<String> getOperatorNames() {
        return names;
    }

    public long getDraws() {
        return draws;
    }

    public long getKept() {
        return kept;
    }
}
