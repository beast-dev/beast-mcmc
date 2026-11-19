/*
 * ListedPairSwapOperator.java
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

import java.util.List;

/**
 * Swaps the values at two indices drawn uniformly from a fixed list of index pairs, in every parameter given
 * (e.g. CTMC rates and their Bayesian-bridge local scales, so that a prior exchangeable over those tuples is
 * unchanged). The list is fixed, so the proposal is symmetric (Hastings ratio 1).
 *
 * @author Marc A. Suchard
 */
public class ListedPairSwapOperator extends SimpleMCMCOperator {

    public static final String OPERATOR_NAME = "listedPairSwapOperator";

    /**
     * @param pairs zero-based indices, pair p is (pairs[2p], pairs[2p + 1])
     */
    public ListedPairSwapOperator(List<Parameter> parameters, int[] pairs, double weight) {
        super(weight);
        this.parameters = parameters;
        this.pairs = pairs;
    }

    @Override
    public double doOperation() {
        int p = MathUtils.nextInt(pairs.length / 2);
        int a = pairs[2 * p];
        int b = pairs[2 * p + 1];
        for (Parameter parameter : parameters) {
            double va = parameter.getParameterValue(a);
            double vb = parameter.getParameterValue(b);
            parameter.setParameterValue(a, vb);
            parameter.setParameterValue(b, va);
        }
        return 0.0;
    }

    @Override
    public String getOperatorName() {
        return OPERATOR_NAME + "(" + parameters.get(0).getParameterName() + ")";
    }

    private final List<Parameter> parameters;
    private final int[] pairs;
}
