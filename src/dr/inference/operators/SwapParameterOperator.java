/*
 * SwapParameterOperator.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
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
 */

package dr.inference.operators;

import dr.inference.model.Bounds;
import dr.inference.model.Parameter;
import dr.inferencexml.operators.SwapOperatorParser;
import dr.math.MathUtils;

import java.util.List;

/**
 * A generic operator swapping all values (dimensions) between two or more parameters
 *
 * @author Marc A. Suchard
 */
public class SwapParameterOperator extends SimpleMCMCOperator {

    public SwapParameterOperator(List<Parameter> parameterList, double weight) {
        this.parameterList = parameterList;

        if (parameterList.size() < 2) {
            throw new IllegalArgumentException("More than 1 parameter is needed");
        }

        int size = parameterList.get(0).getDimension();
        for (int i = 1; i < parameterList.size(); ++i) {
            if (size != parameterList.get(i).getDimension()) {
                throw new IllegalArgumentException("All parameters most be the same size");
            }
        }

        setWeight(weight);
    }


    /**
     * swap all values in two random parameters.
     */
    public final double doOperation() {

        int i = MathUtils.nextInt(parameterList.size());
        int j = i;

        while (j == i) {
            j = MathUtils.nextInt(parameterList.size());
        }

        final Parameter a = parameterList.get(i);
        final Parameter b = parameterList.get(j);

        Bounds<Double> aBounds = a.getBounds();
        Bounds<Double> bBounds = b.getBounds();


        for (int k = 0; k < a.getDimension(); ++k) {
            final double ak = a.getParameterValue(k);
            final double bk = b.getParameterValue(k);

            // Check bk outside of aBounds or ak outside of bBounds
            if (isOutside(aBounds, bk, k) || isOutside(bBounds, ak, k)) {
//                throw new OperatorFailedException("proposed value outside boundaries");
                return Double.NEGATIVE_INFINITY;
            }

            // Swap
            a.setParameterValueQuietly(k, bk);
            b.setParameterValueQuietly(k, ak);
        }

        a.fireParameterChangedEvent();
        b.fireParameterChangedEvent();

        return 0.0;
    }

    private boolean isOutside(final Bounds<Double> bounds, final double x, final int index) {
        return (x < bounds.getLowerLimit(index) || x > bounds.getUpperLimit(index));
    }

    private String getParameterNames() {
        if (parameterNames == null) {
            StringBuilder sb = new StringBuilder();
            for (Parameter p : parameterList) {
                sb.append(p.getParameterName()).append(".");
            }
            parameterNames = sb.toString();
        }
        return parameterNames;
    }

    public String getOperatorName() {
        return SwapOperatorParser.SWAP_OPERATOR + "(" + getParameterNames() + "swap)";
    }

    public String getPerformanceSuggestion() {
        return "No suggestions";
    }

    private final List<Parameter> parameterList;
    private String parameterNames = null;
}
