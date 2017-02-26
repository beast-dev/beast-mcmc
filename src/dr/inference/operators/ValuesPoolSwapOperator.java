/*
 * ValuesPoolSwapOperator.java
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

import dr.inference.model.ValuesPool;
import dr.inference.model.Variable;
import dr.inferencexml.operators.ValuesPoolSwapOperatorParser;
import dr.math.MathUtils;

/**
 * @author Joseph Heled
 *         Date: 8/09/2009
 */
public class ValuesPoolSwapOperator extends SimpleMCMCOperator {

    private final ValuesPool pool;

    public ValuesPoolSwapOperator(ValuesPool pool) {
        this.pool = pool;
    }

    public String getOperatorName() {
        return ValuesPoolSwapOperatorParser.VALUESPOOL_OPERATOR + "(" + pool.getModelName() + ")";
    }

    public double doOperation() {
        final Variable<Double> selector = pool.getSelector();
        final int[] ints = SelectorOperator.intVals(selector);
        int[] c = SelectorOperator.counts_used_m2(ints);

        int n = 0;
        for(int k = 0; k < c.length-1; ++k) {
            if( c[k] == c[k+1] ) {
              ++n;
            }
        }

        if( n == 0 ) {
            throw new RuntimeException("No moves");
        }

        int j = MathUtils.nextInt(n);
        int p;
        for(p = 0; p < c.length-1; ++p) {
            if( c[p] == c[p+1] ) {
                if( j == 0 ) {
                    break;
                }
                j -= 1;
            }
        }

        int ip = -1, ip1 = -1;
        // exchange p's and p+1's
        int count = c[p];
        while( count > 0 ) {
            while( ints[++ip] != p );
            while( ints[++ip1] != p+1 );
            selector.setValue(ip, (double)(p+1));
            selector.setValue(ip1, (double)(p));
            --count;
        }
        final Variable<Double> vals = pool.getPool();
        final Double vp = vals.getValue(p);
        final Double vp1 = vals.getValue(p+1);
        vals.setValue(p, vp1);
        vals.setValue(p+1, vp);

        return 0;
    }

    public String getPerformanceSuggestion() {
        return null;
    }

}
