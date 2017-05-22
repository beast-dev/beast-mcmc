/*
 * SelectorOperator.java
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

import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.inferencexml.operators.SelectorOperatorParser;
import dr.math.MathUtils;

import java.util.ArrayList;
import java.util.List;

/**
 * The code is much more elegant in Python.
 * I definitly don't have to write my own max for a list of integers (yikes).
 *
 * @author Joseph Heled
 *         Date: 4/09/2009
 */
public class SelectorOperator extends SimpleMCMCOperator {
    private final Parameter selector;

    private final int[] np_m1;

    private final int[] np_m2;

    public SelectorOperator(Parameter selector) {
        this.selector = selector;
        final int len = selector.getSize();
        np_m1 = new int[len +1];

        for(int l = 0; l < np_m1.length; ++l) {
            np_m1[l] = npos(len, l);
        }

        np_m2 = new int[len + 1];
        np_m2[0] = 1;
        for(int mx = 0; mx < len; ++mx) {
           np_m2[mx+1] = 0;
           for(int l = 1; l < len+1; ++l) {
               np_m2[mx+1] += npos(l, mx);
           }
        }
    }
    
    public String getOperatorName() {
        return SelectorOperatorParser.SELECTOR_OPERATOR + "(" + selector.getParameterName() + ")";
    }

    public double doOperation() {

        final int[] s = vals();
        final List<Integer> poss = movesFrom_m2(s);
        final int i = MathUtils.nextInt(poss.size()/2);

        final int[] y = new int[s.length];
        System.arraycopy(s, 0, y, 0, s.length);
        final Integer p = poss.get(2 * i);
        y[p] = poss.get(2*i+1);

        double hr = count_sr_m2(s, y);
        hr *= (double)(poss.size()* np_m2[max(s)+1])/(movesFrom_m2(y).size() * np_m2[max(y)+1]);

        selector.setParameterValue(p, y[p]);
        
        return Math.log(hr);
    }

//    public double doOperation_m1() throws OperatorFailedException {
//
//        final int[] s = vals();
//        final List<Integer> poss = movesFrom_m1(s);
//        final int i = MathUtils.nextInt(poss.size()/2);
//
//        final int[] y = new int[s.length];
//        System.arraycopy(s, 0, y, 0, s.length);
//        final Integer p = poss.get(2 * i);
//        y[p] = poss.get(2*i+1);
//
//        double hr = count_sr_m1(s, y);
//        hr *= (double)(poss.size()* np_m1[max(s)])/(movesFrom_m1(y).size() * np_m1[max(y)]);
//
//        selector.setParameterValue(p, y[p]);
//
//        return Math.log(hr);
//    }

    public String getPerformanceSuggestion() {
        return null;
    }

    private int[] vals() {
        return intVals(selector);
    }

    static int[] intVals(Variable<Double> var) {
        int[] v = new int[var.getSize()];
        for(int k = 0; k < v.length; ++k) {
            final double vk = var.getValue(k);
            v[k] = (int)(vk + ((vk>= 0) ? 0.5 : -0.5));
        }
        return v;
    }

    private List<Integer> movesFrom_m2(int[] s) {
        final int mx = max(s);

        final int[] counts = counts_used_m2(s, mx);
        final List<Integer> opt = new ArrayList<Integer>(5);

        for(int k = 0; k < s.length; ++k) {
            final int si = s[k];
            if( si < 0 ) {
                opt.add(k); opt.add(0);
                if( mx >= 0 ) {
                    opt.add(k); opt.add(mx+1);
                }
                for(int x = 1; x < mx+1; ++x) {
                    if( counts[x]+1 <= counts[x-1]) {
                        opt.add(k); opt.add(x);
                    }
                }
            } else {
                if( si < mx && ((counts[si] == 1) || counts[si] == counts[si+1]) ) {
                    // only or breaks order -> no moves
                } else {
                    for(int x = 0; x < mx+1; ++x) {
                        if(x == si) {
                            continue;
                        }
                        if( (x > si && counts[si] - 1 >= counts[x] + 1 && counts[x-1] >= counts[x]+1)
                                ||
                                (x < si && (x > 0 && counts[x]+1 <= counts[x-1] || x == 0)) ) {
                            opt.add(k); opt.add(x);
                        }
                    }
                    if( counts[si] > 1) {
                        opt.add(k); opt.add(mx+1);
                    }
                    opt.add(k); opt.add(-1);
                }
            }
        }
        return opt;
    }

    private List<Integer> movesFrom_m1(int[] s) {
        final int mx = max(s);

        final int[] counts = counts_m1(s, mx);
        final List<Integer> opt = new ArrayList<Integer>(5);

        for(int k = 0; k < s.length; ++k) {
            final int si = s[k];
            if( si<mx && ((counts[si] == 1) || counts[si] == counts[si+1]) ) {
                // only or breaks order -> no moves
            } else {
                for(int x = 0; x < mx+1; ++x) {
                    if(x == si) {
                        continue;
                    }
                    if( (x > si && counts[si] - 1 >= counts[x] + 1 && counts[x-1] >= counts[x]+1)
                            ||
                            (x < si && (x > 0 && counts[x]+1 <= counts[x-1] || x == 0)) ) {
                        opt.add(k);
                        opt.add(x);
                    }
                }
                if( counts[si] > 1) {
                    opt.add(k);
                    opt.add(mx+1);
                }
            }
        }
        return opt;
    }

    private static int npos(int s, int m) {
        return npos(s, m, 1);
    }

    private static int npos(int s, int m, int mn) {
        if( m == 0 || s == 0 ) {
            return 1;
        }

        int tot = 0;
        for(int k = mn; k < 1+s/m; ++k) {
            final int r = s - k*(m+1);
            if(r < 0 ) {
                break;
            }
            tot += npos(r, m-1, 0);
        }
        return tot;
    }

    private static int sum(int[] s) {
        int sum = 0;
        for(int si : s) {
            sum += si;
        }
        return sum;
    }

    private static int max(int[] s) {
        int mx = s[0];

        for(int k = 1; k < s.length; ++k) {
            if( mx < s[k] ) {
                mx = s[k];
            }
        }
        return mx;
    }

    private static int[] counts_m1(int[] s, int mx) {
        int[] c = new int[mx+1];
        for(int si : s) {
            c[si]++;
        }
        return c;
    }

    // Counts in s including unused (-1). indices are shifted by 1
    private static int[] counts_m2(int[] s, int mx) {
        int[] c = new int[mx+2];
        for(int si : s) {
            c[si+1]++;
        }
        return c;
    }

    static int[] counts_m2(int[] s) {
        return counts_m2(s, max(s));
    }

    private static int[] counts_used_m2(int[] s, int mx) {
        int[] c = new int[mx+1];
        for(int si : s) {
            if( si >= 0 ) {
                c[si]++;
            }
        }
        return c;
    }

    static int[] counts_used_m2(int[] s) {
        return  counts_used_m2(s, max(s));
    }

    private static long choose(int n, int k) {
        double r = 1;
        while( n > k ) {
            r *= n;
            r /= (n-k);
            --n;
        }
        return (long)(r+0.5);
    }

   private static long[] countl_m1(int[] ls) {
       int l = sum(ls);
       int i = 0;
       long[] r = new long[ls.length];

       while(l > 0) {
           r[i] = choose(l, ls[i]);

           l -= ls[i];
           i += 1;
       }
       return r;
   }

    private static double count_sr_m1(int[] x, int[] y) {
        long[] r1 = countl_m1(counts_m1(x, max(x)));
        long[] r2 = countl_m1(counts_m1(y, max(y)));

        int k = Math.min(r1.length, r2.length);
        double r = 1;
        for(int i=0; i < k; ++i) {
            r *= r1[i];
            r /= r2[i];
        }
        for(int i=k; i < r1.length; ++i) {
            r *= r1[i];
        }
        for(int i=k; i < r2.length; ++i) {
            r /= r2[i];
        }
        return r;
    }

    private static long[] countl_m2(int[] ls) {
        if( ls.length == 1 ) {
            return new long[]{1};
        }

        int l = sum(ls);

        long[] r = new long[ls.length];
        r[0] = choose(l, ls[0]);
        l -= ls[0];
        int i = 1;

        while( l > 0 ) {
            r[i] = choose(l, ls[i]);

            l -= ls[i];
            i += 1;
        }
        return r;
    }

    private static double count_sr_m2(int[] x, int[] y) {
        long[] r1 = countl_m2(counts_m2(x));
        long[] r2 = countl_m2(counts_m2(y));

        int k = Math.min(r1.length, r2.length);
        double r = 1;
        for(int i=0; i < k; ++i) {
            r *= r1[i];
            r /= r2[i];
        }
        for(int i=k; i < r1.length; ++i) {
            r *= r1[i];
        }
        for(int i=k; i < r2.length; ++i) {
            r /= r2[i];
        }
        return r;
    }

}
