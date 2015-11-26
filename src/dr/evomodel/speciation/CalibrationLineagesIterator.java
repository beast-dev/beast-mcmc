/*
 * CalibrationLineagesIterator.java
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

package dr.evomodel.speciation;

/**
 * Created by IntelliJ IDEA.
 * User: joseph
 * Date: 22/06/11
 * Time: 12:39 PM
 * To change this template use File | Settings | File Templates.
 */
public class CalibrationLineagesIterator {
    final int[][] taxaPartialOrder;
    final int[] cladesFreeLins;
    private final LinsIterator[] iters;
    private int nCurIters;
    private int[][] vals;
    private int nFreeLineages;
    private final int[] maximalClades;

    CalibrationLineagesIterator(int[][] clades, int[][] taxaPartialOrder, boolean[] maximal, int externalNodeCount) {
        cladesFreeLins = new int[clades.length];
        for(int k = 0; k < cladesFreeLins.length; ++k) {
            cladesFreeLins[k] = clades[k].length;
            for( int l : taxaPartialOrder[k] ) {
                cladesFreeLins[k] -= clades[l].length;
            }
            assert cladesFreeLins[k] >= 0;
        }

        this.taxaPartialOrder = taxaPartialOrder;
        iters = new LinsIterator[clades.length+1];
        vals = new int[iters.length][];

        int nMax = 0;
        for(boolean b : maximal) {
           nMax += b ? 1 : 0;
        }
        maximalClades = new int[nMax];
        nFreeLineages = externalNodeCount;

        nMax = 0;
        for(int m = 0; m < maximal.length; ++m) {
            if( maximal[m] ) {
              maximalClades[nMax] = m;
              ++nMax;
              nFreeLineages -= clades[m].length;
            }
        }

        assert nFreeLineages >= 0;
    }

    int setup(int[] ranks) {
        int n = cladesFreeLins.length;

        nCurIters = 0;

        for(int k = 0; k < n; ++k) {
            setOneIterator(ranks, taxaPartialOrder[k], cladesFreeLins[k], ranks[k]);
        }

        if( nFreeLineages > 0 ) {
          setOneIterator(ranks, maximalClades, nFreeLineages, n+1);
        }
        
        for(int k = 0; k < nCurIters-1; ++k) {
            vals[k] = iters[k].next();
        }

        return nCurIters;
    }

    private void setOneIterator(int[] ranks, int[] joinerClades, int nl, int rank) {
        int nSubs = joinerClades.length;

        LinsIterator itr = null;
        if( nSubs == 0 ) {
            itr = new LinsIterator(nl, rank, null);
        } else /*if( nl > 0 || nSubs > 2 ) */ {
            final int[] s = new int[nSubs];
            for(int i = 0; i < nSubs; ++i) {
                s[i] = ranks[joinerClades[i]];
            }
            itr = new LinsIterator(nl, rank, s);
        }

        if( itr != null ) {
            // sorted according to rank
            iters[itr.rank-1] = itr;
            itr.startIter();
            ++nCurIters;
        }
    }

    int[][] next()
    {
        int[] l = iters[nCurIters-1].next();

        if( l != null ) {
            vals[nCurIters-1] = l;
            return vals;
        }

        int i = nCurIters-2;
        for( ; i >= 0; --i) {
            if( (vals[i] = iters[i].next()) != null) {
                break;
            }
        }

        if( i < 0 ) {
            return null;
        }

        ++i;

        for( ; i < nCurIters; ++i) {
            iters[i].startIter();
            vals[i] = iters[i].next();
        }

        return vals;
    }

    public int[][] allJoiners() {
        int[][] joiners = new int[nCurIters][];

        for(int i = 0; i < nCurIters; ++i) {
            joiners[i] = iters[i].ljoins();
        }
        return joiners;
    }

    public int nStart(int i) {
        return iters[i].nStart;
    }

    class LinsIterator {

        private final int rank;
        private final int nStart;
        private final int[] joiners;
        private final int[] aStart;
        private final int[] lins;
        private int lastJoinger;
        private boolean stopIter;

        LinsIterator(int ns, int r, int[] jnr) {
            rank = r;
            nStart = ns;
            joiners = new int [r];

            lastJoinger = -1;

            // 2 for start+end, rank-1 intermediate levels
            aStart = new int [2 + rank-1];
            lins = new int [2 + rank-1];

            for(int k = 0; k < rank; ++k) {
                joiners[k] = 0;
            }

            if( jnr != null ) {
                for (int j : jnr) {
                    joiners[j] = 1;
                    if (lastJoinger < j) {
                        lastJoinger = j;
                    }
                }
            }
            aStart[0] = ns;

            if( lastJoinger <= 0 ) {
                for(int i = 1; i < rank+1; ++i) {
                    aStart[i] = 2;
                }
                if( rank > 1) {
                    // first iteration increments this
                    aStart[rank-1] -= 1;
                }
            } else {
                //assert(rank > 1);

                if( nStart > 0 ) {
                    int i = 1;
                    for(; i < lastJoinger+1; ++i) {
                        aStart[i] = 1;
                    }
                    for(; i < rank+1; ++i) {
                        aStart[i] = 2;
                    }
                } else {
                    int mj = jnr[0];
                    for(int k = 0; k < jnr.length; ++k) {
                        mj = Math.min(mj, jnr[k]);
                    }
                    int i = 1;
                    for(; i < mj+1; ++i) {
                        aStart[i] = 0;
                    }
                    for(; i < lastJoinger+1; ++i) {
                        aStart[i] = 1;
                    }
                    for(; i < rank+1; ++i) {
                        aStart[i] = 2;
                    }
                }
                // first iteration increments this
                aStart[rank-1] -= 1;
            }

        }

        void startIter() {
            for(int i = 0; i < rank+1; ++i) {
                lins[i] = aStart[i];
            }
            stopIter = false;
        }

        final int[] next()
        {
            int i = rank - 1;
            if( lastJoinger <= 0 ) {
                while( i >= 1 && lins[i] == lins[i-1]) {
                    --i;
                }
                if( i == 0 ) {
                    if( rank == 1 ) {
                        if( !stopIter ) {
                            stopIter = true;
                            return lins;
                        }
                    }
                    return null;
                }
                lins[i] += 1;
                ++i;
                while( i < rank ) {
                    lins[i] = 2;
                    ++i;
                }
            } else {

                while( i >= 1 && lins[i] == lins[i-1] + joiners[i-1] ) {
                    --i;
                }
                if( i == 0 ) {
                    return null;
                }
                lins[i] += 1;
                i++;
                while( i < rank ) {
                    lins[i] = (i <= lastJoinger) ? 1 : 2;
                    i++;
                }
            }
            return lins;
        }

        final int[] ljoins()  {
            return joiners;
        }
    }
}


