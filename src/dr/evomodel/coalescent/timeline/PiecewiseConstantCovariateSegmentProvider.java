/*
 * PiecewiseConstantCovariateSegmentProvider.java
 *
 * Copyright (c) 2002-2026 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodel.coalescent.timeline;

import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;

import java.util.Arrays;

/**
 * Splits each tree interval at a fixed set of covariate breakpoints and provides
 * a per-segment covariate vector drawn from a user-supplied array.
 *
 * Segment convention (half-open): segment {@code k} covers {@code [break_k, break_{k+1})}.
 * A breakpoint or event at exactly {@code break_k} belongs to segment {@code k}
 * (the right-hand value), consistent with the survival sub-interval immediately
 * to the left having used the segment-{@code (k-1)} value.
 *
 * This convention is applied identically in {@link #forEachSegment} and
 * {@link #getEventContext}; using different conventions in the two methods would
 * produce a gradient inconsistent with the log-likelihood.
 *
 * @author Filippo Monti
 */
public final class PiecewiseConstantCovariateSegmentProvider extends AbstractModel
        implements CoalescentSegmentProvider {

    /** Sorted breakpoint times (ascending).  Segment k covers [breaks[k], breaks[k+1]). */
    private final double[] breaks;

    /**
     * Covariate vectors per segment.  {@code covariates[k]} is the covariate vector on
     * {@code [breaks[k], breaks[k+1])}.  Length must equal {@code breaks.length - 1}.
     */
    private final double[][] covariates;

    /** Reused mutable context; {@link #forEachSegment} calls set() before each consumer call. */
    private final BasisContext mutableContext;

    /**
     * @param breaks     sorted breakpoint times; must have at least 2 entries
     * @param covariates covariate vectors per segment; length must equal breaks.length - 1
     */
    public PiecewiseConstantCovariateSegmentProvider(double[] breaks, double[][] covariates) {
        super("PiecewiseConstantCovariateSegmentProvider");

        if (breaks.length < 2)
            throw new IllegalArgumentException("breaks must have at least 2 entries");
        if (covariates.length != breaks.length - 1)
            throw new IllegalArgumentException(
                    "covariates.length must equal breaks.length - 1");
        if (covariates.length > 0) {
            int dim = covariates[0].length;
            for (double[] cv : covariates)
                if (cv.length != dim)
                    throw new IllegalArgumentException("all covariate vectors must have the same dimension");
        }

        this.breaks    = breaks.clone();
        this.covariates = covariates.clone();
        for (int i = 0; i < covariates.length; i++)
            this.covariates[i] = covariates[i].clone();

        this.mutableContext = new BasisContext(-1, -1, covariates.length > 0
                ? new double[covariates[0].length] : new double[0]);
    }

    /**
     * Returns the segment index for time {@code t}: the largest {@code k} such that
     * {@code breaks[k] <= t}.  Falls back to 0 for {@code t < breaks[0]} and to
     * {@code nSegments - 1} for {@code t >= breaks[nBreaks - 1]}.
     */
    private int segmentIndexAt(double t) {
        // Binary search for the rightmost break <= t.
        int lo = 0, hi = breaks.length - 2; // hi = nSegments - 1
        if (t < breaks[lo])   return 0;
        if (t >= breaks[hi + 1]) return hi;
        // Standard lower_bound to find first break > t, then step back.
        int pos = Arrays.binarySearch(breaks, t);
        if (pos >= 0) {
            // exact match: t == breaks[pos]; the half-open convention says segment = pos
            // (segment k starts at breaks[k]), but we must not exceed nSegments-1.
            return Math.min(pos, hi);
        } else {
            // insertion point ip = -(pos+1), first break > t: segment = ip - 1
            int ip = -(pos + 1);
            return ip - 1;
        }
    }

    @Override
    public void forEachSegment(int treeIndex, double start, double end,
                               SegmentConsumer consumer) {
        int nSeg = breaks.length - 1;

        // First segment index that starts at or before 'start'.
        int seg = segmentIndexAt(start);

        double segStart = start;
        while (segStart < end) {
            // Current segment covers [breaks[seg], breaks[seg+1]).
            double segEnd = (seg < nSeg - 1) ? Math.min(breaks[seg + 1], end) : end;
            if (segEnd > segStart) {
                mutableContext.set(treeIndex, seg, covariates[seg]);
                consumer.accept(segStart, segEnd, mutableContext);
            }
            if (seg >= nSeg - 1) break;
            seg++;
            segStart = breaks[seg];
        }
    }

    @Override
    public BasisContext getEventContext(int treeIndex, double t) {
        int seg = segmentIndexAt(t);
        mutableContext.set(treeIndex, seg, covariates[seg]);
        return mutableContext;
    }

    /** Number of covariate dimensions. */
    public int getCovariateDimension() {
        return covariates.length > 0 ? covariates[0].length : 0;
    }

    /** Number of segments. */
    public int getSegmentCount() {
        return breaks.length - 1;
    }

    @Override protected void handleModelChangedEvent(Model m, Object o, int i) { fireModelChanged(); }
    @Override protected void handleVariableChangedEvent(Variable v, int i, Parameter.ChangeType t) { fireModelChanged(); }
    @Override protected void storeState() {}
    @Override protected void restoreState() {}
    @Override protected void acceptState() {}
}
