/*
 * BasisContext.java
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

/**
 * Mutable context object passed by {@link CoalescentSegmentProvider} to
 * {@link SegmentConsumer#accept} on each segment.
 *
 * Callers MUST NOT retain a reference to this object after returning from
 * {@link SegmentConsumer#accept}; the provider reuses the same instance and
 * calls {@link #set} before each consumer invocation.
 *
 * @author Filippo Monti
 */
public final class BasisContext {

    public int treeIndex;
    public int segmentIndex;
    /** Covariate values constant on this segment; never null; empty for no-covariate case. */
    public double[] covariates;

    public BasisContext(int treeIndex, int segmentIndex, double[] covariates) {
        set(treeIndex, segmentIndex, covariates);
    }

    public void set(int treeIndex, int segmentIndex, double[] covariates) {
        this.treeIndex    = treeIndex;
        this.segmentIndex = segmentIndex;
        this.covariates   = covariates;
    }

    /** Convenience accessor (same reference as the public field). */
    public double[] getCovariates() { return covariates; }
}
