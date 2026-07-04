/*
 * NoCovariateSegmentProvider.java
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

/**
 * Stateless singleton that passes each tree interval through as a single segment
 * with an empty covariate vector.
 *
 * Use {@code NoCovariateSegmentProvider.INSTANCE} everywhere the likelihood
 * does not require piecewise-constant covariates.
 *
 * @author Filippo Monti
 */
public final class NoCovariateSegmentProvider extends AbstractModel
        implements CoalescentSegmentProvider {

    public static final NoCovariateSegmentProvider INSTANCE = new NoCovariateSegmentProvider();

    private static final BasisContext EMPTY_CONTEXT =
            new BasisContext(-1, -1, new double[0]);

    private NoCovariateSegmentProvider() {
        super("NoCovariateSegmentProvider");
    }

    @Override
    public void forEachSegment(int treeIndex, double start, double end,
                               SegmentConsumer consumer) {
        consumer.accept(start, end, EMPTY_CONTEXT);
    }

    @Override
    public BasisContext getEventContext(int treeIndex, double t) {
        return EMPTY_CONTEXT;
    }

    @Override protected void handleModelChangedEvent(Model m, Object o, int i) {}
    @Override protected void handleVariableChangedEvent(Variable v, int i, Parameter.ChangeType t) {}
    @Override protected void storeState() {}
    @Override protected void restoreState() {}
    @Override protected void acceptState() {}
}
