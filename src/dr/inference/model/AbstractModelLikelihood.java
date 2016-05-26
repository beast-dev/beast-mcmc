/*
 * AbstractModelLikelihood.java
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

package dr.inference.model;

import dr.inference.loggers.LogColumn;
import dr.inference.loggers.NumberColumn;

import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;

/**
 * @author Joseph Heled
 *         Date: 16/04/2009
 */
public abstract class AbstractModelLikelihood extends AbstractModel implements Likelihood {
    /**
     * @param name Model Name
     */
    public AbstractModelLikelihood(String name) {
        super(name);
    }

    public String prettyName() {
        return Likelihood.Abstract.getPrettyName(this);
    }

    @Override
    public Set<Likelihood> getLikelihoodSet() {
        return new HashSet<Likelihood>(Arrays.asList(this));
    }

    @Override
    public boolean isUsed() {
        return isUsed;
    }

    public void setUsed() {
        isUsed = true;
    }

    public boolean evaluateEarly() {
        return false;
    }

    private boolean isUsed = false;

    // **************************************************************
    // Loggable IMPLEMENTATION
    // **************************************************************

    public LogColumn[] getColumns() {
        return new LogColumn[]{
                new LikelihoodColumn(getId())
        };
    }

    protected class LikelihoodColumn extends NumberColumn {
        public LikelihoodColumn(String label) {
            super(label);
        }

        public double getDoubleValue() {
            return getLogLikelihood();
        }
    }
}
