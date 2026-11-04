/*
 * ReciprocalLikelihood.java
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

package dr.inference.model;

import dr.inference.loggers.LogColumn;
import dr.xml.Reportable;

import java.util.Set;

/**
 * @author Marc A. Suchard
 */
public class ReciprocalLikelihood implements Likelihood, Reportable {

    private final Likelihood likelihood;

    public ReciprocalLikelihood(Likelihood likelihood) {
        this.likelihood = likelihood;
    }

    @Override
    public LogColumn[] getColumns() {
        return likelihood.getColumns();
    }

    @Override
    public Model getModel() {
        return likelihood.getModel();
    }

    @Override
    public double getLogLikelihood() {
        return -likelihood.getLogLikelihood();
    }

    @Override
    public void makeDirty() {
        likelihood.makeDirty();
    }

    @Override
    public String prettyName() {
        return likelihood.prettyName();
    }

    @Override
    public Set<Likelihood> getLikelihoodSet() {
        return likelihood.getLikelihoodSet();
    }

    @Override
    public boolean isUsed() {
        return likelihood.isUsed();
    }

    @Override
    public void setUsed() {
        likelihood.setUsed();
    }

    @Override
    public boolean evaluateEarly() {
        return likelihood.evaluateEarly();
    }

    @Override
    public String getId() {
        return likelihood.getId();
    }

    @Override
    public void setId(String id) {
        likelihood.setId(id);
    }

    @Override
    public String getReport() {
        return "Reciprocal likelihood " + getId() + " " + getLogLikelihood();
    }
}
