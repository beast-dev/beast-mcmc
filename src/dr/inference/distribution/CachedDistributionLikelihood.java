/*
 * CachedDistributionLikelihood.java
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

package dr.inference.distribution;

import dr.inference.model.AbstractModelLikelihood;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;

/**
 * @author Marc Suchard
 */
public class CachedDistributionLikelihood extends AbstractModelLikelihood {

    public CachedDistributionLikelihood(String name, AbstractDistributionLikelihood likelihood, Variable variable) {
        super(name);
        this.likelihood = likelihood;
        addVariable(variable);
    }

    protected void handleModelChangedEvent(Model model, Object object, int index) {
        likelihoodKnown = false;
    }

    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        likelihoodKnown = false;
    }

    protected void storeState() {
        storedLikelihoodKnow = likelihoodKnown;
        storedLogLikelihood = logLikelihood;

    }

    protected void restoreState() {
        likelihoodKnown = storedLikelihoodKnow;
        logLikelihood = storedLogLikelihood;

    }

    protected void acceptState() {

    }

    public AbstractDistributionLikelihood getDistributionLikelihood() { return likelihood; }

    public Model getModel() {
        return this;
    }

    public double getLogLikelihood() {
        if (!likelihoodKnown) {
            logLikelihood = calculateLogLikelihood();
            likelihoodKnown = true;
        }
        return logLikelihood;
    }

    private double calculateLogLikelihood() {
        return likelihood.calculateLogLikelihood();
    }

    public void makeDirty() {
        likelihoodKnown = false;
    }

    private boolean likelihoodKnown;
    private boolean storedLikelihoodKnow;
    private double logLikelihood;
    private double storedLogLikelihood;

    private final AbstractDistributionLikelihood likelihood;
}
