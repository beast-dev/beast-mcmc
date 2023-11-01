/*
 * RandomField.java
 *
 * Copyright (c) 2002-2023 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

import dr.inference.model.*;
import dr.math.distributions.RandomFieldDistribution;

/**
 * @author Marc Suchard
 * @author Pratyusa Datta
 * @author Filippo Monti
 * @author Xiang Ji
 */

public class RandomField extends AbstractModelLikelihood {

    private final Parameter field;
    private final RandomFieldDistribution distribution;

    private boolean likelihoodKnown;
    private double logLikelihood;

    private boolean savedLikelihoodKnown;
    private double savedLogLikelihood;

    public RandomField(String name,
                       Parameter field,
                       RandomFieldDistribution distribution) {
        super(name);

        this.field = field;
        this.distribution = distribution;

        addVariable(field);
        addModel(distribution);

        likelihoodKnown = false;
    }

    public Parameter getField() { return field; }

    public RandomFieldDistribution getDistribution() { return distribution; }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
        if (model == distribution) {
            throw new RuntimeException("Not yet implemented");
        } else {
            throw new IllegalArgumentException("Unknown model");
        }
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        if (variable == field) {
            likelihoodKnown = false;
        } else {
            throw new IllegalArgumentException("Unknown variable");
        }
    }

    @Override
    protected void storeState() {
        savedLikelihoodKnown = likelihoodKnown;
        savedLogLikelihood = logLikelihood;
    }

    @Override
    protected void restoreState() {
        likelihoodKnown = savedLikelihoodKnown;
        logLikelihood = savedLogLikelihood;
    }

    @Override
    protected void acceptState() { }

    @Override
    public Model getModel() { return this; }

    @Override
    public double getLogLikelihood() {

        if (!likelihoodKnown) {
            logLikelihood = 0;
            likelihoodKnown = true;
        }

        return logLikelihood;
    }

    @Override
    public void makeDirty() {
        likelihoodKnown = false;
    }
}