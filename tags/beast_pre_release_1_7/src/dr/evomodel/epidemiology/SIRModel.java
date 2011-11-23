/*
 * ExponentialGrowthModel.java
 *
 * Copyright (C) 2002-2009 Alexei Drummond and Andrew Rambaut
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
 * BEAST is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package dr.evomodel.epidemiology;

import dr.evolution.coalescent.DemographicFunction;
import dr.evolution.coalescent.ExponentialGrowth;
import dr.evolution.util.Units;
import dr.evomodel.coalescent.DemographicModel;
import dr.evomodelxml.coalescent.ExponentialGrowthModelParser;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import org.apache.commons.math.FunctionEvaluationException;

/**
 * This class models an exponentially growing (or shrinking) population
 * (Parameters: N0=present-day population size; r=growth rate).
 * This model is nested with the constant-population size model (r=0).
 *
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @version $Id: ExponentialGrowthModel.java,v 1.14 2005/05/24 20:25:57 rambaut Exp $
 */
public class SIRModel extends DemographicModel {

    //
    // Public stuff
    //
    /**
     * Construct demographic model with default settings
     */
    public SIRModel(Parameter transmissionRateParameter,
                    Parameter recoveryRateParameter,
                    Type units) {

        this(SIRModelParser.SIR_MODEL, transmissionRateParameter, recoveryRateParameter, units);
    }

    /**
     * Construct demographic model with default settings
     */
    public SIRModel(String name,
                    Parameter transmissionRateParameter,
                    Parameter recoveryRateParameter,
                    Type units) {

        super(name);

        demographicFunction = new SIRDemographicFunction(units);

        this.transmissionRateParameter = transmissionRateParameter;
        addVariable(transmissionRateParameter);
        transmissionRateParameter.addBounds(new Parameter.DefaultBounds(Double.MAX_VALUE, 0.0, 1));

        this.recoveryRateParameter = recoveryRateParameter;
        addVariable(recoveryRateParameter);
        recoveryRateParameter.addBounds(new Parameter.DefaultBounds(Double.MAX_VALUE, 0.0, 1));

        setUnits(units);
    }

    @Override
    protected void handleVariableChangedEvent(final Variable variable, final int index, final Variable.ChangeType type) {
        demographicFunction.setBeta(transmissionRateParameter.getParameterValue(0));
        demographicFunction.setGamma(recoveryRateParameter.getParameterValue(0));
    }

    public DemographicFunction getDemographicFunction() {
        return demographicFunction;
    }

    Parameter transmissionRateParameter = null;
    Parameter recoveryRateParameter = null;

    boolean functionKnown = false;

    SIRDemographicFunction demographicFunction = null;

    class SIRDemographicFunction extends DemographicFunction.Abstract {
        /**
         * Construct demographic model with default settings
         */
        public SIRDemographicFunction(Type units) {
            super(units);
        }

        double beta = 1.0;
        double gamma = 1.0;

        public void setBeta(final double beta) {
            this.beta = beta;
            functionKnown = false;
        }

        public void setGamma(final double gamma) {
            this.gamma = gamma;
            functionKnown = false;
        }

        public double getDemographic(final double t) {
            if (!functionKnown) {
                calculateDemographicFunction();
                functionKnown = true;
            }
            return 0;
        }

        public double getLogDemographic(final double t) {
            if (!functionKnown) {
                calculateDemographicFunction();
                functionKnown = true;
            }
            return 0;
        }

        private void calculateDemographicFunction() {
            // todo
        }

        public double getIntensity(final double t) {
            throw new RuntimeException("Function not implemented");
        }

        public double getInverseIntensity(final double x) {
            throw new RuntimeException("Function not implemented");
        }

        // ignore the rest:
        public int getNumArguments() {
            return 0;
        }

        public String getArgumentName(final int n) {
            return null;
        }

        public double getArgument(final int n) {
            return 0;
        }

        public void setArgument(final int n, final double value) {
        }

        public double getLowerBound(final int n) {
            return 0;
        }

        public double getUpperBound(final int n) {
            return 0;
        }

        public DemographicFunction getCopy() {
            return null;
        }

        public double getThreshold() {
            return 0;
        }

    }
}
