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
 * This class gives an SIR trajectory and hands off a rate of coalescence at a given point in time.
 * SIR trajectories are simulated deterministically
 * Transmission is frequency-dependent
 *
 * transmissionRateParameter is the number of contacts an infected individual makes per time unit
 * recoveryRateParameter is the number of recovery events an infected individual makes per time unit
 * susceptibleParameter is the number of susceptible hosts at present day (t=0)
 * infectedParameter is the number of infected hosts at present day (t=0)
 * recoveredParameter is the number of recovered hosts at present day (t=0)
 *
 * @author Trevor Bedford
 * @author Tanja Stadler
 * @author Denise Kuehnert
 * @author David Rasmussen
 * @author Sam Lycett
 * @author Erik Volz
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
                    Parameter susceptiblesParameter,
                    Parameter infectedsParameter,
                    Parameter recoveredsParameter,
                    Type units) {

        this(SIRModelParser.SIR_MODEL, transmissionRateParameter, recoveryRateParameter,
            susceptiblesParameter, infectedsParameter, recoveredsParameter, units);
    }

    /**
     * Construct demographic model with default settings
     */
    public SIRModel(String name,
                    Parameter transmissionRateParameter,
                    Parameter recoveryRateParameter,
                    Parameter susceptiblesParameter,
                    Parameter infectedsParameter,
                    Parameter recoveredsParameter,
                    Type units) {

        super(name);

        demographicFunction = new SIRDemographicFunction(units);

        this.transmissionRateParameter = transmissionRateParameter;
        addVariable(transmissionRateParameter);
        transmissionRateParameter.addBounds(new Parameter.DefaultBounds(Double.MAX_VALUE, 0.0, 1));

        this.recoveryRateParameter = recoveryRateParameter;
        addVariable(recoveryRateParameter);
        recoveryRateParameter.addBounds(new Parameter.DefaultBounds(Double.MAX_VALUE, 0.0, 1));

        this.susceptiblesParameter = susceptiblesParameter;
        addVariable(susceptiblesParameter);
        susceptiblesParameter.addBounds(new Parameter.DefaultBounds(Double.MAX_VALUE, 0.0, 1));

        this.infectedsParameter = infectedsParameter;
        addVariable(infectedsParameter);
        infectedsParameter.addBounds(new Parameter.DefaultBounds(Double.MAX_VALUE, 0.0, 1));

        this.recoveredsParameter = recoveredsParameter;
        addVariable(recoveredsParameter);
        recoveredsParameter.addBounds(new Parameter.DefaultBounds(Double.MAX_VALUE, 0.0, 1));

        setUnits(units);
    }

    @Override
    protected void handleVariableChangedEvent(final Variable variable, final int index, final Variable.ChangeType type) {
        demographicFunction.setTransmissionRate(transmissionRateParameter.getParameterValue(0));
        demographicFunction.setRecoveryRate(recoveryRateParameter.getParameterValue(0));
        demographicFunction.setSusceptibles(susceptiblesParameter.getParameterValue(0));
        demographicFunction.setInfecteds(infectedsParameter.getParameterValue(0));
        demographicFunction.setRecovereds(recoveredsParameter.getParameterValue(0));
    }

    public DemographicFunction getDemographicFunction() {
        return demographicFunction;
    }

    Parameter transmissionRateParameter = null;
    Parameter recoveryRateParameter = null;
    Parameter susceptiblesParameter = null;
    Parameter infectedsParameter = null;
    Parameter recoveredsParameter = null;

    boolean functionKnown = false;

    SIRDemographicFunction demographicFunction = null;

    class SIRDemographicFunction extends DemographicFunction.Abstract {
        /**
         * Construct demographic model with default settings
         */
        public SIRDemographicFunction(Type units) {
            super(units);
        }

        double transmissionRate = 1.0;
        double recoveryRate = 1.0;
        double susceptibles = 1.0;
        double infecteds = 1.0;
        double recovereds = 1.0;

        public void setTransmissionRate(final double transmissionRate) {
            this.transmissionRate = transmissionRate;
            functionKnown = false;
        }

        public void setRecoveryRate(final double recoveryRate) {
            this.recoveryRate = recoveryRate;
            functionKnown = false;
        }

        public void setSusceptibles(final double susceptibles) {
            this.susceptibles = susceptibles;
            functionKnown = false;
        }

        public void setInfecteds(final double infecteds) {
            this.infecteds = infecteds;
            functionKnown = false;
        }

        public void setRecovereds(final double recovereds) {
            this.recovereds = recovereds;
            functionKnown = false;
        }

        public double getDemographic(final double t) {
            if (!functionKnown) {
                calculateDemographicFunction();
                functionKnown = true;
            }
            return 1.0;
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
            return 1.0;
        }

        public double getInverseIntensity(final double x) {
            return 1.0;
        }

        public double getIntegral(final double start, final double finish) {
            return 1.0;
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
