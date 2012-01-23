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
import dr.evomodel.coalescent.DemographicModel;
import dr.inference.loggers.LogColumn;
import dr.inference.loggers.NumberColumn;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.inference.model.Likelihood;
import dr.inference.model.Model;
import dr.math.distributions.NormalDistribution;


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
 * @author Andrew Rambaut
 */
public class SIRModel extends DemographicModel implements Likelihood {

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

        demographicFunction = new SIRDemographicFunction(units);
        setUnits(units);

    }

    @Override
    protected void handleVariableChangedEvent(final Variable variable, final int index, final Variable.ChangeType type) {
        demographicFunction.reset();
        fireModelChanged();
        likelihoodKnown = false;
    }

    @Override
    protected void storeState() {
        storedLogLikelihood = logLikelihood;
        demographicFunction.store();
    }

    @Override
    protected void restoreState() {
        logLikelihood = storedLogLikelihood;
        likelihoodKnown = true;
        demographicFunction.restore();
    }

    /* Likelihood methods */

    public String prettyName() {
        return Likelihood.Abstract.getPrettyName(this);
    }

    public final double getLogLikelihood() {
        if (!getLikelihoodKnown()) {
            logLikelihood = calculateLogLikelihood();
            likelihoodKnown = true;
        }
        return logLikelihood;
    }

    protected boolean getLikelihoodKnown() {
        return likelihoodKnown;
    }

    public boolean isUsed() {
        return isUsed;
    }

    public void setUsed() {
        isUsed = true;
    }

    public Model getModel() {
        return this;
    }

    public void makeDirty() {
        likelihoodKnown = false;
    }

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

    public double calculateLogLikelihood() {
        double r = demographicFunction.getRecovereds(5);
        return NormalDistribution.logPdf(r, 0, 100);
    }

    private boolean isUsed = false;
    private boolean likelihoodKnown = false;
    private double logLikelihood;
    private double storedLogLikelihood;

    public DemographicFunction getDemographicFunction() {
        return demographicFunction;
    }

    Parameter transmissionRateParameter = null;
    Parameter recoveryRateParameter = null;
    Parameter susceptiblesParameter = null;
    Parameter infectedsParameter = null;
    Parameter recoveredsParameter = null;

    SIRDemographicFunction demographicFunction = null;

    class SIRDemographicFunction extends DemographicFunction.Abstract {

        DynamicalSystem syst = new DynamicalSystem(0.01);

        public SIRDemographicFunction(Type units) {

            super(units);

            syst.addVariable("susceptibles", susceptiblesParameter.getParameterValue(0));
            syst.addVariable("infecteds", infectedsParameter.getParameterValue(0));
            syst.addVariable("recovereds", recoveredsParameter.getParameterValue(0));
            syst.addVariable("total", susceptiblesParameter.getParameterValue(0) + infectedsParameter.getParameterValue(0)
                    + recoveredsParameter.getParameterValue(0));
            syst.addForce("contact", transmissionRateParameter.getParameterValue(0), new String[]{"infecteds","susceptibles"},
                    new String[]{"total"}, "susceptibles", "infecteds");
            syst.addForce("recovery", recoveryRateParameter.getParameterValue(0), new String[]{"infecteds"},
                    new String[]{}, "infecteds", "recovereds");

        }

        public void reset() {
            syst.resetVar("susceptibles", susceptiblesParameter.getParameterValue(0));
            syst.resetVar("infecteds", infectedsParameter.getParameterValue(0));
            syst.resetVar("recovereds", recoveredsParameter.getParameterValue(0));
            syst.resetVar("total", susceptiblesParameter.getParameterValue(0) + infectedsParameter.getParameterValue(0)
                    + recoveredsParameter.getParameterValue(0));
            syst.resetForce("contact", transmissionRateParameter.getParameterValue(0));
            syst.resetForce("recovery", recoveryRateParameter.getParameterValue(0));
            syst.resetTime();
        }

        public void store () {
            syst.store();
        }

        public void restore () {
            syst.restore();
        }

        // return N(t)
        public double getDemographic(final double t) {

            double inf = syst.getValue("infecteds", t);
            if (inf < 1)
                inf = 1.0;

            double numer = inf * syst.getValue("total", t);
            double denom = 2.0 * transmissionRateParameter.getParameterValue(0) * syst.getValue("susceptibles", t);

            return numer / denom;

        }

        // return log N(t)
        public double getLogDemographic(final double t) {
            return Math.log(getDemographic(t));
        }

         // return S(t)
        public double getSusceptibles(final double t) {
            return syst.getValue("susceptibles", t);
        }

         // return I(t)
        public double getInfecteds(final double t) {
            return syst.getValue("infecteds", t);
        }

         // return R(t)
        public double getRecovereds(final double t) {
            return syst.getValue("recovereds", t);
        }

        // return t/N(t)
        public double getIntensity(final double t) {
            return 1.0;
        }

        // return x*N(t)
        public double getInverseIntensity(final double x) {
            return 1.0;
        }

        // return integral of 1/N(t)
        public double getIntegral(final double start, final double finish) {

            double inf = syst.getIntegral("infecteds", start, finish);
            if (inf < 1)
                inf = 1.0;

            double numer = 2.0 * transmissionRateParameter.getParameterValue(0) * syst.getIntegral("susceptibles", start, finish);
            double denom = inf * syst.getIntegral("total", start, finish);

            double integral = (finish-start)*(numer / denom);

            return integral ;

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