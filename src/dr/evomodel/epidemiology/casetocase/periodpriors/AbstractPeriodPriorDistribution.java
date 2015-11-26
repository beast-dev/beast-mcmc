/*
 * AbstractPeriodPriorDistribution.java
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

package dr.evomodel.epidemiology.casetocase.periodpriors;

import dr.inference.loggers.LogColumn;
import dr.inference.loggers.Loggable;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;

import java.util.ArrayList;


/**
 * Abstract class for the probability of a set of latent or infectious periods being drawn from an unknown probability
 * distribution, given hyperpriors on the parameters of that distribution.
 */

public abstract class AbstractPeriodPriorDistribution extends AbstractModel implements Loggable {

    // are we working on the logarithms of the values?
    protected boolean log;

    protected double logL;
    protected double storedLogL;

    public AbstractPeriodPriorDistribution(String name, boolean log) {
        super(name);
        this.log = log;
    }

    public double getLogLikelihood(double[] values){
        if(!log){
            return calculateLogLikelihood(values);
        } else {
            double[] logValues = new double[values.length];
            for(int i=0; i<values.length; i++){
                logValues[i] = Math.log(values[i]);
            }
            return calculateLogLikelihood(logValues);
        }
    }

    public double getLogPosteriorProbability(double newValue, double minValue){
        if(!log){
            return calculateLogPosteriorProbability(newValue, minValue);
        } else {
            return calculateLogPosteriorProbability(Math.log(newValue), Math.log(minValue));
        }
    }

    public double getLogPosteriorCDF(double limit, boolean upper){
        if(!log){
            return calculateLogPosteriorCDF(limit, upper);
        } else {
            return calculateLogPosteriorCDF(Math.log(limit), upper);
        }
    }

    protected void handleModelChangedEvent(Model model, Object object, int index) {
        //generally nothing to do
    }

    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        //generally nothing to do
    }

    protected void storeState() {
        storedLogL = logL;

    }

    protected void restoreState() {
        logL = storedLogL;
    }

    protected void acceptState() {
        //generally nothing to do
    }


    public LogColumn[] getColumns() {
        ArrayList<LogColumn> columns = new ArrayList<LogColumn>();

        columns.add(new LogColumn.Abstract(getModelName()+"_LL"){
            protected String getFormattedValue() {
                return String.valueOf(logL);
            }
        });
        return columns.toArray(new LogColumn[columns.size()]);
    }

    public abstract void reset();

    public abstract double calculateLogPosteriorProbability(double newValue, double minValue);

    public abstract double calculateLogPosteriorCDF(double limit, boolean upper);

    public abstract double calculateLogLikelihood(double[] values);




}
