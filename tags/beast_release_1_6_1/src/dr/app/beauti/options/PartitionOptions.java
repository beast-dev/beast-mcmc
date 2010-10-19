/*
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

package dr.app.beauti.options;

import dr.app.beauti.enumTypes.PriorScaleType;
import dr.app.beauti.enumTypes.PriorType;

import java.util.List;


/**
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @author Walter Xie
 */
public abstract class PartitionOptions extends ModelOptions {

    protected String partitionName;

	protected abstract void selectParameters(List<Parameter> params);
    protected abstract void selectOperators(List<Operator> ops);

    public abstract String getPrefix();    

    protected void createParameterClockRateUndefinedPrior(PartitionOptions options, String name, String description, PriorScaleType scaleType,
            double initial, double lower, double upper) {
        new Parameter.Builder(name, description).scaleType(scaleType).prior(PriorType.UNDEFINED).initial(initial)
                .lower(lower).upper(upper).partitionOptions(options).build(parameters);
    }

    protected void createParameterClockRateUniform(PartitionOptions options, String name, String description, PriorScaleType scaleType,
            double initial, double lower, double upper) {
        new Parameter.Builder(name, description).scaleType(scaleType).prior(PriorType.UNIFORM_PRIOR).initial(initial)
                .lower(lower).upper(upper).partitionOptions(options).build(parameters);
    }

    protected void createParameterClockRateGamma(PartitionOptions options, String name, String description, PriorScaleType scaleType,
            double initial, double shape, double scale, double lower, double upper) {
        new Parameter.Builder(name, description).scaleType(scaleType).prior(PriorType.GAMMA_PRIOR).initial(initial).
                shape(shape).scale(scale).lower(lower).upper(upper).partitionOptions(options).build(parameters);
    }

    public void createParameterClockRateExponential(PartitionOptions options, String name, String description, PriorScaleType scaleType, 
            double initial, double mean, double offset, double lower, double upper) {
        new Parameter.Builder(name, description).scaleType(scaleType).prior(PriorType.EXPONENTIAL_PRIOR)
                  .initial(initial).mean(mean).offset(offset).lower(lower).upper(upper).partitionOptions(options).build(parameters);
    }


    protected void createParameterTree(PartitionOptions options, String name, String description, boolean isNodeHeight, double value,
            double lower, double upper) {
        new Parameter.Builder(name, description).isNodeHeight(isNodeHeight).scaleType(PriorScaleType.TIME_SCALE)
                .initial(value).lower(lower).upper(upper).partitionOptions(options).build(parameters);        
    }

    protected void createAllMusParameter(PartitionOptions options, String name, String description) {
        new Parameter.Builder(name, description).partitionOptions(options).build(parameters);
    }

    public Parameter getParameter(String name) {

        Parameter parameter = parameters.get(name);

        if (parameter == null) {
            throw new IllegalArgumentException("Parameter with name, " + name + ", is unknown");
        }

        parameter.setPrefix(getPrefix());

        return parameter;
    }

    public Operator getOperator(String name) {

        Operator operator = operators.get(name);

        if (operator == null) throw new IllegalArgumentException("Operator with name, " + name + ", is unknown");

        operator.setPrefix(getPrefix());

        return operator;
    }
 
    public String getName() {
        return partitionName;
    }

    public void setName(String name) {
        this.partitionName = name;
    }

    public String toString() {
        return getName();
    }
}
