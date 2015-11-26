/*
 * LinkedParameter.java
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

package dr.app.beauti.components.linkedparameters;

import dr.app.beauti.options.Operator;
import dr.app.beauti.options.Parameter;
import dr.app.beauti.types.PriorType;

import java.util.List;

/**
 * @author Andrew Rambaut
 * @author Marc A. Suchard
 */
public class LinkedParameter {

    private String name;
    final private Parameter argumentParameter;
    final private Operator argumentOperator;
    final private LinkedParameterComponentOptions options;

    /**
     * A simple data class to store the definition of a linked parameter.
     * @param name the name of the parameter
     * @param argumentParameter the control parameter
     * @param argumentOperator the control parameter's operator
     */
    public LinkedParameter(String name, Parameter argumentParameter, Operator argumentOperator, LinkedParameterComponentOptions options) {
        this.name = name;
        this.argumentParameter = argumentParameter;
        this.argumentOperator = argumentOperator;
        this.options = options;
    }

    public Parameter getArgumentParameter() {
        return argumentParameter;
    }

    public Operator getArgumentOperator() {
        return argumentOperator;
    }

    public List<Parameter> getDependentParameterList() {
        return options.getDependentParameters(this);
    }

    public void linkDependentParameters(List<Parameter> parameterList) {
        options.setDependentParameters(this, parameterList);
    }

    public String getName() {
        return name;
    }

    public void setName(String name) {
        argumentParameter.setName(name);
        argumentOperator.setName(name);
        this.name = name;

    }

}
