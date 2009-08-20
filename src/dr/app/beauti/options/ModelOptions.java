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

import dr.app.beauti.enumTypes.OperatorType;
import dr.app.beauti.enumTypes.PriorScaleType;
import dr.app.beauti.enumTypes.PriorType;
import dr.evolution.util.TaxonList;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

/**
 * @author Alexei Drummond
 * @author Andrew Rambaut
 */
public abstract class ModelOptions {

    HashMap<String, Parameter> parameters = new HashMap<String, Parameter>();    
	HashMap<TaxonList, Parameter> statistics = new HashMap<TaxonList, Parameter>();
    HashMap<String, Operator> operators = new HashMap<String, Operator>();

    public static final String version = "1.5";
    public static final int YEARS = 0;
    public static final int MONTHS = 1;
    public static final int DAYS = 2;
    public static final int FORWARDS = 0;
    public static final int BACKWARDS = 1;
    
    public static final int BIN_SIMPLE = 0;
    public static final int BIN_COVARION = 1;

    public static final int GROWTH_RATE = 0;
    public static final int DOUBLING_TIME = 1;
    public static final int CONSTANT_SKYLINE = 0;
    public static final int LINEAR_SKYLINE = 1;

    public static final int SKYRIDE_UNIFORM_SMOOTHING = 0;
    public static final int SKYRIDE_TIME_AWARE_SMOOTHING = 1;

    public static final double demoTuning = 0.75;
    public static final double demoWeights = 3.0;

	protected static final double branchWeights = 30.0;
	protected static final double treeWeights = 15.0;
	protected static final double rateWeights = 3.0;
	
	private final List<ComponentOptions> components = new ArrayList<ComponentOptions>();
	   
	
	//+++++++++++++++++++ Create Operator ++++++++++++++++++++++++++++++++
    public void createOperator(String parameterName, OperatorType type, double tuning, double weight) {
        Parameter parameter = getParameter(parameterName);
        operators.put(parameterName, new Operator(parameterName, "", parameter, type, tuning, weight));
    }

    public void createOperator(String key, String name, String description, String parameterName, OperatorType type, double tuning, double weight) {
        Parameter parameter = getParameter(parameterName);
        operators.put(key, new Operator(name, description, parameter, type, tuning, weight));
    }

    public void createOperator(String key, String name, String description, Parameter parameter1, Parameter parameter2, OperatorType type, double tuning, double weight) {
//        Parameter parameter1 = getParameter(parameterName1);
//        Parameter parameter2 = getParameter(parameterName2);
        operators.put(key, new Operator(name, description, parameter1, parameter2, type, tuning, weight));
    }
    
    public void createTagOperator(String key, String name, String description, String parameterName, String tag, String idref, 
    		OperatorType type, double tuning, double weight) {
    	Parameter parameter = getParameter(parameterName);
      operators.put(key, new Operator(name, description, parameter, tag, idref, type, tuning, weight));
  }

    public void createScaleOperator(String parameterName, double tuning, double weight) {
        Parameter parameter = getParameter(parameterName);
        operators.put(parameterName, new Operator(parameterName, "", parameter, OperatorType.SCALE, tuning, weight));
    }

    public void createScaleAllOperator(String parameterName, double weight) {
        Parameter parameter = getParameter(parameterName);
        operators.put(parameterName, new Operator(parameterName, "", parameter, OperatorType.SCALE_ALL, 0.75, weight));
    }
    
    public void createUpDownAllOperator(String paraName, String opName, String description, double tuning, double weight) {
    	final Parameter parameter = new Parameter(paraName, description);
        operators.put(paraName, new Operator(opName, description, parameter, OperatorType.UP_DOWN_ALL_RATES_HEIGHTS, tuning, weight));
    }

    //+++++++++++++++++++ Create Parameter ++++++++++++++++++++++++++++++++
    public Parameter createParameter(String name, String description) {
        final Parameter parameter = new Parameter(name, description);
        parameters.put(name, parameter);
        return parameter;
    }

    public Parameter createParameter(String name, String description, PriorScaleType scale, double value, double lower, double upper) {
        final Parameter parameter = new Parameter(name, description, scale, value, lower, upper);
        parameters.put(name, parameter);
        return parameter;
    }

    public void createParameter(String name, String description, boolean isNodeHeight, double value, double lower, double upper) {
        parameters.put(name, new Parameter(name, description, isNodeHeight, value, lower, upper));
    }

    public void createScaleParameter(String name, String description, PriorScaleType scale, double value, double lower, double upper) {
        Parameter p = createParameter(name, description, scale, value, lower, upper);
        p.priorType = PriorType.JEFFREYS_PRIOR;
    }

  //+++++++++++++++++++ Create Statistic ++++++++++++++++++++++++++++++++
    public Parameter createStatistic(String name, String description, boolean isDiscrete) {
        final Parameter parameter = new Parameter(name, description, isDiscrete);
        parameters.put(name, parameter);
        return parameter;
    }

    public void createStatistic(String name, String description, double lower, double upper) {
        parameters.put(name, new Parameter(name, description, lower, upper));
    }

    //+++++++++++++++++++ Methods ++++++++++++++++++++++++++++++++
    public Parameter getParameter(String name) {
        Parameter parameter = parameters.get(name);
        if (parameter == null) {
            for (String key : parameters.keySet()) {
                System.err.println(key);
            }
            throw new IllegalArgumentException("Parameter with name, " + name + ", is unknown");
        }
        return parameter;
    }

    public Parameter getStatistic(TaxonList taxonList) {
        Parameter parameter = statistics.get(taxonList);
        if (parameter == null) {
            for (TaxonList key : statistics.keySet()) {
                System.err.println("Taxon list: " + key.getId());
            }
            throw new IllegalArgumentException("Statistic for taxon list, " + taxonList.getId() + ", is unknown");
        }
        return parameter;
    }

    public Operator getOperator(String name) {
        Operator operator = operators.get(name);
        if (operator == null) throw new IllegalArgumentException("Operator with name, " + name + ", is unknown");
        return operator;
    }

    abstract public String getPrefix();

    protected void addComponent(ComponentOptions component) {
        components.add(component);
        component.createParameters(this);
    }

    public ComponentOptions getComponentOptions(Class<?> theClass) {
        for (ComponentOptions component : components) {
            if (theClass.isAssignableFrom(component.getClass())) {
                return component;
            }
        }

        return null;
    }

    protected void selectComponentParameters(ModelOptions options, List<Parameter> params) {
        for (ComponentOptions component : components) {
            component.selectParameters(options, params);
        }
    }

    protected void selectComponentStatistics(ModelOptions options, List<Parameter> stats) {
        for (ComponentOptions component : components) {
            component.selectStatistics(options, stats);
        }
    }

    protected void selectComponentOperators(ModelOptions options, List<Operator> ops) {
        for (ComponentOptions component : components) {
            component.selectOperators(options, ops);
        }
    }

    public HashMap<String, Parameter> getParameters() {
		return parameters;
	}

	public HashMap<TaxonList, Parameter> getStatistics() {
		return statistics;
	}

	public HashMap<String, Operator> getOperators() {
		return operators;
	}


}
