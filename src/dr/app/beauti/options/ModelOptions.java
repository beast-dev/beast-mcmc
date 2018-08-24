/*
 * ModelOptions.java
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

package dr.app.beauti.options;

import dr.app.beauti.components.ComponentFactory;
import dr.app.beauti.types.OperatorType;
import dr.app.beauti.types.PriorScaleType;
import dr.app.beauti.types.PriorType;
import dr.evolution.util.TaxonList;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * @author Alexei Drummond
 * @author Andrew Rambaut
 */
public abstract class ModelOptions implements Serializable {

    private static final long serialVersionUID = 6199011531067286245L;

    protected final Map<String, Parameter> parameters = new HashMap<String, Parameter>();
    protected final Map<String, Operator> operators = new HashMap<String, Operator>();
    protected final Map<TaxonList, Parameter> statistics = new HashMap<TaxonList, Parameter>();

    public static final double demoTuning = 0.75;
    public static final double demoWeights = 3.0;

    protected static final double branchWeights = 30.0;
    protected static final double treeWeights = 30.0;
    protected static final double rateWeights = 3.0;

    private final List<ComponentOptions> components = new ArrayList<ComponentOptions>();

    //+++++++++++++++++++ Abstract Methods ++++++++++++++++++++++++++++++++

    public abstract void initModelParametersAndOpererators();

    public abstract List<Parameter> selectParameters(List<Parameter> params);

    public abstract List<Operator> selectOperators(List<Operator> ops);

    public abstract String getPrefix();

    //+++++++++++++++++++ Create Parameter ++++++++++++++++++++++++++++++++
    public Parameter createParameter(String name, String description) {
        return new Parameter.Builder(name, description).build(parameters);
    }

    public Parameter createParameter(String name, String description, double initial) {
        return new Parameter.Builder(name, description).initial(initial).isFixed(true).build(parameters);
    }

    public Parameter createZeroOneParameter(String name, String description, double initial) {
        return new Parameter.Builder(name, description)
                .initial(initial).isZeroOne(true).build(parameters);
    }

    public Parameter createZeroOneParameterUniformPrior(String name, String description, double initial) {
        return createZeroOneParameterUniformPrior(name, description, initial, false);
    }

    public Parameter createZeroOneParameterUniformPrior(String name, String description, double initial, boolean amtk) {
        return new Parameter.Builder(name, description).prior(PriorType.UNIFORM_PRIOR)
                .initial(initial).isZeroOne(true).isAdaptiveMultivariateCompatible(amtk).build(parameters);
    }

    public Parameter createNonNegativeParameterDirichletPrior(String name, String description, PartitionOptions options,
                                                              double maintainedSum, boolean amtk) {
        return new Parameter.Builder(name, description).prior(PriorType.DIRICHLET_PRIOR)
                .isNonNegative(true).maintainedSum(maintainedSum)
                .isAdaptiveMultivariateCompatible(amtk)
                .partitionOptions(options).build(parameters);
    }

    public Parameter createNonNegativeParameterDirichletPrior(String name, String description, PartitionOptions options,
                                                              int dimension, double maintainedSum, boolean amtk) {
        return createNonNegativeParameterDirichletPrior(name, description, options, dimension, maintainedSum, false, amtk);
    }

    public Parameter createNonNegativeParameterDirichletPrior(String name, String description, PartitionOptions options,
                                                              int dimension, double maintainedSum, boolean priorFixed, boolean amtk) {
        return new Parameter.Builder(name, description).dimension(dimension).prior(PriorType.DIRICHLET_PRIOR)
                .isNonNegative(true).maintainedSum(maintainedSum).isPriorFixed(priorFixed)
                .isAdaptiveMultivariateCompatible(amtk)
                .partitionOptions(options).build(parameters);
    }

    public Parameter createNonNegativeParameterInfinitePrior(String name, String description, PriorScaleType scaleType, double initial) {
        return createNonNegativeParameterInfinitePrior(name, description, null, scaleType, initial, false);
    }

    public Parameter createNonNegativeParameterInfinitePrior(String name, String description, PriorScaleType scaleType, double initial, boolean amtk) {
        return createNonNegativeParameterInfinitePrior(name, description, null, scaleType, initial, amtk);
    }

    public Parameter createNonNegativeParameterInfinitePrior(String name, String description, PartitionOptions options, PriorScaleType scaleType, double initial, boolean amtk) {
        return new Parameter.Builder(name, description).scaleType(scaleType).prior(PriorType.NONE_IMPROPER).isNonNegative(true)
                .partitionOptions(options).initial(initial).isAdaptiveMultivariateCompatible(amtk).build(parameters);
    }

    public Parameter createNonNegativeParameterUniformPrior(String name, String description, PriorScaleType scaleType, double initial,
                                                            double uniformLower, double uniformUpper) {
        return new Parameter.Builder(name, description).scaleType(scaleType).prior(PriorType.UNIFORM_PRIOR).isNonNegative(true)
                .initial(initial).uniformLower(uniformLower).uniformUpper(uniformUpper).build(parameters);
    }

    public Parameter createParameterUniformPrior(String name, String description, PriorScaleType scaleType, double initial,
                                                 double uniformLower, double uniformUpper) {
        return new Parameter.Builder(name, description).scaleType(scaleType).prior(PriorType.UNIFORM_PRIOR)
                .initial(initial).uniformLower(uniformLower).uniformUpper(uniformUpper).build(parameters);
    }

    public Parameter createParameterGammaPrior(String name, String description, PriorScaleType scaleType, double initial,
                                               double shape, double scale, boolean priorFixed) {
        return createParameterGammaPrior(name, description, scaleType, initial, shape, scale, priorFixed, false);
    }

    public Parameter createParameterGammaPrior(String name, String description, PriorScaleType scaleType, double initial,
                                               double shape, double scale, boolean priorFixed, boolean amtk) {
        return new Parameter.Builder(name, description).scaleType(scaleType).prior(PriorType.GAMMA_PRIOR)
                .initial(initial).shape(shape).scale(scale).isNonNegative(true).isPriorFixed(priorFixed).
                        isAdaptiveMultivariateCompatible(amtk)
                .build(parameters);
    }

    public Parameter createCachedGammaPrior(String name, String description, PriorScaleType scaleType, double initial,
                                            double shape, double scale, boolean priorFixed) {
        return new Parameter.Builder(name, description).scaleType(scaleType).prior(PriorType.GAMMA_PRIOR).initial(initial)
                .shape(shape).scale(scale).isNonNegative(true).isPriorFixed(priorFixed).isCached(true).build(parameters);
    }

    public Parameter createParameterOneOverXPrior(String name, String description, PriorScaleType scaleType, double initial) {
        return new Parameter.Builder(name, description).scaleType(scaleType).prior(PriorType.ONE_OVER_X_PRIOR)
                .initial(initial).isNonNegative(true).build(parameters);
    }

    public Parameter createParameterExponentialPrior(String name, String description, PriorScaleType scaleType, double initial,
                                                     double mean, double offset) {
        return createParameterExponentialPrior(name, description, scaleType, initial, mean, offset, false);
    }

    public Parameter createParameterExponentialPrior(String name, String description, PriorScaleType scaleType, double initial,
                                                     double mean, double offset, boolean amtk) {
        return new Parameter.Builder(name, description).scaleType(scaleType).prior(PriorType.EXPONENTIAL_PRIOR)
                .initial(initial).mean(mean).offset(offset).isNonNegative(true).isAdaptiveMultivariateCompatible(amtk).build(parameters);
    }

    public Parameter createParameterLognormalPrior(String name, String description, PriorScaleType scaleType, double initial,
                                                   double mean, double stdev, double offset) {
        return createParameterLognormalPrior(name, description, scaleType, initial, mean, stdev, offset, false);
    }

    public Parameter createParameterLognormalPrior(String name, String description, PriorScaleType scaleType, double initial,
                                                   double mean, double stdev, double offset, boolean amtk) {
        return new Parameter.Builder(name, description).scaleType(scaleType).prior(PriorType.LOGNORMAL_PRIOR)
                .initial(initial).mean(mean).stdev(stdev).offset(offset).isNonNegative(true)
                .isAdaptiveMultivariateCompatible(amtk)
                .build(parameters);
    }

    public Parameter createParameterNormalPrior(String name, String description, PriorScaleType scaleType, double initial,
                                                double mean, double stdev, double offset) {
        return createParameterNormalPrior(name, description, scaleType, initial, mean, stdev, offset, false);
    }

    public Parameter createParameterNormalPrior(String name, String description, PriorScaleType scaleType, double initial,
                                                double mean, double stdev, double offset, boolean amtk) {
        return new Parameter.Builder(name, description).scaleType(scaleType).prior(PriorType.NORMAL_PRIOR)
                .initial(initial).mean(mean).stdev(stdev).offset(offset).isAdaptiveMultivariateCompatible(amtk)
                .build(parameters);
    }

    public Parameter createParameterLaplacePrior(String name, String description, PriorScaleType scaleType, double initial,
                                                 double mean, double scale) {
        return new Parameter.Builder(name, description).scaleType(scaleType).prior(PriorType.LAPLACE_PRIOR)
                .initial(initial).mean(mean).scale(scale).build(parameters);
    }

    public Parameter createParameterBetaDistributionPrior(String name, String description, double initial,
                                                          double shape, double shapeB, double offset) {
        return new Parameter.Builder(name, description).prior(PriorType.BETA_PRIOR).initial(initial)
                .isZeroOne(true).shape(shape).shapeB(shapeB).offset(offset).build(parameters);
    }

    public Parameter createDuplicate(String name, String description, Parameter source) {
        return new Parameter.Builder(name, description).duplicate(source).build(parameters);
    }

    //+++++++++++++++++++ Create Statistic ++++++++++++++++++++++++++++++++
    public Parameter createDiscreteStatistic(String name, String description) { // Poisson Prior
        return new Parameter.Builder(name, description).isDiscrete(true).isStatistic(true)
                .prior(PriorType.POISSON_PRIOR).mean(Math.log(2)).build(parameters);
    }

    public Parameter createStatistic(String name, String description) {
        return new Parameter.Builder(name, description).isStatistic(true).prior(PriorType.NONE_STATISTIC).build(parameters);
    }

    public Parameter createNonNegativeStatistic(String name, String description) {
        return new Parameter.Builder(name, description).isStatistic(true).prior(PriorType.NONE_STATISTIC).isNonNegative(true).build(parameters);
    }
    //+++++++++++++++++++ Create Operator ++++++++++++++++++++++++++++++++
    public Operator createOperator(String parameterName, OperatorType type, double tuning, double weight) {
        Parameter parameter = getParameter(parameterName);
        return new Operator.Builder(parameterName, parameterName, parameter, type, tuning, weight)
                .build(operators);
    }

    public Operator createOperator(String parameterName, OperatorType type, double tuning, double weight, boolean autoOptimize) {
        Parameter parameter = getParameter(parameterName);
        return new Operator.Builder(parameterName, parameterName, parameter, type, tuning, weight)
                .autoOptimize(autoOptimize)
                .build(operators);
    }

    public Operator createScaleOperator(String parameterName, double tuning, double weight) {
        return createScaleOperator(parameterName, null, tuning, weight);
    }

    public Operator createScaleOperator(String parameterName, String description, double tuning, double weight) {
        Parameter parameter = getParameter(parameterName);
        if (description == null) {
            description = (parameter.getDescription() == null ? parameterName : parameter.getDescription());
        }
        return new Operator.Builder(parameterName, description, parameter, OperatorType.SCALE, tuning, weight).build(operators);
    }

    public Operator createOperator(String key, String name, String description, String parameterName, OperatorType type,
                                   double tuning, double weight) {
        Parameter parameter = null;
        if (parameterName != null) {
            parameter = getParameter(parameterName);
        }
        return operators.put(key, new Operator.Builder(name, description, parameter, type, tuning, weight).build()); // key != name
    }

    public Operator createOperatorUsing2Parameters(String key, String name, String description, String parameterName1, String parameterName2,
                                                   OperatorType type, double tuning, double weight) {
        Parameter parameter1 = getParameter(parameterName1);
        Parameter parameter2 = getParameter(parameterName2);
        return operators.put(key, new Operator.Builder(name, description, parameter1, type, tuning, weight).parameter2(parameter2).build());
    }

    public Operator createUpDownOperator(String key, String name, String description, Parameter parameter1, Parameter parameter2,
                                         OperatorType type, double tuning, double weight) {
        return operators.put(key, new Operator.Builder(name, description, parameter1, type, tuning, weight)
                .parameter2(parameter2).build());
    }

    public Operator createBitFlipInSubstitutionModelOperator(String key, String name, String description, Parameter parameter,
                                                             PartitionOptions options, double tuning, double weight) {
//        Parameter parameter = getParameter(parameterName);
        return operators.put(key, new Operator.Builder(name, description, parameter, OperatorType.BITFLIP_IN_SUBST, tuning, weight)
                .partitionOptions(options).build());
    }

    public Operator createDuplicate(String name, String description, Parameter parameter, Operator source) {
        return new Operator.Builder(name, description, parameter, source.getOperatorType(), source.getTuning(), source.getWeight()).build(operators);
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

    public boolean parameterExists(String name) {
        Parameter parameter = parameters.get(name);
        return parameter != null;
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
        if (operator == null) {
            throw new IllegalArgumentException("Operator with name, " + name + ", is unknown");
        }
        return operator;
    }

    protected void addComponent(ComponentOptions component) {
        components.add(component);
        component.createParameters(this);
    }

    protected boolean hasComponent(ComponentFactory factory) {
        for (ComponentOptions component : components) {
            if (factory.getOptionsClass().isAssignableFrom(component.getClass())) {
                return true;
            }
        }
        return false;
    }

    public ComponentOptions getComponentOptions(Class<?> theClass) {
        for (ComponentOptions component : components) {
            if (theClass.isAssignableFrom(component.getClass())) {
                return component;
            }
        }

        throw new IllegalArgumentException("A component options of class, " + theClass.toString() + ", has not been registered");
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

    public Map<String, Parameter> getParameters() {
        return parameters;
    }

    public Map<TaxonList, Parameter> getStatistics() {
        return statistics;
    }

    public Map<String, Operator> getOperators() {
        return operators;
    }

    @Deprecated
    public String noDuplicatedPrefix1(String a , String b) {
        // @todo AR - this doesn't make sense.
        if (a.equals(b)) {
            return a;
        } else {
            return a + b;
        }
    }
}
