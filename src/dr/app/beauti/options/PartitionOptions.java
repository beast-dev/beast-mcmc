/*
 * PartitionOptions.java
 *
 * Copyright (c) 2002-2012 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

import dr.app.beauti.types.PriorScaleType;
import dr.evolution.datatype.DataType;
import dr.math.MathUtils;

import java.util.List;


/**
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @author Walter Xie
 */
public abstract class PartitionOptions extends ModelOptions {

    protected String partitionName;
    protected final BeautiOptions options;

    protected double[] avgRootAndRate = new double[]{1.0, 1.0};

    public PartitionOptions(BeautiOptions options) {
        this.options = options;
    }

    public PartitionOptions(BeautiOptions options, String name) {
        this.options = options;
        this.partitionName = name;

        initModelParametersAndOpererators();
    }

    protected abstract void initModelParametersAndOpererators();

    protected abstract void selectParameters(List<Parameter> params);

    protected abstract void selectOperators(List<Operator> ops);

    public abstract String getPrefix();

    //    protected void createParameterClockRateUndefinedPrior(PartitionOptions options, String name, String description, PriorScaleType scaleType,
//                                                          double initial, double truncationLower, double truncationUpper) { // it will change to Uniform
//        new Parameter.Builder(name, description).scaleType(scaleType).prior(PriorType.UNDEFINED).initial(initial)
//                .isCMTCRate(true).isNonNegative(true)
//                .truncationLower(truncationLower).truncationUpper(truncationUpper).partitionOptions(options).build(parameters);
//    }
//
//    protected void createParameterClockRateReferencePrior(PartitionOptions options, String name, String description, PriorScaleType scaleType,
//                                                          double initial) { // it will change to Uniform
//        new Parameter.Builder(name, description).scaleType(scaleType).prior(PriorType.CMTC_RATE_REFERENCE_PRIOR).initial(initial)
//                .isCMTCRate(true).isNonNegative(true)
//                .truncationLower(truncationLower).truncationUpper(truncationUpper).partitionOptions(options).build(parameters);
//    }
//
//    protected void createParameterClockRateUniform(PartitionOptions options, String name, String description, PriorScaleType scaleType,
//                                                   double initial, double truncationLower, double truncationUpper) {
//        new Parameter.Builder(name, description).scaleType(scaleType).prior(PriorType.UNIFORM_PRIOR).initial(initial)
//                .isCMTCRate(true).isNonNegative(true)
//                .truncationLower(truncationLower).truncationUpper(truncationUpper).partitionOptions(options).build(parameters);
//    }
//
//    protected void createParameterClockRateGamma(PartitionOptions options, String name, String description, PriorScaleType scaleType,
//                                                 double initial, double shape, double scale) {
//        new Parameter.Builder(name, description).scaleType(scaleType).prior(PriorType.GAMMA_PRIOR).initial(initial)
//                .isCMTCRate(true).isNonNegative(true)
//                .shape(shape).scale(scale).partitionOptions(options).build(parameters);
//    }
//
//    public void createParameterClockRateExponential(PartitionOptions options, String name, String description, PriorScaleType scaleType,
//                                                    double initial, double mean, double offset) {
//        new Parameter.Builder(name, description).scaleType(scaleType).prior(PriorType.EXPONENTIAL_PRIOR)
//                .isCMTCRate(true).isNonNegative(true)
//                .initial(initial).mean(mean).offset(offset).partitionOptions(options).build(parameters);
//    }
//
//
    protected void createParameterTree(PartitionOptions options, String name, String description, boolean isNodeHeight, double value) {
        new Parameter.Builder(name, description).isNodeHeight(isNodeHeight).scaleType(PriorScaleType.TIME_SCALE)
                .isNonNegative(true).initial(value).partitionOptions(options).build(parameters);
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

        autoScale(parameter); // not include clock rate, and treeModel.rootHeight

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

    public DataType getDataType() {
        return options.getDataPartitions(this).get(0).getDataType();
    }

    public double[] getAvgRootAndRate() {
        return avgRootAndRate;
    }

    public void setAvgRootAndRate() {
        this.avgRootAndRate = options.clockModelOptions.calculateInitialRootHeightAndRate(options.getDataPartitions(this));
    }

    protected void autoScale(Parameter param) {
        double avgInitialRootHeight = avgRootAndRate[0];
        double avgInitialRate = avgRootAndRate[1];

//        double growthRateMaximum = 1E6;
        double birthRateMaximum = 1E6;
//        double substitutionRateMaximum = 100;
//        double logStdevMaximum = 10;
//        double substitutionParameterMaximum = 100;

//        if (options.clockModelOptions.getRateOptionClockModel() == FixRateType.FIX_MEAN
//                || options.clockModelOptions.getRateOptionClockModel() == FixRateType.RELATIVE_TO) {
//
//            growthRateMaximum = 1E6 * avgInitialRate;
        birthRateMaximum = 1E6 * avgInitialRate;
//        }

//        if (options.clockModelOptions.getRateOptionClockModel() == FixRateType.FIX_MEAN) {
//            double rate = options.clockModelOptions.getMeanRelativeRate();
//
//            growthRateMaximum = 1E6 * rate;
//            birthRateMaximum = 1E6 * rate;
//
//            if (options.hasData()) {
//                initialRootHeight = meanDistance / rate;
//
//                initialRootHeight = round(initialRootHeight, 2);
//            }
//
//        } else {
//            if (options.maximumTipHeight > 0) {
//                initialRootHeight = options.maximumTipHeight * 10.0;
//            }
//
//            initialRate = round((meanDistance * 0.2) / initialRootHeight, 2);
//        }

//        double timeScaleMaximum = MathUtils.round(avgInitialRootHeight * 1000.0, 2);


//        if (!options.hasData()) param.setPriorEdited(false);

        if (!param.isPriorEdited()) {
            switch (param.scaleType) {
                case TIME_SCALE:
//                        param.lower = Math.max(0.0, param.lower);
//                        param.upper = Math.min(timeScaleMaximum, param.upper);
//                    if (param.isNodeHeight) { //TODO only affecting "treeModel.rootHeight", need to review
//                        param.lower = options.maximumTipHeight;
////                    param.upper = timeScaleMaximum;
////                    param.initial = avgInitialRootHeight;
//                            if (param.getOptions() instanceof PartitionTreeModel) { // move to PartitionTreeModel
//                                param.initial = ((PartitionTreeModel) param.getOptions()).getInitialRootHeight();
//                            }
//                    } else {
                    param.initial = avgInitialRootHeight;
//                    }

                    break;
                case LOG_TIME_SCALE:
                    param.initial = Math.log(avgInitialRootHeight);
                    break;

                case T50_SCALE:
//                        param.lower = Math.max(0.0, param.lower);
                    //param.upper = Math.min(timeScaleMaximum, param.upper);
                    param.initial = avgInitialRootHeight / 5.0;
                    break;

                case GROWTH_RATE_SCALE:
                    param.initial = avgInitialRootHeight / 1000;
                    // use Laplace
                    if (param.getBaseName().startsWith("logistic")) {
                        param.scale = Math.log(1000) / avgInitialRootHeight;
//                            System.out.println("logistic");
                    } else {
                        param.scale = Math.log(10000) / avgInitialRootHeight;
//                            System.out.println("not logistic");
                    }
                    break;

                case BIRTH_RATE_SCALE:
//                    param.uniformLower = Math.max(0.0, param.lower);
//                    param.uniformUpper = Math.min(birthRateMaximum, param.upper);
                    param.initial = MathUtils.round(1 / options.treeModelOptions.getExpectedAvgBranchLength(avgInitialRootHeight), 2);
                    break;
                case ORIGIN_SCALE:
                    param.initial = MathUtils.round(avgInitialRootHeight * 1.1, 2);
                    break;

                case SUBSTITUTION_RATE_SCALE:
//                        param.lower = Math.max(0.0, param.lower);
                    //param.upper = Math.min(substitutionRateMaximum, param.upper);
                    param.initial = avgInitialRate;
                    break;

                case LOG_STDEV_SCALE:
//                        param.lower = Math.max(0.0, param.lower);
                    //param.upper = Math.min(logStdevMaximum, param.upper);
                    break;

                case SUBSTITUTION_PARAMETER_SCALE:
//                        param.lower = Math.max(0.0, param.lower);
                    //param.upper = Math.min(substitutionParameterMaximum, param.upper);
                    break;

                // Now have a field 'isZeroOne'
//                case UNITY_SCALE:
//                    param.lower = 0.0;
//                    param.upper = 1.0;
//                    break;

                case ROOT_RATE_SCALE:
                    param.initial = avgInitialRate;
                    param.shape = 0.5;
                    param.scale = param.initial / 0.5;
                    break;

                case LOG_VAR_SCALE:
                    param.initial = avgInitialRate;
                    param.shape = 2.0;
                    param.scale = param.initial / 2.0;
                    break;

            }

        }
    }

    public BeautiOptions getOptions() {
        return options;
    }
}
