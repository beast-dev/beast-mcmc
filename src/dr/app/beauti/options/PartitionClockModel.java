/*
 * PartitionClockModel.java
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

import dr.app.beauti.types.*;
import dr.evolution.datatype.DataType;
import dr.evolution.util.Taxa;

import java.util.ArrayList;
import java.util.List;

import static dr.app.beauti.options.BeautiOptions.DEFAULT_QUANTILE_RELAXED_CLOCK;

/**
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @author Walter Xie
 * @version $Id$
 */
public class PartitionClockModel extends PartitionOptions {
    private static final long serialVersionUID = -6904595851602060488L;

    private ClockType clockType = ClockType.STRICT_CLOCK;
    private ClockDistributionType clockDistributionType = ClockDistributionType.LOGNORMAL;
    private boolean continuousQuantile = DEFAULT_QUANTILE_RELAXED_CLOCK; // a switch at the top of BeautiOptions
    private boolean performModelAveraging = false;

    private PartitionTreeModel treeModel = null;

    private final AbstractPartitionData partition;
    private final int dataLength;

    public PartitionClockModel(final BeautiOptions options, String name, AbstractPartitionData partition, PartitionTreeModel treeModel) {
        super(options, name);

        dataLength = partition.getSiteCount();

        this.partition = partition;
        this.treeModel = treeModel;

        initModelParametersAndOpererators();
    }

    /**
     * A copy constructor
     *
     * @param options the beauti options
     * @param name    the name of the new model
     * @param source  the source model
     */
    public PartitionClockModel(BeautiOptions options, String name, PartitionClockModel source) {
        super(options, name);

        this.clockType = source.clockType;
        clockDistributionType = source.clockDistributionType;

        dataLength = source.dataLength;

        this.partition = source.partition;
        this.treeModel = source.treeModel;

        initModelParametersAndOpererators();
    }

    public void initModelParametersAndOpererators() {
        double rate = 1.0;

        new Parameter.Builder("clock.rate", "substitution rate")
                .prior(PriorType.CTMC_RATE_REFERENCE_PRIOR).initial(rate)
                .isCMTCRate(true).isNonNegative(true).partitionOptions(this)
                .isAdaptiveMultivariateCompatible(true).build(parameters);

        new Parameter.Builder(ClockType.UCED_MEAN, "uncorrelated exponential relaxed clock mean").
                prior(PriorType.CTMC_RATE_REFERENCE_PRIOR).initial(rate)
                .isCMTCRate(true).isNonNegative(true).partitionOptions(this)
                .isAdaptiveMultivariateCompatible(true).build(parameters);

        new Parameter.Builder(ClockType.UCLD_MEAN, "uncorrelated lognormal relaxed clock mean").
                prior(PriorType.CTMC_RATE_REFERENCE_PRIOR).initial(rate)
                .isCMTCRate(true).isNonNegative(true).partitionOptions(this)
                .isAdaptiveMultivariateCompatible(true).build(parameters);

        new Parameter.Builder(ClockType.UCGD_MEAN, "uncorrelated gamma relaxed clock mean").
                prior(PriorType.CTMC_RATE_REFERENCE_PRIOR).initial(rate)
                .isCMTCRate(true).isNonNegative(true).partitionOptions(this)
                .isAdaptiveMultivariateCompatible(true).build(parameters);

        new Parameter.Builder(ClockType.UCLD_STDEV, "uncorrelated lognormal relaxed clock stdev").
                scaleType(PriorScaleType.LOG_STDEV_SCALE).prior(PriorType.EXPONENTIAL_PRIOR).isNonNegative(true)
                .initial(1.0 / 3.0).mean(1.0 / 3.0).offset(0.0).partitionOptions(this)
                .isAdaptiveMultivariateCompatible(true).build(parameters);

        new Parameter.Builder(ClockType.UCGD_SHAPE, "uncorrelated gamma relaxed clock shape").
                prior(PriorType.EXPONENTIAL_PRIOR).isNonNegative(true)
                .initial(1.0 / 3.0).mean(1.0 / 3.0).offset(0.0).partitionOptions(this)
                .isAdaptiveMultivariateCompatible(true).build(parameters);

        // Random local clock
        createParameterGammaPrior(ClockType.LOCAL_CLOCK + ".relativeRates", "random local clock relative rates",
                PriorScaleType.SUBSTITUTION_RATE_SCALE, 1.0, 0.5, 2.0, false, false);
        createParameter(ClockType.LOCAL_CLOCK + ".changes", "random local clock rate change indicator");

        // Autocorrelated clock
        // TODO: implement AC
//        createParameterGammaPrior("branchRates.var", "autocorrelated lognormal relaxed clock rate variance",
//                PriorScaleType.LOG_VAR_SCALE, 0.1, 1, 0.0001, false);
//
//        createParameterGammaPrior("treeModel.rootRate", "autocorrelated lognormal relaxed clock root rate",
//                PriorScaleType.ROOT_RATE_SCALE, 1.0, 1, 0.0001, false);
//        createParameterUniformPrior("treeModel.nodeRates", "autocorrelated lognormal relaxed clock non-root rates",
//                PriorScaleType.SUBSTITUTION_RATE_SCALE, 1.0, 0.0, Parameter.UNIFORM_MAX_BOUND);
//        createParameterUniformPrior("treeModel.allRates", "autocorrelated lognormal relaxed clock all rates",
//                PriorScaleType.SUBSTITUTION_RATE_SCALE, 1.0, 0.0, Parameter.UNIFORM_MAX_BOUND);
//
//        createScaleOperator("branchRates.var", demoTuning, rateWeights);
//
//        createOperator("scaleRootRate", "treeModel.rootRate", "Scales root rate", "treeModel.rootRate", OperatorType.SCALE, demoTuning, rateWeights);
//        createOperator("scaleOneRate", "treeModel.nodeRates", "Scales one non-root rate", "treeModel.nodeRates",
//                OperatorType.SCALE, demoTuning, branchWeights);
//        createOperator("scaleAllRates", "treeModel.allRates", "Scales all rates simultaneously", "treeModel.allRates",
//                OperatorType.SCALE_ALL, demoTuning, rateWeights);
//        createOperator("scaleAllRatesIndependently", "treeModel.nodeRates", "Scales all non-root rates independently", "treeModel.nodeRates",
//                OperatorType.SCALE_INDEPENDENTLY, demoTuning, rateWeights);

        // Uncorrelated clock
        createParameter("branchRates.categories", "relaxed clock branch rate categories");

        //createZeroOneParameter("branchRates.quantiles", "relaxed clock branch rate quantiles", 0.5);
        createZeroOneParameterUniformPrior("branchRates.quantiles", "relaxed clock branch rate quantiles", 0.5);

        // Model averaging
        //createParameter("branchRates.distributionIndex", "distribution integer index");

        new Parameter.Builder("branchRates.distributionIndex", "distribution integer index").
                prior(PriorType.DISCRETE_UNIFORM_PRIOR).isNonNegative(true)
                .initial(0.0).uniformLower(0.0).uniformUpper(2.0).offset(0.0).partitionOptions(this)
                .isAdaptiveMultivariateCompatible(false).build(parameters);

        createScaleOperator("clock.rate", demoTuning, rateWeights);
        createScaleOperator(ClockType.UCED_MEAN, demoTuning, rateWeights);
        createScaleOperator(ClockType.UCLD_MEAN, demoTuning, rateWeights);
        createScaleOperator(ClockType.UCLD_STDEV, demoTuning, rateWeights);
        createScaleOperator(ClockType.UCGD_MEAN, demoTuning, rateWeights);
        createScaleOperator(ClockType.UCGD_SHAPE, demoTuning, rateWeights);

        // Random local clock
        createScaleOperator(ClockType.LOCAL_CLOCK + ".relativeRates", demoTuning, treeWeights);
        createOperator(ClockType.LOCAL_CLOCK + ".changes", OperatorType.BITFLIP, 1, treeWeights);
        createDiscreteStatistic("rateChanges", "number of random local clocks"); // POISSON_PRIOR

        // A vector of relative rates across all partitions...

        createNonNegativeParameterDirichletPrior("allNus", "relative rates amongst partitions parameter", this, 0, 1.0, true, true);
        createOperator("deltaNus", "allNus",
                "Change partition rates relative to each other maintaining mean", "allNus",
                OperatorType.DELTA_EXCHANGE, 0.01, 3.0);

        createNonNegativeParameterInfinitePrior("allMus", "relative rates amongst partitions parameter", this, PriorScaleType.SUBSTITUTION_PARAMETER_SCALE, 1.0, true);
        createOperator("deltaMus", "allMus",
                "Change partition rates relative to each other maintaining mean", "allMus",
                OperatorType.WEIGHTED_DELTA_EXCHANGE, 0.01, 3.0);


        createOperator("swapBranchRateCategories", "branchRates.categories", "Performs a swap of branch rate categories",
                "branchRates.categories", OperatorType.SWAP, 1, branchWeights / 3);
        createOperator("uniformBranchRateCategories", "branchRates.categories", "Performs an integer uniform draw of branch rate categories",
                "branchRates.categories", OperatorType.INTEGER_UNIFORM, 1, branchWeights / 3);

        if (!options.classicOperatorsAndPriors) {
            createOperator("rwBranchRateQuantiles", "branchRates.quantiles", "Random walk of branch rate quantiles",
                    "branchRates.quantiles", OperatorType.RANDOM_WALK_LOGIT, 1, branchWeights / 3);
        } else {
            createOperator("uniformBranchRateQuantiles", "branchRates.quantiles", "Performs an uniform draw of branch rate quantiles",
                    "branchRates.quantiles", OperatorType.UNIFORM, 0, branchWeights / 3);
        }

        createOperator("uniformBranchRateDistributionIndex", "branchRates.distributionIndex", "Performs a uniform draw of the distribution index",
                "branchRates.distributionIndex", OperatorType.INTEGER_UNIFORM, 0, branchWeights / 3);

        createUpDownOperator("upDownRateHeights", "Substitution rate and heights",
                "Scales substitution rates inversely to node heights of the tree",
                getPartitionTreeModel().getParameter("treeModel.allInternalNodeHeights"),
                getParameter("clock.rate"), OperatorType.UP_DOWN, demoTuning, rateWeights);
        createUpDownOperator("upDownUCEDMeanHeights", "UCED mean and heights",
                "Scales UCED mean inversely to node heights of the tree",
                getPartitionTreeModel().getParameter("treeModel.allInternalNodeHeights"),
                getParameter(ClockType.UCED_MEAN), OperatorType.UP_DOWN, demoTuning, rateWeights);
        createUpDownOperator("upDownUCLDMeanHeights", "UCLD mean and heights",
                "Scales UCLD mean inversely to node heights of the tree",
                getPartitionTreeModel().getParameter("treeModel.allInternalNodeHeights"),
                getParameter(ClockType.UCLD_MEAN), OperatorType.UP_DOWN,
                demoTuning, rateWeights);
        createUpDownOperator("upDownUCGDMeanHeights", "UCGD mean and heights",
                "Scales UCGD mean inversely to node heights of the tree",
                getPartitionTreeModel().getParameter("treeModel.allInternalNodeHeights"),
                getParameter(ClockType.UCGD_MEAN), OperatorType.UP_DOWN, demoTuning, rateWeights);

        createUpDownOperator("microsatUpDownRateHeights", "Substitution rate and heights",
                "Scales substitution rates inversely to node heights of the tree",
                getPartitionTreeModel().getParameter("treeModel.allInternalNodeHeights"),
                getParameter("clock.rate"),OperatorType.MICROSAT_UP_DOWN, demoTuning, branchWeights);
    }

    // From PartitionClockModelTreeModelLink
//    public List<Parameter> selectParameters(List<Parameter> params) {
//        setAvgRootAndRate();
//        getParameter("branchRates.categories");
//        getParameter("treeModel.rootRate");
//        getParameter("treeModel.nodeRates");
//        getParameter("treeModel.allRates");
//
//        if (options.hasData()) {
//            // if not fixed then do mutation rate move and up/down move
//            boolean fixed = !model.isEstimatedRate();
//
//            Parameter rateParam;
//
//            switch (model.getClockType()) {
//                case AUTOCORRELATED:
//                    rateParam = getParameter("treeModel.rootRate");
//                    rateParam.isFixed = fixed;
//                    if (!fixed) params.add(rateParam);
//
//                    params.add(getParameter("branchRates.var"));
//                    break;
//            }
//        }
//        return params;
//    }

    @Override
    public List<Parameter> selectParameters(List<Parameter> params) {
//        setAvgRootAndRate();
        double rate = 1.0;

        if (options.hasData()) {

            switch (clockType) {
                case STRICT_CLOCK:
                    params.add(getClockRateParameter());
                    break;

                case RANDOM_LOCAL_CLOCK:
                    params.add(getClockRateParameter());
                    getParameter(ClockType.LOCAL_CLOCK + ".changes");
                    params.add(getParameter("rateChanges"));
                    params.add(getParameter(ClockType.LOCAL_CLOCK + ".relativeRates"));
                    break;

                case FIXED_LOCAL_CLOCK:
                    params.add(getClockRateParameter());
                    for (Taxa taxonSet : options.taxonSets) {
                        if (options.taxonSetsMono.get(taxonSet)) {
                            String parameterName = taxonSet.getId() + ".rate";
                            if (!hasParameter(parameterName)) {
                                new Parameter.Builder(parameterName, "substitution rate")
                                        .prior(PriorType.CTMC_RATE_REFERENCE_PRIOR)
                                        .initial(rate)
                                        .isCMTCRate(true).isNonNegative(true)
                                        .partitionOptions(this)
                                        .taxonSet(taxonSet)
                                        .build(parameters);
                                createScaleOperator(parameterName, demoTuning, rateWeights);
                            }

                            params.add(getParameter(taxonSet.getId() + ".rate"));
                        }
                    }
                    break;

                case UNCORRELATED:
                    // add the scale parameter (if needed) for the distribution. The location parameter will be added
                    // in getClockRateParameter.
                    switch (clockDistributionType) {
                        case LOGNORMAL:
                            params.add(getClockRateParameter());
                            params.add(getParameter(ClockType.UCLD_STDEV));
                            break;
                        case GAMMA:
                            params.add(getClockRateParameter());
                            params.add(getParameter(ClockType.UCGD_SHAPE));
                            break;
                        case CAUCHY:
                            throw new UnsupportedOperationException("Uncorrelated Cauchy clock not implemented yet");
//                            break;
                        case EXPONENTIAL:
                            params.add(getClockRateParameter());
                            break;
                        case MODEL_AVERAGING:
                            params.add(getClockRateParameter(ClockType.UNCORRELATED, ClockDistributionType.LOGNORMAL));
                            params.add(getParameter(ClockType.UCLD_STDEV));
                            params.add(getClockRateParameter(ClockType.UNCORRELATED, ClockDistributionType.GAMMA));
                            params.add(getParameter(ClockType.UCGD_SHAPE));
                            params.add(getClockRateParameter(ClockType.UNCORRELATED, ClockDistributionType.EXPONENTIAL));
                            params.add(getParameter("branchRates.quantiles"));
                            params.add(getParameter("branchRates.distributionIndex"));
                            break;
                    }
                    break;

                case AUTOCORRELATED:
                    throw new UnsupportedOperationException("Autocorrelated clock not implemented yet");
//                    params.add(getParameter("branchRates.var"));
//                    break;

                default:
                    throw new IllegalArgumentException("Unknown clock model");
            }
            //Parameter rateParam = getClockRateParameter();
            //params.add(rateParam);
        }

        return params;
    }

    public Parameter getClockRateParameter() {
        return getClockRateParameter(clockType, clockDistributionType);
    }

    private Parameter getClockRateParameter(ClockType clockType, ClockDistributionType clockDistributionType) {
        Parameter rateParam = null;
        switch (clockType) {
            case STRICT_CLOCK:
                rateParam = getParameter("clock.rate");
                break;

            case RANDOM_LOCAL_CLOCK:
                rateParam = getParameter("clock.rate");
                break;

            case FIXED_LOCAL_CLOCK:
                rateParam = getParameter("clock.rate");
                break;

            case UNCORRELATED:
                switch (clockDistributionType) {
                    case LOGNORMAL:
                        rateParam = getParameter(ClockType.UCLD_MEAN);
                        break;
                    case GAMMA:
                        rateParam = getParameter(ClockType.UCGD_MEAN);
                        break;
                    case CAUCHY:
                        throw new UnsupportedOperationException("Uncorrelated Cauchy clock not implemented yet");
//                            break;
                    case EXPONENTIAL:
                        rateParam = getParameter(ClockType.UCED_MEAN);
                        break;
                    case MODEL_AVERAGING:
                        break;
                }
                break;

            case AUTOCORRELATED:
                throw new UnsupportedOperationException("Autocorrelated clock not implemented yet");
//                    rateParam = getParameter("treeModel.rootRate");//TODO fix tree?
//                    break;

            default:
                throw new IllegalArgumentException("Unknown clock model");
        }

        if (!rateParam.isPriorEdited()) {
            if (options.treeModelOptions.isNodeCalibrated(partition.treeModel) < 0
                    && !options.clockModelOptions.isTipCalibrated()) {
                rateParam.setFixed(true);
            } else {
                rateParam.priorType = PriorType.CTMC_RATE_REFERENCE_PRIOR;
            }
        }

        return rateParam;
    }

    private Operator getUpDownOperator() {
        switch (clockType) {
            case STRICT_CLOCK:
                return getOperator("upDownRateHeights");
            case RANDOM_LOCAL_CLOCK:
                return getOperator("upDownRateHeights");
            case FIXED_LOCAL_CLOCK:
                return getOperator("upDownRateHeights");

            case UNCORRELATED:
                switch (clockDistributionType) {
                    case LOGNORMAL:
                        return getOperator("upDownUCLDMeanHeights");
                    case GAMMA:
                        return getOperator("upDownUCGDMeanHeights");
                    case CAUCHY:
                        throw new UnsupportedOperationException("Uncorrelated Cauchy clock not implemented yet");
//                            break;
                    case EXPONENTIAL:
                        return getOperator("upDownUCEDMeanHeights");
                }
                break;

            case AUTOCORRELATED:
                throw new UnsupportedOperationException("Autocorrelated clock not implemented yet");
//                    rateParam = getParameter("treeModel.rootRate");//TODO fix tree?
//                    break;

            default:
                throw new IllegalArgumentException("Unknown clock model");
        }
        return null;
    }

    @Override
    public List<Operator> selectOperators(List<Operator> operators) {
        List<Operator> ops = new ArrayList<Operator>();

        if (options.hasData()) {
            Operator op;
            if (getDataType().getType() == DataType.MICRO_SAT) {
                if (getClockType() == ClockType.STRICT_CLOCK) {
                    op = getOperator("microsatUpDownRateHeights");
                    ops.add(op);
                } else {
                    throw new UnsupportedOperationException("Microsatellite only supports strict clock model");
                }

            } else {

                Operator rateOperator = getOperator("clock.rate");
                switch (clockType) {
                    case STRICT_CLOCK:
                        ops.add(rateOperator);
                        break;

                    case RANDOM_LOCAL_CLOCK:
                        ops.add(rateOperator);
                        ops.add(getOperator(ClockType.LOCAL_CLOCK + ".relativeRates"));
                        ops.add(getOperator(ClockType.LOCAL_CLOCK + ".changes"));
                        break;

                    case FIXED_LOCAL_CLOCK:
                        ops.add(rateOperator);
                        for (Taxa taxonSet : options.taxonSets) {
                            if (options.taxonSetsMono.get(taxonSet)) {
                                ops.add(getOperator(taxonSet.getId() + ".rate"));
                            }
                        }
                        break;

                    case UNCORRELATED:
                        switch (clockDistributionType) {
                            case LOGNORMAL:
                                ops.add(rateOperator = getOperator(ClockType.UCLD_MEAN));
                                ops.add(getOperator(ClockType.UCLD_STDEV));
                                isOperatorParameterFixed(ops, rateOperator);
                                break;
                            case GAMMA:
                                ops.add(rateOperator = getOperator(ClockType.UCGD_MEAN));
                                ops.add(getOperator(ClockType.UCGD_SHAPE));
                                isOperatorParameterFixed(ops, rateOperator);
                                break;
                            case CAUCHY:
//                                throw new UnsupportedOperationException("Uncorrelated Couchy clock not implemented yet");
                                break;
                            case EXPONENTIAL:
                                ops.add(rateOperator = getOperator(ClockType.UCED_MEAN));
                                isOperatorParameterFixed(ops, rateOperator);
                                break;
                            case MODEL_AVERAGING:
                                ops.add(getOperator(ClockType.UCLD_MEAN));
                                ops.add(getOperator(ClockType.UCLD_STDEV));
                                ops.add(getOperator(ClockType.UCGD_MEAN));
                                ops.add(getOperator(ClockType.UCGD_SHAPE));
                                ops.add(getOperator(ClockType.UCED_MEAN));

                                if (!getOperator(ClockType.UCLD_MEAN).isParameterFixed()) {
                                    ops.add(getOperator("upDownUCLDMeanHeights"));
                                }
                                if (!getOperator(ClockType.UCGD_MEAN).isParameterFixed()) {
                                    ops.add(getOperator("upDownUCGDMeanHeights"));
                                }
                                if (!getOperator(ClockType.UCED_MEAN).isParameterFixed()) {
                                    ops.add(getOperator("upDownUCEDMeanHeights"));
                                }

                                ops.add(getOperator("uniformBranchRateDistributionIndex"));
                                break;
                        }

                        if (isContinuousQuantile()) {
                            if (!options.classicOperatorsAndPriors) {
                                ops.add(getOperator("rwBranchRateQuantiles"));
                            } else {
                                ops.add(getOperator("uniformBranchRateQuantiles"));
                            }
                        } else {
                            ops.add(getOperator("swapBranchRateCategories"));
                            ops.add(getOperator("uniformBranchRateCategories"));
                        }
                        break;

                    case AUTOCORRELATED:
                        throw new UnsupportedOperationException("Autocorrelated clock not implemented yet");
//                        break;

                    default:
                        throw new IllegalArgumentException("Unknown clock model");
                }

            }
        }

        Parameter allMusNus = getParameter(!options.classicOperatorsAndPriors && options.NEW_RELATIVE_RATE_PARAMETERIZATION ? "allNus" : "allMus");

        if (allMusNus.getSubParameters().size() > 1) {
            ops.add(getOperator(!options.classicOperatorsAndPriors && options.NEW_RELATIVE_RATE_PARAMETERIZATION ? "deltaNus" : "deltaMus"));
        }

        if (options.operatorSetType != OperatorSetType.CUSTOM) {
            // unless a custom mix has been chosen these operators should be off if AMTK is on
            for (Operator op : ops) {
                if (op.getParameter1().isAdaptiveMultivariateCompatible) {
                    op.setUsed(options.operatorSetType != OperatorSetType.ADAPTIVE_MULTIVARIATE);
                }
            }
        }

        operators.addAll(ops);

        return ops;
    }

    private void isOperatorParameterFixed(List<Operator> ops, Operator rateOperator) {
        if (!rateOperator.isParameterFixed()) {
            Operator upDownOperator = getUpDownOperator();
            // need to set the node heights parameter again in case the treeModel has changed and
            upDownOperator.setParameter1(
                    getPartitionTreeModel().getParameter("treeModel.allInternalNodeHeights"));
            ops.add(upDownOperator);
        } else {
            ops.add(getPartitionTreeModel().getOperator("treeModel.allInternalNodeHeights"));
        }
    }

    /////////////////////////////////////////////////////////////
    public void setClockType(ClockType clockType) {
        this.clockType = clockType;
    }

    public ClockType getClockType() {
        return clockType;
    }

    public ClockDistributionType getClockDistributionType() {
        return clockDistributionType;
    }

    public void setClockDistributionType(final ClockDistributionType clockDistributionType) {
        this.clockDistributionType = clockDistributionType;
    }

    public boolean performModelAveraging() {
        return performModelAveraging;
    }

    public void setPerformModelAveraging(boolean performModelAveraging) {
        this.performModelAveraging = performModelAveraging;
    }

    public boolean isContinuousQuantile() {
        return continuousQuantile;
    }

    public void setContinuousQuantile(boolean continuousQuantile) {
        this.continuousQuantile = continuousQuantile;
    }

    public String getPrefix() {
        String prefix = "";
        if (options.getPartitionClockModels().size() > 1) {
            // There is more than one active partition model
            prefix += getName() + ".";
        }
        return prefix;
    }

    public PartitionTreeModel getPartitionTreeModel() {
        return treeModel;
    }

    public void setPartitionTreeModel(PartitionTreeModel treeModel) {
        options.clearDataPartitionCaches();
        this.treeModel = treeModel;
    }

    public void copyFrom(PartitionClockModel source) {
        clockType = source.clockType;
        clockDistributionType = source.clockDistributionType;
        continuousQuantile = source.continuousQuantile;
        performModelAveraging = source.performModelAveraging;
    }
}
