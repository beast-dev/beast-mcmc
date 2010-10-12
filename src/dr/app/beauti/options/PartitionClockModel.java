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

import dr.app.beauti.enumTypes.*;

import java.util.List;

/**
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @author Walter Xie
 * @version $Id$
 */
public class PartitionClockModel extends PartitionOptions {

    // Instance variables
    private final BeautiOptions options;

    private ClockType clockType = ClockType.STRICT_CLOCK;
    private boolean isEstimatedRate = true;
    private double rate = 1.0;

    public PartitionClockModel(BeautiOptions options, PartitionData partition) {
        this.options = options;
        this.partitionName = partition.getName();

        initClockModelParaAndOpers();
    }

    /**
     * A copy constructor
     *
     * @param options the beauti options
     * @param name    the name of the new model
     * @param source  the source model
     */
    public PartitionClockModel(BeautiOptions options, String name, PartitionClockModel source) {
        this.options = options;
        this.partitionName = name;
        this.clockType = source.clockType;

        initClockModelParaAndOpers();
    }

//    public PartitionClockModel(BeautiOptions options, String name) {
//        this.options = options;
//        this.name = name;
//    }

    private void initClockModelParaAndOpers() {

        int dataLength = 0;
        for (PartitionData partitionData : options.getAllPartitionData(this)) {
            dataLength += partitionData.getSiteCount();
        }
        if (dataLength <= 1) { // TODO Discuss threshold

            createParameterClockRateUndefinedPrior(this, "clock.rate", "substitution rate",
                    PriorScaleType.SUBSTITUTION_RATE_SCALE, rate, 0.0, Double.POSITIVE_INFINITY);
            createParameterClockRateUndefinedPrior(this, ClockType.UCED_MEAN, "uncorrelated exponential relaxed clock mean",
                    PriorScaleType.SUBSTITUTION_RATE_SCALE, rate, 0.0, Double.POSITIVE_INFINITY);
            createParameterClockRateUndefinedPrior(this, ClockType.UCLD_MEAN, "uncorrelated lognormal relaxed clock mean",
                    PriorScaleType.SUBSTITUTION_RATE_SCALE, rate, 0.0, Double.POSITIVE_INFINITY);
        } else {

            createParameterClockRateUniform(this, "clock.rate", "substitution rate. ",
                    PriorScaleType.SUBSTITUTION_RATE_SCALE, rate, 0.0, Double.POSITIVE_INFINITY);
            createParameterClockRateUniform(this, ClockType.UCED_MEAN, "uncorrelated exponential relaxed clock mean. ",
                    PriorScaleType.SUBSTITUTION_RATE_SCALE, rate, 0.0, Double.POSITIVE_INFINITY);
            createParameterClockRateUniform(this, ClockType.UCLD_MEAN, "uncorrelated lognormal relaxed clock mean. ",
                    PriorScaleType.SUBSTITUTION_RATE_SCALE, rate, 0.0, Double.POSITIVE_INFINITY);
        }
        createParameterClockRateExponential(this, ClockType.UCLD_STDEV, "uncorrelated lognormal relaxed clock stdev",
                PriorScaleType.LOG_STDEV_SCALE, 1.0 / 3.0, 1.0 / 3.0, 0.0, 0.0, Double.POSITIVE_INFINITY);
        // Random local clock
        createParameterGammaPrior(ClockType.LOCAL_CLOCK + ".relativeRates", "random local clock relative rates",
                PriorScaleType.SUBSTITUTION_RATE_SCALE, 1.0, 0.5, 2.0, 0.0, Double.POSITIVE_INFINITY, false);
        createParameter(ClockType.LOCAL_CLOCK + ".changes", "random local clock rate change indicator");

        createScaleOperator("clock.rate", demoTuning, rateWeights);
        createScaleOperator(ClockType.UCED_MEAN, demoTuning, rateWeights);
        createScaleOperator(ClockType.UCLD_MEAN, demoTuning, rateWeights);
        createScaleOperator(ClockType.UCLD_STDEV, demoTuning, rateWeights);
        // Random local clock
        createScaleOperator(ClockType.LOCAL_CLOCK + ".relativeRates", demoTuning, treeWeights);
        createOperator(ClockType.LOCAL_CLOCK + ".changes", OperatorType.BITFLIP, 1, treeWeights);
        createDiscreteStatistic("rateChanges", "number of random local clocks"); // POISSON_PRIOR
    }

    /**
     * return a list of parameters that are required
     *
     * @param params the parameter list
     */
    public void selectParameters(List<Parameter> params) {
        if (options.hasData()) {
            // if not fixed then do mutation rate move and up/down move
            boolean fixed = !isEstimatedRate;
//
//            List<PartitionClockModel> models = new ArrayList<PartitionClockModel>();
//            models.add(this);
//            double selectedRate = options.clockModelOptions.getSelectedRate(models);

            Parameter rateParam = null;

            switch (clockType) {
                case STRICT_CLOCK:
                    rateParam = getParameter("clock.rate");
                    break;

                case RANDOM_LOCAL_CLOCK:
                    rateParam = getParameter("clock.rate");                     
                    getParameter(ClockType.LOCAL_CLOCK + ".changes");
                    params.add(getParameter("rateChanges"));
                    params.add(getParameter(ClockType.LOCAL_CLOCK + ".relativeRates"));
                    break;

                case UNCORRELATED_EXPONENTIAL:
                    rateParam = getParameter(ClockType.UCED_MEAN);
                    break;

                case UNCORRELATED_LOGNORMAL:
                    rateParam = getParameter(ClockType.UCLD_MEAN);
                    params.add(getParameter(ClockType.UCLD_STDEV));
                    break;

                case AUTOCORRELATED_LOGNORMAL:
//                    rateParam = getParameter("treeModel.rootRate");//TODO fix tree?
//                    params.add(getParameter("branchRates.var"));
                    break;

                default:
                    throw new IllegalArgumentException("Unknown clock model");
            }

//            if (this.getAllPartitionData().get(0) instanceof TraitData) {
//                rateParam.priorType = PriorType.ONE_OVER_X_PRIOR; // 1/location.clock.rate
//            }

            rateParam.isFixed = fixed;
            if (fixed && options.clockModelOptions.getRateOptionClockModel() == FixRateType.RELATIVE_TO
                    && rateParam.initial != rate) {
                rateParam.initial = rate;
                rateParam.setPriorEdited(true);
            }
//            if (options.clockModelOptions.getRateOptionClockModel() == FixRateType.FIX_MEAN
//                     || options.clockModelOptions.getRateOptionClockModel() == FixRateType.RELATIVE_TO) {
//
//                rateParam.priorEdited = true; // important
//            }
//
//            if (!rateParam.priorEdited) {
//                rateParam.initial = selectedRate;
//            }

            if (!fixed) params.add(rateParam);
        }
    }

    /**
     * return a list of operators that are required
     *
     * @param ops the operator list
     */
    public void selectOperators(List<Operator> ops) {
        if (options.hasData()) {

            if ((!(options.clockModelOptions.getRateOptionClockModel() == FixRateType.FIX_MEAN))
                    && isEstimatedRate) {
                switch (clockType) {
                    case STRICT_CLOCK:
                        ops.add(getOperator("clock.rate"));
                        break;

                    case RANDOM_LOCAL_CLOCK:
                        ops.add(getOperator("clock.rate"));
                        addRandomLocalClockOperators(ops);
                        break;

                    case UNCORRELATED_EXPONENTIAL:
                        ops.add(getOperator(ClockType.UCED_MEAN));
                        break;

                    case UNCORRELATED_LOGNORMAL:
                        ops.add(getOperator(ClockType.UCLD_MEAN));
                        ops.add(getOperator(ClockType.UCLD_STDEV));
                        break;

                    case AUTOCORRELATED_LOGNORMAL:
                        //TODO
                        break;

                    default:
                        throw new IllegalArgumentException("Unknown clock model");
                }
            } else {
                switch (clockType) {
                    case STRICT_CLOCK:
                    case UNCORRELATED_EXPONENTIAL:
                    case AUTOCORRELATED_LOGNORMAL:
                        // no parameter to operator on
                        break;

                    case RANDOM_LOCAL_CLOCK:
                        addRandomLocalClockOperators(ops);
                        break;

                    case UNCORRELATED_LOGNORMAL:
                        ops.add(getOperator(ClockType.UCLD_STDEV));
                        break;

                    default:
                        throw new IllegalArgumentException("Unknown clock model");
                }
            }
        }
    }

    private void addRandomLocalClockOperators(List<Operator> ops) {
        ops.add(getOperator(ClockType.LOCAL_CLOCK + ".relativeRates"));
        ops.add(getOperator(ClockType.LOCAL_CLOCK + ".changes"));
    }

    // +++++++++++++++++++++++++++ *BEAST ++++++++++++++++++++++++++++++++++++

    public void iniClockRateStarBEAST() {
        Parameter rateParam = null;
        switch (clockType) {
            case STRICT_CLOCK:
                rateParam = getParameter("clock.rate");
                break;

            case RANDOM_LOCAL_CLOCK:
                rateParam = getParameter("clock.rate");
                getParameter(ClockType.LOCAL_CLOCK + ".relativeRates");
                getParameter(ClockType.LOCAL_CLOCK + ".changes");
                break;

            case UNCORRELATED_EXPONENTIAL:
                rateParam = getParameter(ClockType.UCED_MEAN);
                break;

            case UNCORRELATED_LOGNORMAL:
                rateParam = getParameter(ClockType.UCLD_MEAN);
                break;

            case AUTOCORRELATED_LOGNORMAL:
//                rateParam = getParameter("treeModel.rootRate");//TODO fix tree?
                break;

            default:
                throw new IllegalArgumentException("Unknown clock model");
        }

        double selectedRate = options.clockModelOptions.getSelectedRate(options.getNonTraitsDataList());

        rateParam.priorType = PriorType.GAMMA_PRIOR;
        rateParam.initial = selectedRate;
        rateParam.shape = 0.1;
        rateParam.scale = 10 * selectedRate;
        rateParam.offset = 0;
    }

    /////////////////////////////////////////////////////////////
    public void setClockType(ClockType clockType) {
        this.clockType = clockType;
    }

    public ClockType getClockType() {
        return clockType;
    }

    public void setEstimatedRate(boolean isEstimatedRate) {
        this.isEstimatedRate = isEstimatedRate;
    }

    public boolean isEstimatedRate() {
        return isEstimatedRate;
    }

    public double getRate() {
        return rate;
    }

    public void setRate(double rate) {
        this.rate = rate;
    }

    public String getPrefix() {
        String prefix = "";
        if (options.getPartitionClockModels().size() > 1) { //|| options.isSpeciesAnalysis()
            // There is more than one active partition model
            prefix += getName() + ".";
        }
        return prefix;
    }

}
