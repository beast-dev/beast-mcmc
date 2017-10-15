/*
 * PartitionTreePrior.java
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
import dr.evomodel.coalescent.VariableDemographicModel;
import dr.evomodel.speciation.CalibrationPoints;
import dr.evomodelxml.speciation.BirthDeathEpidemiologyModelParser;
import dr.evomodelxml.speciation.BirthDeathModelParser;
import dr.evomodelxml.speciation.BirthDeathSerialSamplingModelParser;
import dr.math.MathUtils;

import java.util.List;

/**
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @author Walter Xie
 */
public class PartitionTreePrior extends PartitionOptions {

    private static final long serialVersionUID = 8222609132259262723L;

    private TreePriorType nodeHeightPrior = TreePriorType.CONSTANT;
    private TreePriorParameterizationType parameterization = TreePriorParameterizationType.GROWTH_RATE;
    private int skylineGroupCount = 10;
    private TreePriorParameterizationType skylineModel = TreePriorParameterizationType.CONSTANT_SKYLINE;
    private TreePriorParameterizationType skyrideSmoothing = TreePriorParameterizationType.TIME_AWARE_SKYRIDE;
    private int skyGridCount = 50;
    private double skyGridInterval = Double.NaN;
    private VariableDemographicModel.Type extendedSkylineModel = VariableDemographicModel.Type.LINEAR;
    private double birthDeathSamplingProportion = 1.0;
    private PopulationSizeModelType populationSizeModel = PopulationSizeModelType.CONTINUOUS_CONSTANT;
    private CalibrationPoints.CorrectionType calibCorrectionType = CalibrationPoints.CorrectionType.EXACT;
    private boolean fixedTree = false;

    public PartitionTreePrior(BeautiOptions options, PartitionTreeModel treeModel) {
        super(options, treeModel.getName());

        initModelParametersAndOpererators();
    }

    /**
     * A copy constructor
     *
     * @param options the beauti options
     * @param name    the name of the new model
     * @param source  the source model
     */
    public PartitionTreePrior(BeautiOptions options, String name, PartitionTreePrior source) {
        super(options, name);

        this.nodeHeightPrior = source.nodeHeightPrior;
        this.parameterization = source.parameterization;
        this.skylineGroupCount = source.skylineGroupCount;
        this.skylineModel = source.skylineModel;
        this.skyrideSmoothing = source.skyrideSmoothing;
        this.extendedSkylineModel = source.extendedSkylineModel;
        this.birthDeathSamplingProportion = source.birthDeathSamplingProportion;
        this.populationSizeModel = source.populationSizeModel;
        this.calibCorrectionType = source.calibCorrectionType;
        this.fixedTree = source.fixedTree;

        initModelParametersAndOpererators();
    }

    @Override
    public void initModelParametersAndOpererators() {

        createParameterOneOverXPrior("constant.popSize", "coalescent population size parameter",
                PriorScaleType.TIME_SCALE, 1.0);

        createParameterOneOverXPrior("exponential.popSize", "coalescent population size parameter",
                PriorScaleType.TIME_SCALE, 1.0);
        createParameterLaplacePrior("exponential.growthRate", "coalescent growth rate parameter",
                PriorScaleType.GROWTH_RATE_SCALE, 0.0, 0.0, 1.0);
        createParameterGammaPrior("exponential.doublingTime", "coalescent doubling time parameter",
                PriorScaleType.NONE, 100.0, 0.001, 1000, true);

        createParameterOneOverXPrior("logistic.popSize", "coalescent population size parameter",
                PriorScaleType.TIME_SCALE, 1.0);
        createParameterLaplacePrior("logistic.growthRate", "coalescent logistic growth rate parameter",
                PriorScaleType.GROWTH_RATE_SCALE, 0.0, 0.0, 1.0);
        createParameterGammaPrior("logistic.doublingTime", "coalescent doubling time parameter",
                PriorScaleType.NONE, 100.0, 0.001, 1000, true);
        createParameterGammaPrior("logistic.t50", "logistic shape parameter",
                PriorScaleType.NONE, 1.0, 0.001, 1000, true);

        createParameterOneOverXPrior("expansion.popSize", "coalescent population size parameter",
                PriorScaleType.TIME_SCALE, 1.0);
        createParameterLaplacePrior("expansion.growthRate", "coalescent expansion growth rate parameter",
                PriorScaleType.GROWTH_RATE_SCALE, 0.0, 0.0, 1.0);
        createParameterGammaPrior("expansion.doublingTime", "coalescent doubling time parameter",
                PriorScaleType.NONE, 100.0, 0.001, 1000, true);
        createZeroOneParameterUniformPrior("expansion.ancestralProportion", "ancestral population proportion", 0.1);

        createNonNegativeParameterUniformPrior("skyline.popSize", "Bayesian Skyline population sizes",
                PriorScaleType.TIME_SCALE, 1.0, 0.0, Parameter.UNIFORM_MAX_BOUND);
        createParameter("skyline.groupSize", "Bayesian Skyline group sizes");
        // skyride.logPopSize is log unit unlike other popSize
        createParameterUniformPrior("skyride.logPopSize", "GMRF Bayesian skyride population sizes (log unit)",
                PriorScaleType.LOG_TIME_SCALE, 1.0, -Parameter.UNIFORM_MAX_BOUND, Parameter.UNIFORM_MAX_BOUND);
        createParameter("skyride.groupSize", "GMRF Bayesian skyride group sizes (for backward compatibility)");
        createParameterGammaPrior("skyride.precision", "GMRF Bayesian skyride precision",
                PriorScaleType.NONE, 1.0, 0.001, 1000, true);

        createParameterUniformPrior("skygrid.logPopSize", "GMRF Bayesian SkyGrid population sizes (log unit)",
                PriorScaleType.LOG_TIME_SCALE, 1.0, -Parameter.UNIFORM_MAX_BOUND, Parameter.UNIFORM_MAX_BOUND);
        createParameterGammaPrior("skygrid.precision", "GMRF Bayesian SkyGrid precision",
                PriorScaleType.NONE, 0.1, 0.001, 1000, true);
        createParameterUniformPrior("skygrid.numGridPoints", "GMRF Bayesian SkyGrid number of grid points)",
                PriorScaleType.NONE, 1.0, -Parameter.UNIFORM_MAX_BOUND, Parameter.UNIFORM_MAX_BOUND);
        createParameterUniformPrior("skygrid.cutOff", "GMRF Bayesian SkyGrid cut-off time",
                PriorScaleType.TIME_SCALE, 1.0, 0.0, Parameter.UNIFORM_MAX_BOUND);

        createNonNegativeParameterUniformPrior("demographic.popSize", "Extended Bayesian Skyline population sizes",
                PriorScaleType.TIME_SCALE, 1.0, 0.0, Parameter.UNIFORM_MAX_BOUND);
        createParameter("demographic.indicators", "Extended Bayesian Skyline population switch", 0.0);
        createParameterOneOverXPrior("demographic.populationMean", "Extended Bayesian Skyline population prior mean",
                PriorScaleType.TIME_SCALE, 1);

        createDiscreteStatistic("demographic.populationSizeChanges", "Average number of population change points"); // POISSON_PRIOR

        createNonNegativeParameterUniformPrior("yule.birthRate", "Yule speciation process birth rate",
                PriorScaleType.BIRTH_RATE_SCALE, 1.0, 0.0, Parameter.UNIFORM_MAX_BOUND);

        createNonNegativeParameterUniformPrior(BirthDeathModelParser.MEAN_GROWTH_RATE_PARAM_NAME, "Birth-Death speciation process rate",
                PriorScaleType.BIRTH_RATE_SCALE, 0.01, 0.0, 100000.0);
        createNonNegativeParameterUniformPrior(BirthDeathModelParser.RELATIVE_DEATH_RATE_PARAM_NAME, "Birth-Death speciation process relative death rate",
                PriorScaleType.NONE, 0.5, 0.0, 1.0);
        createParameterBetaDistributionPrior(BirthDeathModelParser.BIRTH_DEATH + "." + BirthDeathModelParser.SAMPLE_PROB,
                "Birth-Death the proportion of taxa sampled from birth-death tree",
                0.01, 1.0, 1.0, 0.0);
        createNonNegativeParameterUniformPrior(BirthDeathSerialSamplingModelParser.BDSS + "."
                        + BirthDeathSerialSamplingModelParser.LAMBDA,
                "Birth-Death speciation process rate", PriorScaleType.BIRTH_RATE_SCALE,
                2.0, 0.0, 100000.0);
        createZeroOneParameterUniformPrior(BirthDeathSerialSamplingModelParser.BDSS + "."
                        + BirthDeathSerialSamplingModelParser.RELATIVE_MU,
                "Birth-Death relative death rate", 0.5);
        //Issue 656
//        createParameterBetaDistributionPrior(BirthDeathSerialSamplingModelParser.BDSS + "."
//                + BirthDeathSerialSamplingModelParser.SAMPLE_PROBABILITY,
//                "Birth-Death the proportion of taxa sampled from birth death tree",
//                0.01, 1.0, 1.0, 0.0);
        createNonNegativeParameterUniformPrior(BirthDeathSerialSamplingModelParser.BDSS + "."
                        + BirthDeathSerialSamplingModelParser.PSI,
                "Birth-Death rate of sampling taxa through time", PriorScaleType.NONE,
                0.05, 0.0, 100.0);
        createNonNegativeParameterUniformPrior(BirthDeathSerialSamplingModelParser.BDSS + "."
                        + BirthDeathSerialSamplingModelParser.ORIGIN,
                "Birth-Death the time of the lineage originated (must > root height)", PriorScaleType.ORIGIN_SCALE,
                1.0, 0.0, Parameter.UNIFORM_MAX_BOUND);
//        createParameter(BirthDeathSerialSamplingModelParser.BDSS + "." + BirthDeathSerialSamplingModelParser.R,
//                "Birth-Death the probabilty that a sampled individual continues being infectious after sample event",
//                1.0); // fixed to 1
//        createNonNegativeParameterUniformPrior(BirthDeathSerialSamplingModelParser.BDSS + "."
//                + BirthDeathSerialSamplingModelParser.HAS_FINAL_SAMPLE,
//                "Birth-Death the time in the past when the process starts with the first individual", PriorScaleType.NONE,
//                80.0, 0.0, 1000.0);
        createNonNegativeParameterUniformPrior(BirthDeathEpidemiologyModelParser.ORIGIN,
                "The origin of the infection, x0 > tree.rootHeight", PriorScaleType.ORIGIN_SCALE,
                1.0, 0.0, Parameter.UNIFORM_MAX_BOUND);
        createParameterLognormalPrior(BirthDeathEpidemiologyModelParser.R0, "R0",
                PriorScaleType.NONE, 2.0, 1.0, 1.25, 0.0);
        createNonNegativeParameterUniformPrior(BirthDeathEpidemiologyModelParser.RECOVERY_RATE,
                "recoveryRate", PriorScaleType.NONE,
                0.05, 0.0, 100.0);
        createParameterBetaDistributionPrior(BirthDeathEpidemiologyModelParser.SAMPLING_PROBABILITY,
                "samplingProbability",
                0.01, 1.0, 1.0, 0.0);

        createScaleOperator("constant.popSize", demoTuning, demoWeights);
        createScaleOperator("exponential.popSize", demoTuning, demoWeights);
        createOperator("exponential.growthRate", OperatorType.RANDOM_WALK, 1.0, demoWeights);
        createScaleOperator("exponential.doublingTime", demoTuning, demoWeights);
        createScaleOperator("logistic.popSize", demoTuning, demoWeights);
        createOperator("logistic.growthRate", OperatorType.RANDOM_WALK, 1.0, demoWeights);
//        createScaleOperator("logistic.growthRate", demoTuning, demoWeights);
        createScaleOperator("logistic.doublingTime", demoTuning, demoWeights);
        createScaleOperator("logistic.t50", demoTuning, demoWeights);
        createScaleOperator("expansion.popSize", demoTuning, demoWeights);
        createOperator("expansion.growthRate", OperatorType.RANDOM_WALK, 1.0, demoWeights);
//        createScaleOperator("expansion.growthRate", demoTuning, demoWeights);
        createScaleOperator("expansion.doublingTime", demoTuning, demoWeights);
        createScaleOperator("expansion.ancestralProportion", demoTuning, demoWeights);
        createScaleOperator("skyline.popSize", demoTuning, demoWeights * 5);
        createOperator("skyline.groupSize", OperatorType.INTEGER_DELTA_EXCHANGE, 1.0, demoWeights * 2);
        createOperator("demographic.populationMean", OperatorType.SCALE, 0.9, demoWeights);
        createOperator("demographic.indicators", OperatorType.BITFLIP, 1, 2 * treeWeights);

        // hack pass distribution in name
        createOperatorUsing2Parameters("demographic.popSize", "demographic.populationMeanDist", "", "demographic.popSize",
                "demographic.indicators", OperatorType.SAMPLE_NONACTIVE, 1, 5 * demoWeights);
        createOperatorUsing2Parameters("demographic.scaleActive", "demographic.scaleActive", "", "demographic.popSize",
                "demographic.indicators", OperatorType.SCALE_WITH_INDICATORS, 0.5, 2 * demoWeights);
        createOperatorUsing2Parameters("gmrfGibbsOperator", "gmrfGibbsOperator", "Gibbs sampler for GMRF Skyride", "skyride.logPopSize",
                "skyride.precision", OperatorType.GMRF_GIBBS_OPERATOR, 2, 2);
        createOperatorUsing2Parameters("gmrfSkyGridGibbsOperator", "gmrfGibbsOperator", "Gibbs sampler for Bayesian SkyGrid", "skygrid.logPopSize",
                "skygrid.precision", OperatorType.SKY_GRID_GIBBS_OPERATOR, 1.0, 2);
        createScaleOperator("skygrid.precision", "description", 0.75, 1.0);

        createScaleOperator("yule.birthRate", demoTuning, demoWeights);

        createScaleOperator(BirthDeathModelParser.MEAN_GROWTH_RATE_PARAM_NAME, demoTuning, demoWeights);
        createScaleOperator(BirthDeathModelParser.RELATIVE_DEATH_RATE_PARAM_NAME, demoTuning, demoWeights);
        createScaleOperator(BirthDeathModelParser.BIRTH_DEATH + "." + BirthDeathModelParser.SAMPLE_PROB, demoTuning, demoWeights);
        createScaleOperator(BirthDeathSerialSamplingModelParser.BDSS + "."
                + BirthDeathSerialSamplingModelParser.LAMBDA, demoTuning, 1);
        createScaleOperator(BirthDeathSerialSamplingModelParser.BDSS + "."
                + BirthDeathSerialSamplingModelParser.RELATIVE_MU, demoTuning, 1);
        //Issue 656
//        createScaleOperator(BirthDeathSerialSamplingModelParser.BDSS + "."
//                + BirthDeathSerialSamplingModelParser.SAMPLE_PROBABILITY, demoTuning, 1);
        createScaleOperator(BirthDeathSerialSamplingModelParser.BDSS + "."
                + BirthDeathSerialSamplingModelParser.PSI, demoTuning, 1);   // todo random worl op ?
        createScaleOperator(BirthDeathSerialSamplingModelParser.BDSS + "."
                + BirthDeathSerialSamplingModelParser.ORIGIN, demoTuning, 1);
//        createScaleOperator(BirthDeathSerialSamplingModelParser.BDSS + "."
//                + BirthDeathSerialSamplingModelParser.R, demoTuning, 1);
//        createScaleOperator(BirthDeathSerialSamplingModelParser.BDSS + "."
//                + BirthDeathSerialSamplingModelParser.HAS_FINAL_SAMPLE, demoTuning, 1);
//        createOperator(BirthDeathSerialSamplingModelParser.BDSS + "."
//                + BirthDeathSerialSamplingModelParser.SAMPLE_BECOMES_NON_INFECTIOUS,
//                OperatorType.RANDOM_WALK, 1.0, demoWeights);
        createScaleOperator(BirthDeathEpidemiologyModelParser.ORIGIN, demoTuning, 1);
        createScaleOperator(BirthDeathEpidemiologyModelParser.R0, demoTuning, 1);
        createScaleOperator(BirthDeathEpidemiologyModelParser.RECOVERY_RATE, demoTuning, 1);
        createScaleOperator(BirthDeathEpidemiologyModelParser.SAMPLING_PROBABILITY, demoTuning, 1);
    }

    @Override
    public List<Parameter> selectParameters(List<Parameter> params) {
//        setAvgRootAndRate();

        if (nodeHeightPrior == TreePriorType.CONSTANT) {
            params.add(getParameter("constant.popSize"));
        } else if (nodeHeightPrior == TreePriorType.EXPONENTIAL) {
            params.add(getParameter("exponential.popSize"));
            if (parameterization == TreePriorParameterizationType.GROWTH_RATE) {
                params.add(getParameter("exponential.growthRate"));
            } else {
                params.add(getParameter("exponential.doublingTime"));
            }
        } else if (nodeHeightPrior == TreePriorType.LOGISTIC) {
            params.add(getParameter("logistic.popSize"));
            if (parameterization == TreePriorParameterizationType.GROWTH_RATE) {
                params.add(getParameter("logistic.growthRate"));
            } else {
                params.add(getParameter("logistic.doublingTime"));
            }
            params.add(getParameter("logistic.t50"));
        } else if (nodeHeightPrior == TreePriorType.EXPANSION) {
            params.add(getParameter("expansion.popSize"));
            if (parameterization == TreePriorParameterizationType.GROWTH_RATE) {
                params.add(getParameter("expansion.growthRate"));
            } else {
                params.add(getParameter("expansion.doublingTime"));
            }
            params.add(getParameter("expansion.ancestralProportion"));
        } else if (nodeHeightPrior == TreePriorType.SKYLINE) {
            params.add(getParameter("skyline.popSize"));
        } else if (nodeHeightPrior == TreePriorType.EXTENDED_SKYLINE) {
            params.add(getParameter("demographic.populationSizeChanges"));
            params.add(getParameter("demographic.populationMean"));
        } else if (nodeHeightPrior == TreePriorType.GMRF_SKYRIDE) {
//            params.add(getParameter("skyride.popSize")); // force user to use GMRF, not allowed to change
            params.add(getParameter("skyride.precision"));
        } else if (nodeHeightPrior == TreePriorType.SKYGRID) {
//            params.add(getParameter("skyride.popSize")); // force user to use GMRF, not allowed to change
            params.add(getParameter("skygrid.precision"));
        } else if (nodeHeightPrior == TreePriorType.YULE || nodeHeightPrior == TreePriorType.YULE_CALIBRATION) {
            params.add(getParameter("yule.birthRate"));
        } else if (nodeHeightPrior == TreePriorType.BIRTH_DEATH || nodeHeightPrior == TreePriorType.BIRTH_DEATH_INCOMPLETE_SAMPLING) {
            params.add(getParameter(BirthDeathModelParser.MEAN_GROWTH_RATE_PARAM_NAME));
            params.add(getParameter(BirthDeathModelParser.RELATIVE_DEATH_RATE_PARAM_NAME));
            if (nodeHeightPrior == TreePriorType.BIRTH_DEATH_INCOMPLETE_SAMPLING)
                params.add(getParameter(BirthDeathModelParser.BIRTH_DEATH + "." + BirthDeathModelParser.SAMPLE_PROB));

        } else if (nodeHeightPrior == TreePriorType.BIRTH_DEATH_SERIAL_SAMPLING) {
            params.add(getParameter(BirthDeathSerialSamplingModelParser.BDSS + "."
                    + BirthDeathSerialSamplingModelParser.LAMBDA));
            params.add(getParameter(BirthDeathSerialSamplingModelParser.BDSS + "."
                    + BirthDeathSerialSamplingModelParser.RELATIVE_MU));
            Parameter psi = getParameter(BirthDeathSerialSamplingModelParser.BDSS + "."
                    + BirthDeathSerialSamplingModelParser.PSI);
            if (options.maximumTipHeight > 0) {
                psi.setInitial(MathUtils.round(1 / options.maximumTipHeight, 4));
            }
            params.add(psi);
            params.add(getParameter(BirthDeathSerialSamplingModelParser.BDSS + "."
                    + BirthDeathSerialSamplingModelParser.ORIGIN));
//            params.add(getParameter(BirthDeathSerialSamplingModelParser.BDSS + "."
//                    + BirthDeathSerialSamplingModelParser.SAMPLE_PROBABILITY));

        } else if (nodeHeightPrior == TreePriorType.BIRTH_DEATH_BASIC_REPRODUCTIVE_NUMBER) {
            params.add(getParameter(BirthDeathEpidemiologyModelParser.ORIGIN));
            params.add(getParameter(BirthDeathEpidemiologyModelParser.R0));
            params.add(getParameter(BirthDeathEpidemiologyModelParser.RECOVERY_RATE));
            params.add(getParameter(BirthDeathEpidemiologyModelParser.SAMPLING_PROBABILITY));

        }
        return params;
    }

    @Override
    public List<Operator> selectOperators(List<Operator> ops) {

        if (nodeHeightPrior == TreePriorType.CONSTANT) {
            ops.add(getOperator("constant.popSize"));
        } else if (nodeHeightPrior == TreePriorType.EXPONENTIAL) {
            ops.add(getOperator("exponential.popSize"));
            if (parameterization == TreePriorParameterizationType.GROWTH_RATE) {
                ops.add(getOperator("exponential.growthRate"));
            } else {
                ops.add(getOperator("exponential.doublingTime"));
            }
        } else if (nodeHeightPrior == TreePriorType.LOGISTIC) {
            ops.add(getOperator("logistic.popSize"));
            if (parameterization == TreePriorParameterizationType.GROWTH_RATE) {
                ops.add(getOperator("logistic.growthRate"));
            } else {
                ops.add(getOperator("logistic.doublingTime"));
            }
            ops.add(getOperator("logistic.t50"));
        } else if (nodeHeightPrior == TreePriorType.EXPANSION) {
            ops.add(getOperator("expansion.popSize"));
            if (parameterization == TreePriorParameterizationType.GROWTH_RATE) {
                ops.add(getOperator("expansion.growthRate"));
            } else {
                ops.add(getOperator("expansion.doublingTime"));
            }
            ops.add(getOperator("expansion.ancestralProportion"));
        } else if (nodeHeightPrior == TreePriorType.SKYLINE) {
            ops.add(getOperator("skyline.popSize"));
            ops.add(getOperator("skyline.groupSize"));
        } else if (nodeHeightPrior == TreePriorType.GMRF_SKYRIDE) {
            ops.add(getOperator("gmrfGibbsOperator"));
        } else if (nodeHeightPrior == TreePriorType.SKYGRID) {
            ops.add(getOperator("gmrfSkyGridGibbsOperator"));
            ops.add(getOperator("skygrid.precision"));
        } else if (nodeHeightPrior == TreePriorType.EXTENDED_SKYLINE) {
            ops.add(getOperator("demographic.populationMean"));
            ops.add(getOperator("demographic.popSize"));
            ops.add(getOperator("demographic.indicators"));
            ops.add(getOperator("demographic.scaleActive"));
        } else if (nodeHeightPrior == TreePriorType.YULE || nodeHeightPrior == TreePriorType.YULE_CALIBRATION) {
            ops.add(getOperator("yule.birthRate"));
        } else if (nodeHeightPrior == TreePriorType.BIRTH_DEATH || nodeHeightPrior == TreePriorType.BIRTH_DEATH_INCOMPLETE_SAMPLING) {
            ops.add(getOperator(BirthDeathModelParser.MEAN_GROWTH_RATE_PARAM_NAME));
            ops.add(getOperator(BirthDeathModelParser.RELATIVE_DEATH_RATE_PARAM_NAME));
            if (nodeHeightPrior == TreePriorType.BIRTH_DEATH_INCOMPLETE_SAMPLING)
                ops.add(getOperator(BirthDeathModelParser.BIRTH_DEATH + "." + BirthDeathModelParser.SAMPLE_PROB));
        } else if (nodeHeightPrior == TreePriorType.BIRTH_DEATH_SERIAL_SAMPLING) {
            ops.add(getOperator(BirthDeathSerialSamplingModelParser.BDSS + "."
                    + BirthDeathSerialSamplingModelParser.LAMBDA));
            ops.add(getOperator(BirthDeathSerialSamplingModelParser.BDSS + "."
                    + BirthDeathSerialSamplingModelParser.RELATIVE_MU));
            ops.add(getOperator(BirthDeathSerialSamplingModelParser.BDSS + "."
                    + BirthDeathSerialSamplingModelParser.PSI));
            ops.add(getOperator(BirthDeathSerialSamplingModelParser.BDSS + "."
                    + BirthDeathSerialSamplingModelParser.ORIGIN));
//            ops.add(getOperator(BirthDeathSerialSamplingModelParser.BDSS + "."
//                    + BirthDeathSerialSamplingModelParser.SAMPLE_PROBABILITY));

//            if (nodeHeightPrior == TreePriorType.BIRTH_DEATH_BASIC_REPRODUCTIVE_NUMBER) {
//                ops.add(getOperator(BirthDeathSerialSamplingModelParser.BDSS + "."
//                + BirthDeathSerialSamplingModelParser.R));
//                ops.add(getOperator(BirthDeathSerialSamplingModelParser.BDSS + "."
//                + BirthDeathSerialSamplingModelParser.HAS_FINAL_SAMPLE));
//            }
        } else if (nodeHeightPrior == TreePriorType.BIRTH_DEATH_BASIC_REPRODUCTIVE_NUMBER) {
            ops.add(getOperator(BirthDeathEpidemiologyModelParser.ORIGIN));
            ops.add(getOperator(BirthDeathEpidemiologyModelParser.R0));
            ops.add(getOperator(BirthDeathEpidemiologyModelParser.RECOVERY_RATE));
            ops.add(getOperator(BirthDeathEpidemiologyModelParser.SAMPLING_PROBABILITY));
        }
        return ops;
    }


    //////////////////////////////////////////////////////
    public String getPrefix() {
        String prefix = "";
        if (options.getPartitionTreePriors().size() > 1) {//|| options.isSpeciesAnalysis()
            // There is more than one active partition model, or doing species analysis
            prefix += getName() + ".";
        }
        return prefix;
    }

    /////////////////////////////////////////////////////////////////////////

    public TreePriorType getNodeHeightPrior() {
        return nodeHeightPrior;
    }

    public void setNodeHeightPrior(TreePriorType nodeHeightPrior) {
        this.nodeHeightPrior = nodeHeightPrior;
    }

    public TreePriorParameterizationType getParameterization() {
        return parameterization;
    }

    public void setParameterization(TreePriorParameterizationType parameterization) {
        this.parameterization = parameterization;
    }

    public int getSkylineGroupCount() {
        return skylineGroupCount;
    }

    public void setSkylineGroupCount(int skylineGroupCount) {
        this.skylineGroupCount = skylineGroupCount;
    }

    public TreePriorParameterizationType getSkylineModel() {
        return skylineModel;
    }

    public void setSkylineModel(TreePriorParameterizationType skylineModel) {
        this.skylineModel = skylineModel;
    }

    public TreePriorParameterizationType getSkyrideSmoothing() {
        return skyrideSmoothing;
    }

    public void setSkyrideSmoothing(TreePriorParameterizationType skyrideSmoothing) {
        this.skyrideSmoothing = skyrideSmoothing;
    }

    public int getSkyGridCount() {
        return skyGridCount;
    }

    public void setSkyGridCount(int count) {
        this.skyGridCount = count;
    }

    public double getSkyGridInterval() {
        return skyGridInterval;
    }

    public void setSkyGridInterval(double x) {
        this.skyGridInterval = x;
    }

    public double getBirthDeathSamplingProportion() {
        return birthDeathSamplingProportion;
    }

    public void setBirthDeathSamplingProportion(double birthDeathSamplingProportion) {
        this.birthDeathSamplingProportion = birthDeathSamplingProportion;
    }

    public boolean isFixedTree() {
        return fixedTree;
    }

    public void setFixedTree(boolean fixedTree) {
        this.fixedTree = fixedTree;
    }

    public void setExtendedSkylineModel(VariableDemographicModel.Type extendedSkylineModel) {
        this.extendedSkylineModel = extendedSkylineModel;
    }

    public VariableDemographicModel.Type getExtendedSkylineModel() {
        return extendedSkylineModel;
    }

    public PopulationSizeModelType getPopulationSizeModel() {
        return populationSizeModel;
    }

    public void setPopulationSizeModel(PopulationSizeModelType populationSizeModel) {
        this.populationSizeModel = populationSizeModel;
    }

    public CalibrationPoints.CorrectionType getCalibCorrectionType() {
        return calibCorrectionType;
    }

    public void setCalibCorrectionType(CalibrationPoints.CorrectionType calibCorrectionType) {
        this.calibCorrectionType = calibCorrectionType;
    }

    public BeautiOptions getOptions() {
        return options;
    }

}