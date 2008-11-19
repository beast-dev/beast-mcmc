/*
 * BeautiOptions.java
 *
 * Copyright (C) 2002-2007 Alexei Drummond and Andrew Rambaut
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

import dr.app.beauti.priorsPanel.PriorType;
import dr.evolution.datatype.DataType;
import dr.evolution.tree.Tree;
import dr.evolution.util.Date;
import dr.evolution.util.Taxa;
import dr.evolution.util.TaxonList;
import dr.evolution.util.Units;
import dr.evomodel.coalescent.VariableDemographicModel;
import dr.evomodelxml.BirthDeathModelParser;
import dr.inference.operators.OperatorSchedule;
import dr.util.NumberFormatter;
import org.jdom.Document;
import org.jdom.Element;

import java.text.NumberFormat;
import java.text.ParseException;
import java.util.*;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 */
public class BeautiOptions extends ModelOptions {

    public BeautiOptions() {
        double demoWeights = 3.0;
        double branchWeights = 30.0;
        double treeWeights = 15.0;
        double rateWeights = 3.0;

        createParameter("tree", "The tree");
        createParameter("treeModel.internalNodeHeights", "internal node heights of the tree (except the root)");
        createParameter("treeModel.allInternalNodeHeights", "internal node heights of the tree");
        createParameter("treeModel.rootHeight", "root height of the tree", true, 1.0, 0.0, Double.POSITIVE_INFINITY);

        // A vector of relative rates across all partitions...
        createParameter("allMus", "All the relative rates");

        createParameter("clock.rate", "substitution rate", SUBSTITUTION_RATE_SCALE, 1.0, 0.0, Double.POSITIVE_INFINITY);
        createParameter("uced.mean", "uncorrelated exponential relaxed clock mean", SUBSTITUTION_RATE_SCALE, 1.0, 0.0, Double.POSITIVE_INFINITY);
        createParameter("ucld.mean", "uncorrelated lognormal relaxed clock mean", SUBSTITUTION_RATE_SCALE, 1.0, 0.0, Double.POSITIVE_INFINITY);
        createParameter("ucld.stdev", "uncorrelated lognormal relaxed clock stdev", LOG_STDEV_SCALE, 0.1, 0.0, Double.POSITIVE_INFINITY);
        createParameter("branchRates.categories", "relaxed clock branch rate categories");
        createParameter("localClock.rates", "random local clock rates", SUBSTITUTION_RATE_SCALE, 1.0, 0.0, Double.POSITIVE_INFINITY);
        createParameter("localClock.changes", "random local clock rate change indicator");

        {
            final Parameter p = createParameter("treeModel.rootRate", "autocorrelated lognormal relaxed clock root rate", ROOT_RATE_SCALE, 1.0, 0.0, Double.POSITIVE_INFINITY);
            p.priorType = PriorType.GAMMA_PRIOR;
            p.gammaAlpha = 1;
            p.gammaBeta = 0.0001;
        }
        createParameter("treeModel.nodeRates", "autocorrelated lognormal relaxed clock non-root rates", SUBSTITUTION_RATE_SCALE, 1.0, 0.0, Double.POSITIVE_INFINITY);
        createParameter("treeModel.allRates", "autocorrelated lognormal relaxed clock all rates", SUBSTITUTION_RATE_SCALE, 1.0, 0.0, Double.POSITIVE_INFINITY);
        {
            final Parameter p = createParameter("branchRates.var", "autocorrelated lognormal relaxed clock rate variance ", LOG_VAR_SCALE, 0.1, 0.0, Double.POSITIVE_INFINITY);
            p.priorType = PriorType.GAMMA_PRIOR;
            p.gammaAlpha = 1;
            p.gammaBeta = 0.0001;
        }

        createParameter("errorModel.ageRate", "age dependent sequence error rate", SUBSTITUTION_RATE_SCALE, 1.0E-8, 0.0, Double.POSITIVE_INFINITY);
        createParameter("errorModel.baseRate", "base sequence error rate", SUBSTITUTION_RATE_SCALE, 1.0E-8, 0.0, Double.POSITIVE_INFINITY);

        createScaleParameter("constant.popSize", "coalescent population size parameter", TIME_SCALE, 1.0, 0.0, Double.POSITIVE_INFINITY);

        createScaleParameter("exponential.popSize", "coalescent population size parameter", TIME_SCALE, 1.0, 0.0, Double.POSITIVE_INFINITY);
        createParameter("exponential.growthRate", "coalescent growth rate parameter", GROWTH_RATE_SCALE, 0.0, Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY);
        createParameter("exponential.doublingTime", "coalescent doubling time parameter", TIME_SCALE, 0.5, 0.0, Double.POSITIVE_INFINITY);
        createScaleParameter("logistic.popSize", "coalescent population size parameter", TIME_SCALE, 1.0, 0.0, Double.POSITIVE_INFINITY);
        createParameter("logistic.growthRate", "coalescent logistic growth rate parameter", GROWTH_RATE_SCALE, 0.001, 0.0, Double.POSITIVE_INFINITY);
        createParameter("logistic.doublingTime", "coalescent doubling time parameter", TIME_SCALE, 0.5, 0.0, Double.POSITIVE_INFINITY);
        createParameter("logistic.t50", "logistic shape parameter", T50_SCALE, 0.1, 0.0, Double.POSITIVE_INFINITY);
        createScaleParameter("expansion.popSize", "coalescent population size parameter", TIME_SCALE, 1.0, 0.0, Double.POSITIVE_INFINITY);
        createParameter("expansion.growthRate", "coalescent logistic growth rate parameter", GROWTH_RATE_SCALE, 0.001, 0.0, Double.POSITIVE_INFINITY);
        createParameter("expansion.doublingTime", "coalescent doubling time parameter", TIME_SCALE, 0.5, 0.0, Double.POSITIVE_INFINITY);
        createParameter("expansion.ancestralProportion", "ancestral population proportion", NONE, 0.1, 0.0, 1.0);
        createParameter("skyline.popSize", "Bayesian Skyline population sizes", TIME_SCALE, 1.0, 0.0, Double.POSITIVE_INFINITY);
        createParameter("skyline.groupSize", "Bayesian Skyline group sizes");

        createParameter("skyride.popSize", "GMRF Bayesian skyride population sizes", TIME_SCALE, 1.0, Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY);
        createParameter("skyride.groupSize", "GMRF Bayesian skyride group sizes (for backward compatibility)");
        {
            final Parameter p = createParameter("skyride.precision", "GMRF Bayesian skyride precision", NONE, 1.0, 0.0, Double.POSITIVE_INFINITY);
            p.priorType = PriorType.GAMMA_PRIOR;
            p.gammaAlpha = 0.001;
            p.gammaBeta = 1000;
            p.priorFixed = true;
        }

        createParameter("demographic.popSize", "Extended Bayesian Skyline population sizes", TIME_SCALE, 1.0, 0.0, Double.POSITIVE_INFINITY);
        createParameter("demographic.indicators", "Extended Bayesian Skyline population switch");
        createScaleParameter("demographic.populationMean", "Extended Bayesian Skyline population prior mean", TIME_SCALE, 1, 0, Double.POSITIVE_INFINITY);
        {
            final Parameter p = createStatistic("demographic.populationSizeChanges", "Average number of population change points", true);
            p.priorType = PriorType.POISSON_PRIOR;
            p.poissonMean = Math.log(2);
        }
        createParameter("yule.birthRate", "Yule speciation process birth rate", BIRTH_RATE_SCALE, 1.0, 0.0, Double.POSITIVE_INFINITY);

        createParameter(BirthDeathModelParser.BIRTHDIFF_RATE_PARAM_NAME, "Birth-Death speciation process rate", BIRTH_RATE_SCALE, 1.0, 0.0, Double.POSITIVE_INFINITY);
        createParameter(BirthDeathModelParser.RELATIVE_DEATH_RATE_PARAM_NAME, "Death/Birth speciation process relative death rate", BIRTH_RATE_SCALE, 0.5, 0.0, 1.0);

        createScaleOperator("constant.popSize", demoWeights);
        createScaleOperator("exponential.popSize", demoWeights);
        createOperator("exponential.growthRate", OperatorType.RANDOM_WALK, 1.0, demoWeights);
        createScaleOperator("exponential.doublingTime", demoWeights);
        createScaleOperator("logistic.popSize", demoWeights);
        createScaleOperator("logistic.growthRate", demoWeights);
        createScaleOperator("logistic.doublingTime", demoWeights);
        createScaleOperator("logistic.t50", demoWeights);
        createScaleOperator("expansion.popSize", demoWeights);
        createScaleOperator("expansion.growthRate", demoWeights);
        createScaleOperator("expansion.doublingTime", demoWeights);
        createScaleOperator("expansion.ancestralProportion", demoWeights);
        createScaleOperator("skyline.popSize", demoWeights * 5);
        createOperator("skyline.groupSize", OperatorType.INTEGER_DELTA_EXCHANGE, 1.0, demoWeights * 2);

        createOperator("demographic.populationMean", OperatorType.SCALE, 0.9, demoWeights);
        createOperator("demographic.indicators", OperatorType.BITFLIP, 1, 2 * treeWeights);
        // hack pass distribution in name
        createOperator("demographic.popSize", "demographic.populationMeanDist", "", "demographic.popSize",
                "demographic.indicators", OperatorType.SAMPLE_NONACTIVE, 1, 5 * demoWeights);
        createOperator("demographic.scaleActive", "demographic.scaleActive", "", "demographic.popSize",
                "demographic.indicators", OperatorType.SCALE_WITH_INDICATORS, 0.5, 2 * demoWeights);

        createOperator("gmrfGibbsOperator", "gmrfGibbsOperator", "Gibbs sampler for GMRF", "skyride.popSize",
                "skyride.precision", OperatorType.GMRF_GIBBS_OPERATOR, 0.5, 2);

        createScaleOperator("yule.birthRate", demoWeights);

        createScaleOperator(BirthDeathModelParser.BIRTHDIFF_RATE_PARAM_NAME, demoWeights);
        createScaleOperator(BirthDeathModelParser.RELATIVE_DEATH_RATE_PARAM_NAME, demoWeights);

        // These are statistics which could have priors on...
        createStatistic("meanRate", "The mean rate of evolution over the whole tree", 0.0, Double.POSITIVE_INFINITY);
        createStatistic("coefficientOfVariation", "The variation in rate of evolution over the whole tree",
                0.0, Double.POSITIVE_INFINITY);
        createStatistic("covariance",
                "The covariance in rates of evolution on each lineage with their ancestral lineages",
                Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY);

        // This only works if the partitions are of the same size...
//        createOperator("centeredMu", "Relative rates",
//                "Scales codon position rates relative to each other maintaining mean", "allMus",
//                OperatorType.CENTERED_SCALE, 0.75, rateWeights);
        createOperator("deltaMu", "Relative rates",
                "Changes partition relative rates relative to each other maintaining their mean", "allMus",
                OperatorType.DELTA_EXCHANGE, 0.75, rateWeights);

        createScaleOperator("errorModel.ageRate", rateWeights);
        createScaleOperator("errorModel.baseRate", rateWeights);

        createScaleOperator("clock.rate", rateWeights);
        createScaleOperator("uced.mean", rateWeights);
        createScaleOperator("ucld.mean", rateWeights);
        createScaleOperator("ucld.stdev", rateWeights);


        createOperator("scaleRootRate", "treeModel.rootRate",
                "Scales root rate", "treeModel.rootRate",
                OperatorType.SCALE, 0.75, rateWeights);
        createOperator("scaleOneRate", "treeModel.nodeRates",
                "Scales one non-root rate", "treeModel.nodeRates",
                OperatorType.SCALE, 0.75, branchWeights);
        createOperator("scaleAllRates", "treeModel.allRates",
                "Scales all rates simultaneously", "treeModel.allRates",
                OperatorType.SCALE_ALL, 0.75, rateWeights);
        createOperator("scaleAllRatesIndependently", "treeModel.nodeRates",
                "Scales all non-root rates independently", "treeModel.nodeRates",
                OperatorType.SCALE_INDEPENDENTLY, 0.75, rateWeights);

        createOperator("upDownAllRatesHeights", "All rates and heights",
                "Scales all rates inversely to node heights of the tree", "treeModel.allRates",
                "treeModel.allInternalNodeHeights", OperatorType.UP_DOWN, 0.75, branchWeights);
        createScaleOperator("branchRates.var", rateWeights);

        createOperator("swapBranchRateCategories", "branchRates.categories",
                "Performs a swap of branch rate categories", "branchRates.categories",
                OperatorType.SWAP, 1, branchWeights / 3);
        createOperator("randomWalkBranchRateCategories", "branchRates.categories",
                "Performs an integer random walk of branch rate categories", "branchRates.categories",
                OperatorType.INTEGER_RANDOM_WALK, 1, branchWeights / 3);
        createOperator("unformBranchRateCategories", "branchRates.categories",
                "Performs an integer uniform draw of branch rate categories", "branchRates.categories",
                OperatorType.INTEGER_UNIFORM, 1, branchWeights / 3);

        createScaleOperator("localClock.rates", treeWeights);
        createOperator("localClock.changes", OperatorType.BITFLIP, 1, treeWeights);
        createOperator("treeBitMove", "Tree", "Swaps the rates and change locations of local clocks", "tree",
                OperatorType.TREE_BIT_MOVE, -1.0, treeWeights);

        createScaleOperator("treeModel.rootHeight", demoWeights);
        createOperator("uniformHeights", "Internal node heights", "Draws new internal node heights uniformally",
                "treeModel.internalNodeHeights", OperatorType.UNIFORM, -1, branchWeights);

        createOperator("upDownRateHeights", "Substitution rate and heights",
                "Scales substitution rates inversely to node heights of the tree", "clock.rate",
                "treeModel.allInternalNodeHeights", OperatorType.UP_DOWN, 0.75, rateWeights);
        createOperator("upDownUCEDMeanHeights", "UCED mean and heights",
                "Scales UCED mean inversely to node heights of the tree", "uced.mean",
                "treeModel.allInternalNodeHeights", OperatorType.UP_DOWN, 0.75, rateWeights);
        createOperator("upDownUCLDMeanHeights", "UCLD mean and heights",
                "Scales UCLD mean inversely to node heights of the tree", "ucld.mean",
                "treeModel.allInternalNodeHeights", OperatorType.UP_DOWN, 0.75, rateWeights);

        createOperator("subtreeSlide", "Tree", "Performs the subtree-slide rearrangement of the tree", "tree",
                OperatorType.SUBTREE_SLIDE, 1.0, treeWeights);
        createOperator("narrowExchange", "Tree", "Performs local rearrangements of the tree", "tree",
                OperatorType.NARROW_EXCHANGE, -1, treeWeights);
        createOperator("wideExchange", "Tree", "Performs global rearrangements of the tree", "tree",
                OperatorType.WIDE_EXCHANGE, -1, demoWeights);
        createOperator("wilsonBalding", "Tree", "Performs the Wilson-Balding rearrangement of the tree", "tree",
                OperatorType.WILSON_BALDING, -1, demoWeights);
    }

    public void addPartitionModel(PartitionModel model) {

        if (!models.contains(model)) {
            models.add(model);

            // update delta mu opertor weight
            Operator deltaMuOperator = getOperator("deltaMu");
            deltaMuOperator.weight = getActiveModels().size();
        }
    }

    /**
     * @return a list of all partition models, whether or not they are used
     */
    public List<PartitionModel> getPartitionModels() {
        return models;
    }

    /**
     * @return a list of partition models that are of the given data type
     */
    public List<PartitionModel> getPartitionModels(DataType dataType) {
        List<PartitionModel> models = new ArrayList<PartitionModel>();
        for (PartitionModel model : getPartitionModels()) {
            if (model.dataType == dataType) {
                models.add(model);
            }
        }
        return models;
    }

    private double round(double value, int sf) {
        NumberFormatter formatter = new NumberFormatter(sf);
        try {
            return NumberFormat.getInstance().parse(formatter.format(value)).doubleValue();
        } catch (ParseException e) {
            return value;
        }
    }

    public Parameter getParameter(String name, PartitionModel model) {
        return model.getParameter(name);
    }

    /**
     * return an list of parameters that are required
     *
     * @return the parameter list
     */
    public ArrayList<Parameter> selectParameters() {

        ArrayList<Parameter> parameters = new ArrayList<Parameter>();

        selectParameters(parameters);
        selectStatistics(parameters);

        boolean multiplePartitions = getTotalActivePartitionCount() > 1;

        for (PartitionModel model : getActiveModels()) {
            parameters.addAll(model.getParameters(multiplePartitions));
        }


        double growthRateMaximum = 1E6;
        double birthRateMaximum = 1E6;
        double substitutionRateMaximum = 100;
        double logStdevMaximum = 10;
        double substitutionParameterMaximum = 100;
        double initialRootHeight = 1;
        double initialRate = 1;


        if (isFixedSubstitutionRate()) {
            double rate = getMeanSubstitutionRate();

            growthRateMaximum = 1E6 * rate;
            birthRateMaximum = 1E6 * rate;

            if (hasData()) {
                initialRootHeight = meanDistance / rate;

                initialRootHeight = round(initialRootHeight, 2);
            }

        } else {
            if (maximumTipHeight > 0) {
                initialRootHeight = maximumTipHeight * 10.0;
            }

            initialRate = round((meanDistance * 0.2) / initialRootHeight, 2);
        }

        double timeScaleMaximum = round(initialRootHeight * 1000.0, 2);

        for (Parameter param : parameters) {
            if (dataReset) param.priorEdited = false;

            if (!param.priorEdited) {
                switch (param.scale) {
                    case TIME_SCALE:
                        param.uniformLower = Math.max(0.0, param.lower);
                        param.uniformUpper = Math.min(timeScaleMaximum, param.upper);
                        param.initial = initialRootHeight;
                        break;
                    case T50_SCALE:
                        param.uniformLower = Math.max(0.0, param.lower);
                        param.uniformUpper = Math.min(timeScaleMaximum, param.upper);
                        param.initial = initialRootHeight / 5.0;
                        break;
                    case GROWTH_RATE_SCALE:
                        param.uniformLower = Math.max(-growthRateMaximum, param.lower);
                        param.uniformUpper = Math.min(growthRateMaximum, param.upper);
                        break;
                    case BIRTH_RATE_SCALE:
                        param.uniformLower = Math.max(0.0, param.lower);
                        param.uniformUpper = Math.min(birthRateMaximum, param.upper);
                        break;
                    case SUBSTITUTION_RATE_SCALE:
                        param.uniformLower = Math.max(0.0, param.lower);
                        param.uniformUpper = Math.min(substitutionRateMaximum, param.upper);
                        param.initial = initialRate;
                        break;
                    case LOG_STDEV_SCALE:
                        param.uniformLower = Math.max(0.0, param.lower);
                        param.uniformUpper = Math.min(logStdevMaximum, param.upper);
                        break;
                    case SUBSTITUTION_PARAMETER_SCALE:
                        param.uniformLower = Math.max(0.0, param.lower);
                        param.uniformUpper = Math.min(substitutionParameterMaximum, param.upper);
                        break;

                    case UNITY_SCALE:
                        param.uniformLower = 0.0;
                        param.uniformUpper = 1.0;
                        break;

                    case ROOT_RATE_SCALE:
                        param.initial = initialRate;
                        param.gammaAlpha = 0.5;
                        param.gammaBeta = param.initial / 0.5;
                        break;

                    case LOG_VAR_SCALE:
                        param.initial = initialRate;
                        param.gammaAlpha = 2.0;
                        param.gammaBeta = param.initial / 2.0;
                        break;

                }
                if (param.isNodeHeight) {
                    param.lower = maximumTipHeight;
                    param.uniformLower = maximumTipHeight;
                    param.uniformUpper = timeScaleMaximum;
                    param.initial = initialRootHeight;
                }
            }
        }

        dataReset = false;

        return parameters;
    }

    public boolean isFixedSubstitutionRate() {
        return fixedSubstitutionRate;
    }

    public double getMeanSubstitutionRate() {
        return meanSubstitutionRate;
    }

    public List<PartitionModel> getActiveModels() {

        Set<PartitionModel> models = new HashSet<PartitionModel>();

        for (DataPartition partition : dataPartitions) {
            models.add(partition.getPartitionModel());
        }

        // I have changed this to a list to ensure that the order is kept the same as in the table.
        List<PartitionModel> activeModels = new ArrayList<PartitionModel>();
        for (PartitionModel model : getPartitionModels()) {
            if (models.contains(model)) {
                activeModels.add(model);
            }
        }

        return activeModels;
    }

    public int getTotalActivePartitionCount() {
        int totalPartitionCount = 0;
        for (PartitionModel model : getActiveModels()) {
            totalPartitionCount += model.getCodonPartitionCount();
        }
        return totalPartitionCount;
    }

    /**
     * This returns an integer vector of the number of sites in each partition (including any codon partitions). These
     * are strictly in the same order as the 'mu' relative rates are listed.
     */
    public int[] getPartitionWeights() {
        int[] weights = new int[getTotalActivePartitionCount()];

        int k = 0;
        for (PartitionModel model : getActiveModels()) {
            for (DataPartition partition : dataPartitions) {
                if (partition.getPartitionModel() == model) {
                    model.addWeightsForPartition(partition, weights, k);
                }
            }
            k += model.getCodonPartitionCount();
        }

        assert (k == weights.length);

        return weights;
    }

    /**
     * return an list of operators that are required
     *
     * @return the operator list
     */
    public ArrayList<Operator> selectOperators() {

        ArrayList<Operator> ops = new ArrayList<Operator>();

        selectOperators(ops);

        boolean multiplePartitions = getTotalActivePartitionCount() > 1;

        for (PartitionModel model : getActiveModels()) {
            ops.addAll(model.getOperators());
        }

        if (multiplePartitions) {
            Operator deltaMuOperator = getOperator("deltaMu");
            ops.add(deltaMuOperator);
        }

        double initialRootHeight = 1;

        if (isFixedSubstitutionRate()) {
            double rate = getMeanSubstitutionRate();

            if (hasData()) {
                initialRootHeight = meanDistance / rate;
                initialRootHeight = round(initialRootHeight, 2);
            }

        } else {
            if (maximumTipHeight > 0) {
                initialRootHeight = maximumTipHeight * 10.0;
            }
        }

        Operator op = getOperator("subtreeSlide");
        if (!op.tuningEdited) {
            op.tuning = initialRootHeight / 10.0;
        }

        return ops;
    }

    /**
     * return a list of parameters that are required
     *
     * @param params the parameter list
     */
    private void selectParameters(ArrayList<Parameter> params) {

        if (hasData()) {

            // if not fixed then do mutation rate move and up/down move
            boolean fixed = isFixedSubstitutionRate();
            Parameter rateParam;

            switch (clockType) {
                case STRICT_CLOCK:
                    rateParam = getParameter("clock.rate");
                    if (!fixed) params.add(rateParam);
                    break;

                case UNCORRELATED_EXPONENTIAL:
                    rateParam = getParameter("uced.mean");
                    if (!fixed) params.add(rateParam);
                    break;

                case UNCORRELATED_LOGNORMAL:
                    rateParam = getParameter("ucld.mean");
                    if (!fixed) params.add(rateParam);
                    params.add(getParameter("ucld.stdev"));
                    break;

                case AUTOCORRELATED_LOGNORMAL:
                    rateParam = getParameter("treeModel.rootRate");
                    if (!fixed) params.add(rateParam);
                    params.add(getParameter("branchRates.var"));
                    break;

                case RANDOM_LOCAL_CLOCK:
                    rateParam = getParameter("clock.rate");
                    if (!fixed) params.add(rateParam);
                    break;

                default:
                    throw new IllegalArgumentException("Unknown clock model");
            }

            /*if (clockType == ClockType.STRICT_CLOCK || clockType == ClockType.RANDOM_LOCAL_CLOCK) {
				rateParam = getParameter("clock.rate");
				if (!fixed) params.add(rateParam);
			} else {
				if (clockType == ClockType.UNCORRELATED_EXPONENTIAL) {
					rateParam = getParameter("uced.mean");
					if (!fixed) params.add(rateParam);
				} else if (clockType == ClockType.UNCORRELATED_LOGNORMAL) {
					rateParam = getParameter("ucld.mean");
					if (!fixed) params.add(rateParam);
					params.add(getParameter("ucld.stdev"));
				} else {
					throw new IllegalArgumentException("Unknown clock model");
				}
			}*/

            rateParam.isFixed = fixed;

            if (errorModelType == ErrorType.AGE_ALL || errorModelType == ErrorType.AGE_TRANSITIONS) {
                params.add(getParameter("errorModel.ageRate"));
            } else if (errorModelType == ErrorType.BASE_ALL || errorModelType == ErrorType.BASE_TRANSITIONS) {
                params.add(getParameter("errorModel.baseRate"));
            }
        }

        if (nodeHeightPrior == TreePrior.CONSTANT) {
            params.add(getParameter("constant.popSize"));
        } else if (nodeHeightPrior == TreePrior.EXPONENTIAL) {
            params.add(getParameter("exponential.popSize"));
            if (parameterization == GROWTH_RATE) {
                params.add(getParameter("exponential.growthRate"));
            } else {
                params.add(getParameter("exponential.doublingTime"));
            }
        } else if (nodeHeightPrior == TreePrior.LOGISTIC) {
            params.add(getParameter("logistic.popSize"));
            if (parameterization == GROWTH_RATE) {
                params.add(getParameter("logistic.growthRate"));
            } else {
                params.add(getParameter("logistic.doublingTime"));
            }
            params.add(getParameter("logistic.t50"));
        } else if (nodeHeightPrior == TreePrior.EXPANSION) {
            params.add(getParameter("expansion.popSize"));
            if (parameterization == GROWTH_RATE) {
                params.add(getParameter("expansion.growthRate"));
            } else {
                params.add(getParameter("expansion.doublingTime"));
            }
            params.add(getParameter("expansion.ancestralProportion"));
        } else if (nodeHeightPrior == TreePrior.SKYLINE) {
            params.add(getParameter("skyline.popSize"));
        } else if (nodeHeightPrior == TreePrior.EXTENDED_SKYLINE) {
            params.add(getParameter("demographic.populationSizeChanges"));
            params.add(getParameter("demographic.populationMean"));
        } else if (nodeHeightPrior == TreePrior.GMRF_SKYRIDE) {
//            params.add(getParameter("skyride.popSize"));
            params.add(getParameter("skyride.precision"));
        } else if (nodeHeightPrior == TreePrior.YULE) {
            params.add(getParameter("yule.birthRate"));
        } else if (nodeHeightPrior == TreePrior.BIRTH_DEATH) {
            params.add(getParameter(BirthDeathModelParser.BIRTHDIFF_RATE_PARAM_NAME));
            params.add(getParameter(BirthDeathModelParser.RELATIVE_DEATH_RATE_PARAM_NAME));
        }

        params.add(getParameter("treeModel.rootHeight"));
    }

    private void selectStatistics(ArrayList<Parameter> params) {

        if (taxonSets != null) {
            for (Taxa taxonSet : taxonSets) {
                Parameter statistic = statistics.get(taxonSet);
                if (statistic == null) {
                    statistic = new Parameter(taxonSet, "tMRCA for taxon set ");
                    statistics.put(taxonSet, statistic);
                }
                System.err.println("Adding statistic " + statistic.getName());
                params.add(statistic);
            }
        } else {
            System.err.println("TaxonSets are null");
        }

        if (clockType == ClockType.RANDOM_LOCAL_CLOCK) {
            if (localClockRateChangesStatistic == null) {
                localClockRateChangesStatistic = new Parameter("rateChanges", "number of random local clocks", true);
                localClockRateChangesStatistic.priorType = PriorType.POISSON_PRIOR;
                localClockRateChangesStatistic.poissonMean = 1.0;
                localClockRateChangesStatistic.poissonOffset = 0.0;
            }
            if (localClockRatesStatistic == null) {
                localClockRatesStatistic = new Parameter("localClock.rates", "random local clock rates", false);

                localClockRatesStatistic.priorType = PriorType.GAMMA_PRIOR;
                localClockRatesStatistic.gammaAlpha = 0.5;
                localClockRatesStatistic.gammaBeta = 2.0;
            }
            params.add(localClockRatesStatistic);
            params.add(localClockRateChangesStatistic);
        }

        if (clockType != ClockType.STRICT_CLOCK) {
            params.add(getParameter("meanRate"));
            params.add(getParameter("coefficientOfVariation"));
            params.add(getParameter("covariance"));
        }
    }

    /**
     * return a list of operators that are required
     *
     * @param ops the operator list
     */
    private void selectOperators(ArrayList<Operator> ops) {

        if (hasData()) {

            if (!isFixedSubstitutionRate()) {
                switch (clockType) {
                    case STRICT_CLOCK:
                        ops.add(getOperator("clock.rate"));
                        ops.add(getOperator("upDownRateHeights"));
                        break;

                    case UNCORRELATED_EXPONENTIAL:
                        ops.add(getOperator("uced.mean"));
                        ops.add(getOperator("upDownUCEDMeanHeights"));
                        ops.add(getOperator("swapBranchRateCategories"));
                        ops.add(getOperator("randomWalkBranchRateCategories"));
                        ops.add(getOperator("unformBranchRateCategories"));
                        break;

                    case UNCORRELATED_LOGNORMAL:
                        ops.add(getOperator("ucld.mean"));
                        ops.add(getOperator("ucld.stdev"));
                        ops.add(getOperator("upDownUCLDMeanHeights"));
                        ops.add(getOperator("swapBranchRateCategories"));
                        ops.add(getOperator("randomWalkBranchRateCategories"));
                        ops.add(getOperator("unformBranchRateCategories"));
                        break;

                    case AUTOCORRELATED_LOGNORMAL:
                        ops.add(getOperator("scaleRootRate"));
                        ops.add(getOperator("scaleOneRate"));
                        ops.add(getOperator("scaleAllRates"));
                        ops.add(getOperator("scaleAllRatesIndependently"));
                        ops.add(getOperator("upDownAllRatesHeights"));
                        ops.add(getOperator("branchRates.var"));
                        break;

                    case RANDOM_LOCAL_CLOCK:
                        ops.add(getOperator("clock.rate"));
                        ops.add(getOperator("upDownRateHeights"));
                        ops.add(getOperator("localClock.rates"));
                        ops.add(getOperator("localClock.changes"));
                        ops.add(getOperator("treeBitMove"));
                        break;

                    default:
                        throw new IllegalArgumentException("Unknown clock model");
                }
            } else {
                switch (clockType) {
                    case STRICT_CLOCK:
                        // no parameter to operator on
                        break;

                    case UNCORRELATED_EXPONENTIAL:
                        ops.add(getOperator("swapBranchRateCategories"));
                        ops.add(getOperator("randomWalkBranchRateCategories"));
                        ops.add(getOperator("unformBranchRateCategories"));
                        break;

                    case UNCORRELATED_LOGNORMAL:
                        ops.add(getOperator("ucld.stdev"));
                        ops.add(getOperator("swapBranchRateCategories"));
                        ops.add(getOperator("randomWalkBranchRateCategories"));
                        ops.add(getOperator("unformBranchRateCategories"));
                        break;

                    case AUTOCORRELATED_LOGNORMAL:
                        ops.add(getOperator("scaleOneRate"));
                        ops.add(getOperator("scaleAllRatesIndependently"));
                        ops.add(getOperator("branchRates.var"));
                        break;

                    case RANDOM_LOCAL_CLOCK:
                        ops.add(getOperator("localClock.rates"));
                        ops.add(getOperator("localClock.changes"));
                        ops.add(getOperator("treeBitMove"));
                        break;

                    default:
                        throw new IllegalArgumentException("Unknown clock model");
                }
            }

            // if not fixed then do mutation rate move and up/down move
            /*if (!isFixedSubstitutionRate()) {
                   if (clockType == ClockType.STRICT_CLOCK) {
                       ops.add(getOperator("clock.rate"));
                       ops.add(getOperator("upDownRateHeights"));
                   } else if (clockType == ClockType.RANDOM_LOCAL_CLOCK) {
                       ops.add(getOperator("clock.rate"));
                       ops.add(getOperator("upDownRateHeights"));
                       ops.add(getOperator("localClock.rates"));
                       ops.add(getOperator("localClock.changes"));
                       ops.add(getOperator("treeBitMove"));
                   } else {
                       if (clockType == ClockType.UNCORRELATED_EXPONENTIAL) {
                           ops.add(getOperator("uced.mean"));
                           ops.add(getOperator("upDownUCEDMeanHeights"));
                       } else if (clockType == ClockType.UNCORRELATED_LOGNORMAL) {
                           ops.add(getOperator("ucld.mean"));
                           ops.add(getOperator("ucld.stdev"));
                           ops.add(getOperator("upDownUCLDMeanHeights"));
                       } else {
                           throw new IllegalArgumentException("Unknown clock model");
                       }
                       ops.add(getOperator("swapBranchRateCategories"));
                       ops.add(getOperator("randomWalkBranchRateCategories"));
                       ops.add(getOperator("unformBranchRateCategories"));
                   }
               } else {
                   if (clockType == ClockType.STRICT_CLOCK) {
                       // no parameter to operator on
                   } else if (clockType == ClockType.RANDOM_LOCAL_CLOCK) {
                       ops.add(getOperator("localClock.rates"));
                       ops.add(getOperator("localClock.changes"));
                       ops.add(getOperator("treeBitMove"));
                   } else {
                       if (clockType == ClockType.UNCORRELATED_EXPONENTIAL) {
                           // no parameter to operator on
                       } else if (clockType == ClockType.UNCORRELATED_LOGNORMAL) {
                           ops.add(getOperator("ucld.stdev"));
                       } else {
                           throw new IllegalArgumentException("Unknown clock model");
                       }
                       ops.add(getOperator("swapBranchRateCategories"));
                       ops.add(getOperator("randomWalkBranchRateCategories"));
                       ops.add(getOperator("unformBranchRateCategories"));
                   }
               }*/

            if (errorModelType == ErrorType.AGE_ALL || errorModelType == ErrorType.AGE_TRANSITIONS) {
                ops.add(getOperator("errorModel.ageRate"));
            } else if (errorModelType == ErrorType.BASE_ALL || errorModelType == ErrorType.BASE_TRANSITIONS) {
                ops.add(getOperator("errorModel.baseRate"));
            }
        }

        if (nodeHeightPrior == TreePrior.CONSTANT) {
            ops.add(getOperator("constant.popSize"));
        } else if (nodeHeightPrior == TreePrior.EXPONENTIAL) {
            ops.add(getOperator("exponential.popSize"));
            if (parameterization == GROWTH_RATE) {
                ops.add(getOperator("exponential.growthRate"));
            } else {
                ops.add(getOperator("exponential.doublingTime"));
            }
        } else if (nodeHeightPrior == TreePrior.LOGISTIC) {
            ops.add(getOperator("logistic.popSize"));
            if (parameterization == GROWTH_RATE) {
                ops.add(getOperator("logistic.growthRate"));
            } else {
                ops.add(getOperator("logistic.doublingTime"));
            }
            ops.add(getOperator("logistic.t50"));
        } else if (nodeHeightPrior == TreePrior.EXPANSION) {
            ops.add(getOperator("expansion.popSize"));
            if (parameterization == GROWTH_RATE) {
                ops.add(getOperator("expansion.growthRate"));
            } else {
                ops.add(getOperator("expansion.doublingTime"));
            }
            ops.add(getOperator("expansion.ancestralProportion"));
        } else if (nodeHeightPrior == TreePrior.SKYLINE) {
            ops.add(getOperator("skyline.popSize"));
            ops.add(getOperator("skyline.groupSize"));
        } else if (nodeHeightPrior == TreePrior.GMRF_SKYRIDE) {
            ops.add(getOperator("gmrfGibbsOperator"));
        } else if (nodeHeightPrior == TreePrior.EXTENDED_SKYLINE) {
            ops.add(getOperator("demographic.populationMean"));
            ops.add(getOperator("demographic.popSize"));
            ops.add(getOperator("demographic.indicators"));
            ops.add(getOperator("demographic.scaleActive"));
        } else if (nodeHeightPrior == TreePrior.YULE) {
            ops.add(getOperator("yule.birthRate"));
        } else if (nodeHeightPrior == TreePrior.BIRTH_DEATH) {
            ops.add(getOperator(BirthDeathModelParser.BIRTHDIFF_RATE_PARAM_NAME));
            ops.add(getOperator(BirthDeathModelParser.RELATIVE_DEATH_RATE_PARAM_NAME));
        }

        ops.add(getOperator("treeModel.rootHeight"));
        ops.add(getOperator("uniformHeights"));

        // if not a fixed tree then sample tree space
        if (!fixedTree) {
            ops.add(getOperator("subtreeSlide"));
            ops.add(getOperator("narrowExchange"));
            ops.add(getOperator("wideExchange"));
            ops.add(getOperator("wilsonBalding"));
        }
    }

    /**
     * Read options from a file
     *
     * @param guessDates guess dates?
     * @return the Document
     */
    public Document create(boolean guessDates) {

        Element root = new Element("beauti");
        root.setAttribute("version", version);

        Element dataElement = new Element("data");

        //dataElement.addContent(createChild("fileNameStem", fileNameStem));

        dataElement.addContent(createChild("datesUnits", datesUnits));
        dataElement.addContent(createChild("datesDirection", datesDirection));
        dataElement.addContent(createChild("translation", translation));
        dataElement.addContent(createChild("startingTreeType", startingTreeType.name()));

        dataElement.addContent(createChild("guessDates", guessDates));
        dataElement.addContent(createChild("guessDateFromOrder", dateGuesser.guessDateFromOrder));
        dataElement.addContent(createChild("fromLast", dateGuesser.fromLast));
        dataElement.addContent(createChild("order", dateGuesser.order));
        dataElement.addContent(createChild("prefix", dateGuesser.prefix));
        dataElement.addContent(createChild("offset", dateGuesser.offset));
        dataElement.addContent(createChild("unlessLessThan", dateGuesser.unlessLessThan));
        dataElement.addContent(createChild("offset2", dateGuesser.offset2));

        root.addContent(dataElement);

        Element taxaElement = new Element("taxa");

        for (Taxa taxonSet : taxonSets) {
            Element taxonSetElement = new Element("taxonSet");
            taxonSetElement.addContent(createChild("id", taxonSet.getId()));
            taxonSetElement.addContent(createChild("enforceMonophyly",
                    taxonSetsMono.get(taxonSet) ? "true" : "false"));
            for (int j = 0; j < taxonSet.getTaxonCount(); j++) {
                Element taxonElement = new Element("taxon");
                taxonElement.addContent(createChild("id", taxonSet.getTaxon(j).getId()));
                taxonSetElement.addContent(taxonElement);
            }
            taxaElement.addContent(taxonSetElement);
        }

        root.addContent(taxaElement);

        for (PartitionModel model : models) {

            Element modelElement = new Element("model");

            /*modelElement.addContent(createChild("nucSubstitutionModel", nucSubstitutionModel));
                           modelElement.addContent(createChild("aaSubstitutionModel", aaSubstitutionModel));
                           modelElement.addContent(createChild("binarySubstitutionModel", binarySubstitutionModel));
                           modelElement.addContent(createChild("frequencyPolicy", frequencyPolicy));
                           modelElement.addContent(createChild("gammaHetero", gammaHetero));
                           modelElement.addContent(createChild("gammaCategories", gammaCategories));
                           modelElement.addContent(createChild("invarHetero", invarHetero));
                           modelElement.addContent(createChild("codonHeteroPattern", codonHeteroPattern));
                           modelElement.addContent(createChild("maximumTipHeight", maximumTipHeight));
                           modelElement.addContent(createChild("hasSetFixedSubstitutionRate", hasSetFixedSubstitutionRate));
                           modelElement.addContent(createChild("meanSubstitutionRate", meanSubstitutionRate));
                           modelElement.addContent(createChild("fixedSubstitutionRate", fixedSubstitutionRate));
                           modelElement.addContent(createChild("unlinkedSubstitutionModel", unlinkedSubstitutionModel));
                           modelElement.addContent(createChild("unlinkedHeterogeneityModel", unlinkedHeterogeneityModel));
                           modelElement.addContent(createChild("unlinkedFrequencyModel", unlinkedFrequencyModel));
                           modelElement.addContent(createChild("clockModel", clockModel));
                           modelElement.addContent(createChild("nodeHeightPrior", nodeHeightPrior));
                           modelElement.addContent(createChild("parameterization", parameterization));
                           modelElement.addContent(createChild("skylineGroupCount", skylineGroupCount));
                           modelElement.addContent(createChild("skylineModel", skylineModel));
                           modelElement.addContent(createChild("fixedTree", fixedTree)); */

            root.addContent(modelElement);
        }

        Element priorsElement = new Element("priors");

        for (String name : parameters.keySet()) {
            Parameter parameter = parameters.get(name);
            Element e = new Element(name);
            e.addContent(createChild("initial", parameter.initial));
            e.addContent(createChild("priorType", parameter.priorType));
            e.addContent(createChild("priorEdited", parameter.priorEdited));
            e.addContent(createChild("uniformLower", parameter.uniformLower));
            e.addContent(createChild("uniformUpper", parameter.uniformUpper));
            e.addContent(createChild("exponentialMean", parameter.exponentialMean));
            e.addContent(createChild("exponentialOffset", parameter.exponentialOffset));
            e.addContent(createChild("normalMean", parameter.normalMean));
            e.addContent(createChild("normalStdev", parameter.normalStdev));
            e.addContent(createChild("logNormalMean", parameter.logNormalMean));
            e.addContent(createChild("logNormalStdev", parameter.logNormalStdev));
            e.addContent(createChild("logNormalOffset", parameter.logNormalOffset));
            e.addContent(createChild("gammaAlpha", parameter.gammaAlpha));
            e.addContent(createChild("gammaBeta", parameter.gammaBeta));
            e.addContent(createChild("gammaOffset", parameter.gammaOffset));
            priorsElement.addContent(e);
        }

        for (Taxa taxonSet : taxonSets) {
            Parameter statistic = statistics.get(taxonSet);
            Element e = new Element(statistic.getXMLName());
            e.addContent(createChild("initial", statistic.initial));
            e.addContent(createChild("priorType", statistic.priorType));
            e.addContent(createChild("priorEdited", statistic.priorEdited));
            e.addContent(createChild("uniformLower", statistic.uniformLower));
            e.addContent(createChild("uniformUpper", statistic.uniformUpper));
            e.addContent(createChild("exponentialMean", statistic.exponentialMean));
            e.addContent(createChild("exponentialOffset", statistic.exponentialOffset));
            e.addContent(createChild("normalMean", statistic.normalMean));
            e.addContent(createChild("normalStdev", statistic.normalStdev));
            e.addContent(createChild("logNormalMean", statistic.logNormalMean));
            e.addContent(createChild("logNormalStdev", statistic.logNormalStdev));
            e.addContent(createChild("logNormalOffset", statistic.logNormalOffset));
            e.addContent(createChild("gammaAlpha", statistic.gammaAlpha));
            e.addContent(createChild("gammaBeta", statistic.gammaBeta));
            e.addContent(createChild("gammaOffset", statistic.gammaOffset));
            priorsElement.addContent(e);
        }

        root.addContent(priorsElement);

        Element operatorsElement = new Element("operators");

        operatorsElement.addContent(createChild("autoOptimize", autoOptimize));
        for (String name : operators.keySet()) {
            Operator operator = operators.get(name);
            Element e = new Element(name);
            e.addContent(createChild("tuning", operator.tuning));
            e.addContent(createChild("tuningEdited", operator.tuningEdited));
            e.addContent(createChild("weight", operator.weight));
            e.addContent(createChild("inUse", operator.inUse));
            operatorsElement.addContent(e);
        }

        root.addContent(operatorsElement);

        Element mcmcElement = new Element("mcmc");

        mcmcElement.addContent(createChild("chainLength", chainLength));
        mcmcElement.addContent(createChild("logEvery", logEvery));
        mcmcElement.addContent(createChild("echoEvery", echoEvery));
        //if (logFileName != null) mcmcElement.addContent(createChild("logFileName", logFileName));
        //if (treeFileName != null) mcmcElement.addContent(createChild("treeFileName", treeFileName));
        //mcmcElement.addContent(createChild("mapTreeLog", mapTreeLog));
        //if (mapTreeFileName != null) mcmcElement.addContent(createChild("mapTreeFileName", mapTreeFileName));
        mcmcElement.addContent(createChild("substTreeLog", substTreeLog));
        //if (substTreeFileName != null) mcmcElement.addContent(createChild("substTreeFileName", substTreeFileName));

        root.addContent(mcmcElement);

        return new Document(root);
    }

    private Element createChild(String name, String value) {
        Element e = new Element(name);
        if (value != null) {
            e.setText(value);
        }
        return e;
    }

    private Element createChild(String name, int value) {
        Element e = new Element(name);
        e.setText(Integer.toString(value));
        return e;
    }

    private Element createChild(String name, PriorType value) {
        Element e = new Element(name);
        e.setText(value.name());
        return e;
    }

    private Element createChild(String name, double value) {
        Element e = new Element(name);
        e.setText(Double.toString(value));
        return e;
    }

    private Element createChild(String name, boolean value) {
        Element e = new Element(name);
        e.setText(value ? "true" : "false");
        return e;
    }

    /**
     * Read options from a file
     *
     * @param document the Document
     * @throws dr.xml.XMLParseException if there is a problem with XML parsing
     */
    public void parse(Document document) throws dr.xml.XMLParseException {

        Element root = document.getRootElement();
        if (!root.getName().equals("beauti")) {
            throw new dr.xml.XMLParseException("This document does not appear to be a BEAUti file");
        }

        Element taxaElement = root.getChild("taxa");
        Element modelElement = root.getChild("model");
        Element priorsElement = root.getChild("priors");
        Element operatorsElement = root.getChild("operators");
        Element mcmcElement = root.getChild("mcmc");
        /*
                  if (taxaElement != null) {
                      for (Object ts : taxaElement.getChildren("taxonSet")) {
                          Element taxonSetElement = (Element) ts;

                          String id = getStringChild(taxonSetElement, "id", "");
                          final Taxa taxonSet = new Taxa(id);

                          Boolean enforceMonophyly = Boolean.valueOf(getStringChild(taxonSetElement, "enforceMonophyly", "false"));
                          for (Object o : taxonSetElement.getChildren("taxon")) {
                              Element taxonElement = (Element) o;
                              String taxonId = getStringChild(taxonElement, "id", "");
                              int index = taxonList.getTaxonIndex(taxonId);
                              if (index != -1) {
                                  taxonSet.addTaxon(taxonList.getTaxon(index));
                              }
                          }
                          taxonSets.add(taxonSet);
                          taxonSetsMono.put(taxonSet, enforceMonophyly);
                      }
                  }

                  if (modelElement != null) {
                      nucSubstitutionModel = getIntegerChild(modelElement, "nucSubstitutionModel", HKY);
                      aaSubstitutionModel = getIntegerChild(modelElement, "aaSubstitutionModel", BLOSUM_62);
                      binarySubstitutionModel = getIntegerChild(modelElement, "binarySubstitutionModel", BIN_SIMPLE);
                      frequencyPolicy = getIntegerChild(modelElement, "frequencyPolicy", ESTIMATED);
                      gammaHetero = getBooleanChild(modelElement, "gammaHetero", false);
                      gammaCategories = getIntegerChild(modelElement, "gammaCategories", 5);
                      invarHetero = getBooleanChild(modelElement, "invarHetero", false);
                      codonHeteroPattern = (getBooleanChild(modelElement, "codonHetero", false) ? "123" : null);
                      codonHeteroPattern = getStringChild(modelElement, "codonHeteroPattern", null);
                      maximumTipHeight = getDoubleChild(modelElement, "maximumTipHeight", 0.0);
                      fixedSubstitutionRate = getBooleanChild(modelElement, "fixedSubstitutionRate", false);
                      hasSetFixedSubstitutionRate = getBooleanChild(modelElement, "hasSetFixedSubstitutionRate", false);
                      meanSubstitutionRate = getDoubleChild(modelElement, "meanSubstitutionRate", 1.0);
                      unlinkedSubstitutionModel = getBooleanChild(modelElement, "unlinkedSubstitutionModel", false);
                      unlinkedHeterogeneityModel = getBooleanChild(modelElement, "unlinkedHeterogeneityModel", false);
                      unlinkedFrequencyModel = getBooleanChild(modelElement, "unlinkedFrequencyModel", false);

                      clockModel = getIntegerChild(modelElement, "clockModel", clockModel);

                      // the old name was "coalescentModel" so try to read this first
                      nodeHeightPrior = getIntegerChild(modelElement, "coalescentModel", CONSTANT);
                      nodeHeightPrior = getIntegerChild(modelElement, "nodeHeightPrior", nodeHeightPrior);
                      // we don't allow no nodeHeightPrior in BEAUti so switch it to Yule:
                      if (nodeHeightPrior == NONE) nodeHeightPrior = YULE;

                      parameterization = getIntegerChild(modelElement, "parameterization", GROWTH_RATE);
                      skylineGroupCount = getIntegerChild(modelElement, "skylineGroupCount", 10);
                      skylineModel = getIntegerChild(modelElement, "skylineModel", CONSTANT_SKYLINE);
                      fixedTree = getBooleanChild(modelElement, "fixedTree", false);
                  }

                  if (operatorsElement != null) {
                      autoOptimize = getBooleanChild(operatorsElement, "autoOptimize", true);
                      for (String name : operators.keySet()) {
                          Operator operator = operators.get(name);
                          Element e = operatorsElement.getChild(name);
                          if (e == null) {
                              throw new XMLParseException("Operators element, " + name + " missing");
                          }

                          operator.tuning = getDoubleChild(e, "tuning", 1.0);
                          operator.tuningEdited = getBooleanChild(e, "tuningEdited", false);
                          operator.weight = getDoubleChild(e, "weight", 1);
                          operator.inUse = getBooleanChild(e, "inUse", true);
                      }
                  }

                  if (priorsElement != null) {
                      for (String name : parameters.keySet()) {
                          Parameter parameter = parameters.get(name);
                          Element e = priorsElement.getChild(name);
                          if (e == null) {
                              throw new XMLParseException("Priors element, " + name + " missing");
                          }

                          parameter.initial = getDoubleChild(e, "initial", 1.0);
                          parameter.priorType = PriorType.valueOf(getStringChild(e, "priorType", PriorType.UNIFORM_PRIOR.name()));
                          parameter.priorEdited = getBooleanChild(e, "priorEdited", false);
                          parameter.uniformLower = Math.max(getDoubleChild(e, "uniformLower", parameter.uniformLower), parameter.lower);
                          parameter.uniformUpper = Math.min(getDoubleChild(e, "uniformUpper", parameter.uniformUpper), parameter.upper);
                          parameter.exponentialMean = getDoubleChild(e, "exponentialMean", parameter.exponentialMean);
                          parameter.exponentialOffset = getDoubleChild(e, "exponentialOffset", parameter.exponentialOffset);
                          parameter.normalMean = getDoubleChild(e, "normalMean", parameter.normalMean);
                          parameter.normalStdev = getDoubleChild(e, "normalStdev", parameter.normalStdev);
                          parameter.logNormalMean = getDoubleChild(e, "logNormalMean", parameter.logNormalMean);
                          parameter.logNormalStdev = getDoubleChild(e, "logNormalStdev", parameter.logNormalStdev);
                          parameter.logNormalOffset = getDoubleChild(e, "logNormalOffset", parameter.logNormalOffset);
                          parameter.gammaAlpha = getDoubleChild(e, "gammaAlpha", parameter.gammaAlpha);
                          parameter.gammaBeta = getDoubleChild(e, "gammaBeta", parameter.gammaBeta);
                          parameter.gammaOffset = getDoubleChild(e, "gammaOffset", parameter.gammaOffset);
                      }

                      for (Taxa taxonSet : taxonSets) {
                          Parameter statistic = statistics.get(taxonSet);
                          if (statistic == null) {
                              statistic = new Parameter(this, taxonSet, "tMRCA for taxon set ");
                              statistics.put(taxonSet, statistic);
                          }
                          Element e = priorsElement.getChild(statistic.getXMLName());
                          statistic.initial = getDoubleChild(e, "initial", 1.0);
                          statistic.priorType = PriorType.valueOf(getStringChild(e, "priorType", PriorType.UNIFORM_PRIOR.name()));
                          statistic.priorEdited = getBooleanChild(e, "priorEdited", false);
                          statistic.uniformLower = getDoubleChild(e, "uniformLower", statistic.uniformLower);
                          statistic.uniformUpper = getDoubleChild(e, "uniformUpper", statistic.uniformUpper);
                          statistic.exponentialMean = getDoubleChild(e, "exponentialMean", statistic.exponentialMean);
                          statistic.exponentialOffset = getDoubleChild(e, "exponentialOffset", statistic.exponentialOffset);
                          statistic.normalMean = getDoubleChild(e, "normalMean", statistic.normalMean);
                          statistic.normalStdev = getDoubleChild(e, "normalStdev", statistic.normalStdev);
                          statistic.logNormalMean = getDoubleChild(e, "logNormalMean", statistic.logNormalMean);
                          statistic.logNormalStdev = getDoubleChild(e, "logNormalStdev", statistic.logNormalStdev);
                          statistic.logNormalOffset = getDoubleChild(e, "logNormalOffset", statistic.logNormalOffset);
                          statistic.gammaAlpha = getDoubleChild(e, "gammaAlpha", statistic.gammaAlpha);
                          statistic.gammaBeta = getDoubleChild(e, "gammaBeta", statistic.gammaBeta);
                          statistic.gammaOffset = getDoubleChild(e, "gammaOffset", statistic.gammaOffset);
                      }

                  }


                  if (mcmcElement != null) {
                      upgmaStartingTree = getBooleanChild(mcmcElement, "upgmaStartingTree", true);
                      chainLength = getIntegerChild(mcmcElement, "chainLength", 100000000);
                      logEvery = getIntegerChild(mcmcElement, "logEvery", 1000);
                      echoEvery = getIntegerChild(mcmcElement, "echoEvery", 1000);
                      logFileName = getStringChild(mcmcElement, "logFileName", null);
                      treeFileName = getStringChild(mcmcElement, "treeFileName", null);
                      mapTreeLog = getBooleanChild(mcmcElement, "mapTreeLog", false);
                      mapTreeFileName = getStringChild(mcmcElement, "mapTreeFileName", null);
                      substTreeLog = getBooleanChild(mcmcElement, "substTreeLog", false);
                      substTreeFileName = getStringChild(mcmcElement, "substTreeFileName", null);
                  }      */
    }

    private String getStringChild(Element element, String childName, String defaultValue) {
        String value = element.getChildTextTrim(childName);
        if (value == null || value.length() == 0) return defaultValue;
        return value;
    }

    private int getIntegerChild(Element element, String childName, int defaultValue) {
        String value = element.getChildTextTrim(childName);
        if (value == null) return defaultValue;
        return Integer.parseInt(value);
    }

    private double getDoubleChild(Element element, String childName, double defaultValue) {
        String value = element.getChildTextTrim(childName);
        if (value == null) return defaultValue;
        return Double.parseDouble(value);
    }

    private boolean getBooleanChild(Element element, String childName, boolean defaultValue) {
        String value = element.getChildTextTrim(childName);
        if (value == null) return defaultValue;
        return value.equals("true");
    }

    private Date createDate(double timeValue, Units.Type units, boolean backwards, double origin) {
        if (backwards) {
            return Date.createTimeAgoFromOrigin(timeValue, units, origin);
        } else {
            return Date.createTimeSinceOrigin(timeValue, units, origin);
        }
    }

    public boolean hasData() {
        return dataPartitions.size() > 0;
    }

    public String fileNameStem = "untitled";
    public String logFileName = null;
    public String treeFileName = null;
    public boolean mapTreeLog = false;
    public String mapTreeFileName = null;
    public boolean substTreeLog = false;
    public String substTreeFileName = null;

    // Data options
    public DataType dataType = null;
    public boolean dataReset = true;

    public TaxonList taxonList = null;
    public DateGuesser dateGuesser = new DateGuesser();

    public List<Taxa> taxonSets = new ArrayList<Taxa>();
    public Map<Taxa, Boolean> taxonSetsMono = new HashMap<Taxa, Boolean>();
    public List<DataPartition> dataPartitions = new ArrayList<DataPartition>();
    public List<Tree> trees = new ArrayList<Tree>();
    public double meanDistance = 1.0;
    public int datesUnits = YEARS;
    public int datesDirection = FORWARDS;
    public double maximumTipHeight = 0.0;
    public int translation = 0;
    public StartingTreeType startingTreeType = StartingTreeType.RANDOM;
    public Tree userStartingTree = null;

    // Model options
    List<PartitionModel> models = new ArrayList<PartitionModel>();

    public boolean fixedSubstitutionRate = true;
    public double meanSubstitutionRate = 1.0;
    public boolean unlinkPartitionRates = true;

    public TreePrior nodeHeightPrior = TreePrior.CONSTANT;
    public int parameterization = GROWTH_RATE;
    public int skylineGroupCount = 10;
    public int skylineModel = CONSTANT_SKYLINE;
    public int skyrideSmoothing = SKYRIDE_TIME_AWARE_SMOOTHING;
    public String extendedSkylineModel = VariableDemographicModel.LINEAR;
    public boolean multiLoci = false;
    public double birthDeathSamplingProportion = 1.0;
    public boolean fixedTree = false;

    public Units.Type units = Units.Type.SUBSTITUTIONS;
    public ClockType clockType = ClockType.STRICT_CLOCK;
    public ErrorType errorModelType = ErrorType.NO_ERROR;

    // Operator schedule options
    public int coolingSchedule = OperatorSchedule.DEFAULT_SCHEDULE;

    // MCMC options
    public int chainLength = 10000000;
    public int logEvery = 1000;
    public int echoEvery = 1000;
    public int burnIn = 100000;
    public String fileName = null;
    public boolean autoOptimize = true;
    public boolean performTraceAnalysis = false;
    public boolean generateCSV = true;  // until/if a button
    public boolean samplePriorOnly = false;

    public Parameter localClockRateChangesStatistic = null;
    public Parameter localClockRatesStatistic = null;

}