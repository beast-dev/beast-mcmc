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

import dr.app.beauti.priorsPanel.PriorType;
import dr.evomodel.tree.RateStatistic;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @author Walter Xie
 * @version $Id$
 */
public class PartitionClockModel extends ModelOptions {

    // Instance variables

    private final BeautiOptions options;
    private String name;

    private List<PartitionData> allPartitionData;

    private ClockType clockType = ClockType.STRICT_CLOCK;

    public PartitionClockModel(BeautiOptions options, PartitionData partition) {
        this.options = options;
        this.name = partition.getName();

        allPartitionData = new ArrayList<PartitionData>();
        addPartitionData(partition);

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
        this.name = name;

        this.allPartitionData = source.allPartitionData;

        this.clockType = source.clockType;

        initClockModelParaAndOpers();
    }

//    public PartitionClockModel(BeautiOptions options, String name) {
//        this.options = options;
//        this.name = name;
//    }

    public void initClockModelParaAndOpers() {
        double rateWeights = 3.0;

        createParameter("clock.rate", "substitution rate", SUBSTITUTION_RATE_SCALE, 1.0, 0.0, Double.POSITIVE_INFINITY);
        createParameter(ClockType.UCED_MEAN, "uncorrelated exponential relaxed clock mean", SUBSTITUTION_RATE_SCALE, 1.0, 0.0, Double.POSITIVE_INFINITY);
        createParameter(ClockType.UCLD_MEAN, "uncorrelated lognormal relaxed clock mean", SUBSTITUTION_RATE_SCALE, 1.0, 0.0, Double.POSITIVE_INFINITY);
        createParameter(ClockType.UCLD_STDEV, "uncorrelated lognormal relaxed clock stdev", LOG_STDEV_SCALE, 0.1, 0.0, Double.POSITIVE_INFINITY);

        {
            final Parameter p = createParameter("branchRates.var", "autocorrelated lognormal relaxed clock rate variance ", LOG_VAR_SCALE, 0.1, 0.0, Double.POSITIVE_INFINITY);
            p.priorType = PriorType.GAMMA_PRIOR;
            p.gammaAlpha = 1;
            p.gammaBeta = 0.0001;
        }

        createScaleOperator("clock.rate", demoTuning, rateWeights);
        createScaleOperator(ClockType.UCED_MEAN, demoTuning, rateWeights);
        createScaleOperator(ClockType.UCLD_MEAN, demoTuning, rateWeights);
        createScaleOperator(ClockType.UCLD_STDEV, demoTuning, rateWeights);

        createScaleOperator("branchRates.var", demoTuning, rateWeights);

//        for (PartitionTreeModel tree : options.getPartitionTreeModels(getAllPartitionData())) { // borrow the method in BeautiOption        	
//            createOperator(tree.getPrefix() + "upDownRateHeights", "Substitution rate and heights",
//                    "Scales substitution rates inversely to node heights of the tree", super.getParameter("clock.rate"),
//                    tree.getParameter("treeModel.allInternalNodeHeights"), OperatorType.UP_DOWN, 0.75, rateWeights);
//            createOperator(tree.getPrefix() + "upDownUCEDMeanHeights", "UCED mean and heights",
//                    "Scales UCED mean inversely to node heights of the tree", super.getParameter(ClockType.UCED_MEAN),
//                    tree.getParameter("treeModel.allInternalNodeHeights"), OperatorType.UP_DOWN, 0.75, rateWeights);
//            createOperator(tree.getPrefix() + "upDownUCLDMeanHeights", "UCLD mean and heights",
//                    "Scales UCLD mean inversely to node heights of the tree", super.getParameter(ClockType.UCLD_MEAN),
//                    tree.getParameter("treeModel.allInternalNodeHeights"), OperatorType.UP_DOWN, 0.75, rateWeights);
//            
//            // #meanRate = #Relaxed Clock Model * #Tree Model
//            createStatistic(tree.getPrefix() + "meanRate", "The mean rate of evolution over the whole tree", 0.0, Double.POSITIVE_INFINITY);
//            // #covariance = #Relaxed Clock Model * #Tree Model
//            createStatistic(tree.getPrefix() + "covariance", "The covariance in rates of evolution on each lineage with their ancestral lineages",
//                    Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY);
//            
//        }  move to PartitionClockModelTreeModelLink

        // These are statistics which could have priors on...        
        // #COEFFICIENT_OF_VARIATION = #Uncorrelated Clock Model
        createStatistic(RateStatistic.COEFFICIENT_OF_VARIATION, "The variation in rate of evolution over the whole tree",
                0.0, Double.POSITIVE_INFINITY);

    }

    /**
     * return a list of parameters that are required
     *
     * @param params the parameter list
     */
    public void selectParameters(List<Parameter> params) {
        if (options.hasData()) {

            // if not fixed then do mutation rate move and up/down move
            boolean fixed = options.isFixedSubstitutionRate();
            Parameter rateParam;

            switch (clockType) {
                case STRICT_CLOCK:
                    rateParam = getParameter("clock.rate");
                    rateParam.isFixed = fixed;
                    if (!fixed) params.add(rateParam);
                    break;

                case UNCORRELATED_EXPONENTIAL:
                    rateParam = getParameter(ClockType.UCED_MEAN);
                    rateParam.isFixed = fixed;
                    if (!fixed) params.add(rateParam);
                    break;

                case UNCORRELATED_LOGNORMAL:
                    rateParam = getParameter(ClockType.UCLD_MEAN);
                    rateParam.isFixed = fixed;
                    if (!fixed) params.add(rateParam);
                    params.add(getParameter(ClockType.UCLD_STDEV));
                    break;

                case AUTOCORRELATED_LOGNORMAL:
                    for (PartitionTreeModel tree : options.getPartitionTreeModels(getAllPartitionData())) {
                        rateParam = tree.getParameter("treeModel.rootRate");
                        rateParam.isFixed = fixed;
                        if (!fixed) params.add(rateParam);
                    }
                    params.add(getParameter("branchRates.var"));
                    break;

                case RANDOM_LOCAL_CLOCK:
                    rateParam = getParameter("clock.rate");
                    rateParam.isFixed = fixed;
                    if (!fixed) params.add(rateParam);
                    break;

                default:
                    throw new IllegalArgumentException("Unknown clock model");
            }

            selectStatistics(params);
        }
    }

    private void selectStatistics(List<Parameter> params) {
        // Statistics
        if (clockType != ClockType.STRICT_CLOCK) {
            params.add(getParameter(RateStatistic.COEFFICIENT_OF_VARIATION));
        }
    }

    /**
     * return a list of operators that are required
     *
     * @param ops the operator list
     */
    public void selectOperators(List<Operator> ops) {
        if (options.hasData()) {

            if (!options.isFixedSubstitutionRate()) {
                switch (clockType) {
                    case STRICT_CLOCK:
                        ops.add(getOperator("clock.rate"));
//                        ops.add(getOperator("upDownRateHeights"));
                        break;

                    case UNCORRELATED_EXPONENTIAL:
                        ops.add(getOperator(ClockType.UCED_MEAN));
//                        ops.add(getOperator("upDownUCEDMeanHeights"));
                        addBranchRateCategories(ops);
                        break;

                    case UNCORRELATED_LOGNORMAL:
                        ops.add(getOperator(ClockType.UCLD_MEAN));
                        ops.add(getOperator(ClockType.UCLD_STDEV));
//                        ops.add(getOperator("upDownUCLDMeanHeights"));
                        addBranchRateCategories(ops);
                        break;

                    case AUTOCORRELATED_LOGNORMAL:
                        for (PartitionTreeModel tree : options.getPartitionTreeModels(getAllPartitionData())) {
                            ops.add(tree.getOperator("scaleRootRate"));
                            ops.add(tree.getOperator("scaleOneRate"));
                            ops.add(tree.getOperator("scaleAllRates"));
                            ops.add(tree.getOperator("scaleAllRatesIndependently"));
                            ops.add(tree.getOperator("upDownAllRatesHeights"));
//                            ops.add(tree.getOperator("upDownNodeRatesHeights"));
                        }
                        ops.add(getOperator("branchRates.var"));
                        break;

                    case RANDOM_LOCAL_CLOCK:
                        ops.add(getOperator("clock.rate"));
//                        ops.add(getOperator("upDownRateHeights"));
                        addRandomLocalClockOperators(ops);
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
                        addBranchRateCategories(ops);
                        break;

                    case UNCORRELATED_LOGNORMAL:
                        addBranchRateCategories(ops);
                        ops.add(getOperator(ClockType.UCLD_STDEV));
                        break;

                    case AUTOCORRELATED_LOGNORMAL:
                        for (PartitionTreeModel tree : options.getPartitionTreeModels(getAllPartitionData())) {
                            ops.add(tree.getOperator("scaleOneRate"));
                            ops.add(tree.getOperator("scaleAllRatesIndependently"));
                        }
                        ops.add(getOperator("branchRates.var"));
                        break;

                    case RANDOM_LOCAL_CLOCK:
                        addRandomLocalClockOperators(ops);
                        break;

                    default:
                        throw new IllegalArgumentException("Unknown clock model");
                }
            }
        }
    }

    private void addBranchRateCategories(List<Operator> ops) {
        for (PartitionTreeModel tree : options.getPartitionTreeModels(getAllPartitionData())) {
            ops.add(tree.getOperator("swapBranchRateCategories"));
            ops.add(tree.getOperator("randomWalkBranchRateCategories"));
            ops.add(tree.getOperator("uniformBranchRateCategories"));
        }
    }

    private void addRandomLocalClockOperators(List<Operator> ops) {
        for (PartitionTreeModel tree : options.getPartitionTreeModels(getAllPartitionData())) {
            ops.add(tree.getOperator(ClockType.LOCAL_CLOCK + ".rates"));
            ops.add(tree.getOperator(ClockType.LOCAL_CLOCK + ".changes"));
            ops.add(tree.getOperator("treeBitMove"));
        }
    }

    /////////////////////////////////////////////////////////////

    public List<PartitionData> getAllPartitionData() {
        return allPartitionData;
    }

    public void clearAllPartitionData() {
        this.allPartitionData.clear();
    }

    public void addPartitionData(PartitionData partition) {
        allPartitionData.add(partition);
    }

    public boolean removePartitionData(PartitionData partition) {
        return allPartitionData.remove(partition);
    }

    public void setClockType(ClockType clockType) {
        this.clockType = clockType;
    }

    public ClockType getClockType() {
        return clockType;
    }


    public String getName() {
        return name;
    }

    public void setName(String name) {
        this.name = name;
    }

    public String toString() {
        return getName();
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

    public String getPrefix() {
        String prefix = "";
        if (options.getPartitionTreeModels().size() > 1) { //|| options.isSpeciesAnalysis()
            // There is more than one active partition model
            prefix += getName() + ".";
        }
        return prefix;
    }
}
