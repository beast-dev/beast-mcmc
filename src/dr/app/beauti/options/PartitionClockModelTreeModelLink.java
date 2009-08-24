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

import java.util.List;

import dr.app.beauti.enumTypes.ClockType;
import dr.app.beauti.enumTypes.FixRateType;
import dr.app.beauti.enumTypes.OperatorType;
import dr.app.beauti.enumTypes.PriorScaleType;
import dr.app.beauti.enumTypes.PriorType;
import dr.evomodel.tree.RateStatistic;

/**
 * @author Andrew Rambaut
 * @author Walter Xie
 * @version $Id$
 */
public class PartitionClockModelTreeModelLink extends ModelOptions {

    // Instance variables
    public Parameter localClockRateChangesStatistic = null;
    public Parameter localClockRatesStatistic = null;

    private final BeautiOptions options;
    private final PartitionClockModel model;    
	private final PartitionTreeModel tree;

    public PartitionClockModelTreeModelLink(BeautiOptions options, PartitionClockModel model, PartitionTreeModel tree) {
        this.options = options;
        this.model = model;
        this.tree = tree;

        initClockModelParaAndOpers();
    }

    private void initClockModelParaAndOpers() {
        {
            final Parameter p = createParameter("branchRates.var", "autocorrelated lognormal relaxed clock rate variance ", PriorScaleType.LOG_VAR_SCALE, 0.1, 0.0, Double.POSITIVE_INFINITY);
            p.priorType = PriorType.GAMMA_PRIOR;
            p.gammaAlpha = 1;
            p.gammaBeta = 0.0001;
        }

        createParameter("branchRates.categories", "relaxed clock branch rate categories");
        createParameter(ClockType.LOCAL_CLOCK + "." + "rates", "random local clock rates", PriorScaleType.SUBSTITUTION_RATE_SCALE, 1.0, 0.0, Double.POSITIVE_INFINITY);
        createParameter(ClockType.LOCAL_CLOCK + "." + "changes", "random local clock rate change indicator");

        {
            final Parameter p = createParameter("treeModel.rootRate", "autocorrelated lognormal relaxed clock root rate", PriorScaleType.ROOT_RATE_SCALE, 1.0, 0.0, Double.POSITIVE_INFINITY);
            p.priorType = PriorType.GAMMA_PRIOR;
            p.gammaAlpha = 1;
            p.gammaBeta = 0.0001;
        }
        createParameter("treeModel.nodeRates", "autocorrelated lognormal relaxed clock non-root rates", PriorScaleType.SUBSTITUTION_RATE_SCALE, 1.0, 0.0, Double.POSITIVE_INFINITY);
        createParameter("treeModel.allRates", "autocorrelated lognormal relaxed clock all rates", PriorScaleType.SUBSTITUTION_RATE_SCALE, 1.0, 0.0, Double.POSITIVE_INFINITY);

        
        createScaleOperator("branchRates.var", demoTuning, rateWeights);

        createOperator("scaleRootRate", "treeModel.rootRate", "Scales root rate", "treeModel.rootRate", OperatorType.SCALE, demoTuning, rateWeights);
        createOperator("scaleOneRate", "treeModel.nodeRates", "Scales one non-root rate", "treeModel.nodeRates", 
        		OperatorType.SCALE, demoTuning, branchWeights);
        createOperator("scaleAllRates", "treeModel.allRates", "Scales all rates simultaneously", "treeModel.allRates",
                OperatorType.SCALE_ALL, demoTuning, rateWeights);
        createOperator("scaleAllRatesIndependently", "treeModel.nodeRates", "Scales all non-root rates independently", "treeModel.nodeRates",
                OperatorType.SCALE_INDEPENDENTLY, demoTuning, rateWeights);
        
        createOperator("swapBranchRateCategories", "branchRates.categories", "Performs a swap of branch rate categories", 
        		"branchRates.categories", OperatorType.SWAP, 1, branchWeights / 3);
        createOperator("randomWalkBranchRateCategories", "branchRates.categories", "Performs an integer random walk of branch rate categories", 
        		"branchRates.categories", OperatorType.INTEGER_RANDOM_WALK, 1, branchWeights / 3);
        createOperator("uniformBranchRateCategories", "branchRates.categories", "Performs an integer uniform draw of branch rate categories", 
        		"branchRates.categories", OperatorType.INTEGER_UNIFORM, 1, branchWeights / 3);

        createScaleOperator(ClockType.LOCAL_CLOCK + "." + "rates", demoTuning, treeWeights);
        createOperator(ClockType.LOCAL_CLOCK + "." + "changes", OperatorType.BITFLIP, 1, treeWeights);
         
        createOperator("upDownRateHeights", "Substitution rate and heights",
                "Scales substitution rates inversely to node heights of the tree", model.getParameter("clock.rate"),
                tree.getParameter("treeModel.allInternalNodeHeights"), OperatorType.UP_DOWN, demoTuning, rateWeights);
        createOperator("upDownUCEDMeanHeights", "UCED mean and heights",
                "Scales UCED mean inversely to node heights of the tree", model.getParameter(ClockType.UCED_MEAN),
                tree.getParameter("treeModel.allInternalNodeHeights"), OperatorType.UP_DOWN, demoTuning, rateWeights);
        createOperator("upDownUCLDMeanHeights", "UCLD mean and heights",
                "Scales UCLD mean inversely to node heights of the tree", model.getParameter(ClockType.UCLD_MEAN),
                tree.getParameter("treeModel.allInternalNodeHeights"), OperatorType.UP_DOWN, demoTuning, rateWeights);
        

        // These are statistics which could have priors on...        
        // #meanRate = #Relaxed Clock Model * #Tree Model
        createStatistic("meanRate", "The mean rate of evolution over the whole tree", 0.0, Double.POSITIVE_INFINITY);
        // #covariance = #Relaxed Clock Model * #Tree Model
        createStatistic("covariance", "The covariance in rates of evolution on each lineage with their ancestral lineages",
                Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY);
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
	    getParameter("branchRates.categories");
		getParameter(ClockType.LOCAL_CLOCK + "." + "rates");
		getParameter(ClockType.LOCAL_CLOCK + "." + "changes");
		getParameter("treeModel.rootRate");
		getParameter("treeModel.nodeRates");
		getParameter("treeModel.allRates");
		
		if (options.hasData()) {
            // if not fixed then do mutation rate move and up/down move            
            boolean fixed = !model.isEstimatedRate();
                       
            Parameter rateParam;

            switch (model.getClockType()) { 
            	case AUTOCORRELATED_LOGNORMAL:
					rateParam = getParameter("treeModel.rootRate");
			        rateParam.isFixed = fixed;
			        if (!fixed) params.add(rateParam);
        
			        params.add(getParameter("branchRates.var"));
			        break;
            }
		}	
    }
 
    /**
     * return a list of operators that are required
     *
     * @param ops the operator list
     */
    public void selectOperators(List<Operator> ops) {
        if (options.hasData()) {
        	// for isEstimatedRate() = false, write nothing on up part of upDownOp
            if (options.clockModelOptions.getRateOptionClockModel() == FixRateType.RElATIVE_TO) {//&& model.isEstimatedRate()) {
            	switch (model.getClockType()) {
	                case STRICT_CLOCK:
	                    ops.add(getOperator("upDownRateHeights"));
	                    break;
	
	                case UNCORRELATED_EXPONENTIAL:
	                    ops.add(getOperator("upDownUCEDMeanHeights"));
	                    addBranchRateCategories(ops);
	                    break;
	
	                case UNCORRELATED_LOGNORMAL:
	                    ops.add(getOperator("upDownUCLDMeanHeights"));
	                    addBranchRateCategories(ops);
	                    break;
	
	                case AUTOCORRELATED_LOGNORMAL:
	                	ops.add(getOperator("scaleRootRate"));
	                    ops.add(getOperator("scaleOneRate"));
	                    ops.add(getOperator("scaleAllRates"));
	                    ops.add(getOperator("scaleAllRatesIndependently"));
	                    ops.add(getOperator("branchRates.var"));
//	                	ops.add(getOperator("upDownAllRatesHeights"));
	                    break;
	
	                case RANDOM_LOCAL_CLOCK:
	                    ops.add(getOperator("upDownRateHeights"));
	                    addRandomLocalClockOperators(ops);
	                    break;
	
	                default:
	                    throw new IllegalArgumentException("Unknown clock model");
	            }
            } else {
            	switch (model.getClockType()) {
	                case STRICT_CLOCK:
	                    // no parameter to operator on
	                    break;
	
	                case UNCORRELATED_EXPONENTIAL:  
	                case UNCORRELATED_LOGNORMAL:
	                    addBranchRateCategories(ops);                       
	                    break;
	
	                case AUTOCORRELATED_LOGNORMAL:                        
	                    ops.add(getOperator("scaleOneRate"));
	                    ops.add(getOperator("scaleAllRatesIndependently"));                        
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
        ops.add(getOperator("swapBranchRateCategories"));
        ops.add(getOperator("randomWalkBranchRateCategories"));
        ops.add(getOperator("uniformBranchRateCategories"));        
    }

    private void addRandomLocalClockOperators(List<Operator> ops) {
    	ops.add(getOperator(ClockType.LOCAL_CLOCK + ".rates"));
        ops.add(getOperator(ClockType.LOCAL_CLOCK + ".changes"));
                   
//TODO    ops.add(tree.getOperator("treeBitMove"));
        
    }

    /**
     * return a list of parameters that are required
     *
     * @param params the parameter list
     */
    public void selectStatistics(List<Parameter> params) {

//        if (options.taxonSets != null) {
//            for (Taxa taxonSet : options.taxonSets) {
//                Parameter statistic = statistics.get(taxonSet);
//                if (statistic == null) {
//                    statistic = new Parameter(taxonSet, "tMRCA for taxon set ");
//                    statistics.put(taxonSet, statistic);
//                }
//                params.add(statistic);
//            }
//        } else {
//            System.err.println("TaxonSets are null");
//        }
        
    	// Statistics
        if (model.getClockType() != ClockType.STRICT_CLOCK) {
            params.add(getParameter("meanRate"));
            params.add(getParameter("covariance"));
            params.add(getParameter(RateStatistic.COEFFICIENT_OF_VARIATION));
        }
                
    	//TODO ?        
        if (model.getClockType() == ClockType.RANDOM_LOCAL_CLOCK) {
            if (localClockRateChangesStatistic == null) {
                localClockRateChangesStatistic = new Parameter("rateChanges", "number of random local clocks", true);
                localClockRateChangesStatistic.priorType = PriorType.POISSON_PRIOR;
                localClockRateChangesStatistic.poissonMean = 1.0;
                localClockRateChangesStatistic.poissonOffset = 0.0;
            }
            if (localClockRatesStatistic == null) {
                localClockRatesStatistic = new Parameter(ClockType.LOCAL_CLOCK + "." + "rates", "random local clock rates", false);

                localClockRatesStatistic.priorType = PriorType.GAMMA_PRIOR;
                localClockRatesStatistic.gammaAlpha = 0.5;
                localClockRatesStatistic.gammaBeta = 2.0;
            }

            localClockRateChangesStatistic.setPrefix(getPrefix());
            params.add(localClockRateChangesStatistic);
            localClockRatesStatistic.setPrefix(getPrefix());
            params.add(localClockRatesStatistic);
        }

//	        if (clock.getClockType() != ClockType.STRICT_CLOCK) {
//	            params.add(getParameter("meanRate"));
//	            params.add(getParameter(RateStatistic.COEFFICIENT_OF_VARIATION));
//	            params.add(getParameter("covariance"));
//	        }
                
    }
    /////////////////////////////////////////////////////////////
    public PartitionClockModel getPartitionClockModel() {
		return model;
	}

	public PartitionTreeModel getPartitionTreeTree() {
		return tree;
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
        return model.getPrefix() + tree.getPrefix();
    }
}