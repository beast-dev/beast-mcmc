package dr.app.beauti.options;

import java.util.List;

import dr.app.beauti.priorsPanel.PriorType;
import dr.evomodel.tree.RateStatistic;

/**
 * @author Andrew Rambaut
 * @author Walter Xie
 * @version $Id$
 */
public class PartitionClockModel extends ModelOptions {
	
	// Instance variables

    private final BeautiOptions options;
    private String name;
    private PartitionData partition;
    
    private ClockType clockType = ClockType.STRICT_CLOCK;
    
    public Parameter localClockRateChangesStatistic = null;
    public Parameter localClockRatesStatistic = null;
    
    
	public PartitionClockModel(BeautiOptions options, PartitionData partition) {
		this.options = options;
		this.partition = partition;
		this.name = partition.getName();
		
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
        
        String treePrefix = partition.getPartitionTreeModel().getPrefix();        

        createOperator("upDownRateHeights", "Substitution rate and heights",
                "Scales substitution rates inversely to node heights of the tree", "clock.rate",
                treePrefix + "treeModel.allInternalNodeHeights", OperatorType.UP_DOWN, 0.75, rateWeights);
        createOperator("upDownUCEDMeanHeights", "UCED mean and heights",
                "Scales UCED mean inversely to node heights of the tree", ClockType.UCED_MEAN,
                treePrefix + "treeModel.allInternalNodeHeights", OperatorType.UP_DOWN, 0.75, rateWeights);
        createOperator("upDownUCLDMeanHeights", "UCLD mean and heights",
                "Scales UCLD mean inversely to node heights of the tree", ClockType.UCLD_MEAN,
                treePrefix + "treeModel.allInternalNodeHeights", OperatorType.UP_DOWN, 0.75, rateWeights);
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
                    if (!fixed) params.add(rateParam);
                    break;

                case UNCORRELATED_EXPONENTIAL:
                    rateParam = getParameter(ClockType.UCED_MEAN);
                    if (!fixed) params.add(rateParam);
                    break;

                case UNCORRELATED_LOGNORMAL:
                    rateParam = getParameter(ClockType.UCLD_MEAN);
                    if (!fixed) params.add(rateParam);
                    params.add(getParameter(ClockType.UCLD_STDEV));
                    break;

                case AUTOCORRELATED_LOGNORMAL:
//                	rateParam = partition.getPartitionTreeModel().getParameter("treeModel.rootRate");
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
                switch (options.clockType) {
                    case STRICT_CLOCK:
                        ops.add(getOperator("clock.rate"));
                        ops.add(getOperator("upDownRateHeights"));
                        break;

                    case UNCORRELATED_EXPONENTIAL:
                        ops.add(getOperator(ClockType.UCED_MEAN));
                        ops.add(getOperator("upDownUCEDMeanHeights"));
                        ops.add(getOperator("swapBranchRateCategories"));
                        ops.add(getOperator("randomWalkBranchRateCategories"));
                        ops.add(getOperator("unformBranchRateCategories"));
                        break;

                    case UNCORRELATED_LOGNORMAL:
                        ops.add(getOperator(ClockType.UCLD_MEAN));
                        ops.add(getOperator(ClockType.UCLD_STDEV));
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
                        ops.add(getOperator(ClockType.LOCAL_CLOCK + "." + "rates"));
                        ops.add(getOperator(ClockType.LOCAL_CLOCK + "." + "changes"));
                        ops.add(getOperator("treeBitMove"));
                        break;

                    default:
                        throw new IllegalArgumentException("Unknown clock model");
                }
            } else {
                switch (options.clockType) {
                    case STRICT_CLOCK:
                        // no parameter to operator on
                        break;

                    case UNCORRELATED_EXPONENTIAL:
                        ops.add(getOperator("swapBranchRateCategories"));
                        ops.add(getOperator("randomWalkBranchRateCategories"));
                        ops.add(getOperator("unformBranchRateCategories"));
                        break;

                    case UNCORRELATED_LOGNORMAL:
                        ops.add(getOperator(ClockType.UCLD_STDEV));
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
                        ops.add(getOperator(ClockType.LOCAL_CLOCK + "." + "rates"));
                        ops.add(getOperator(ClockType.LOCAL_CLOCK + "." + "changes"));
                        ops.add(getOperator("treeBitMove"));
                        break;

                    default:
                        throw new IllegalArgumentException("Unknown clock model");
                }
            }
        }
    }
    
    // use override method getParameter(String name) and getOperator(String name) in PartitionModel containing prefix
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

        if (clockType == ClockType.RANDOM_LOCAL_CLOCK) {
            if (this.localClockRateChangesStatistic == null) {
            	this.localClockRateChangesStatistic = new Parameter("rateChanges", "number of random local clocks", true);
            	this.localClockRateChangesStatistic.priorType = PriorType.POISSON_PRIOR;
            	this.localClockRateChangesStatistic.poissonMean = 1.0;
            	this.localClockRateChangesStatistic.poissonOffset = 0.0;
            }
            if (this.localClockRatesStatistic == null) {
            	this.localClockRatesStatistic = new Parameter(ClockType.LOCAL_CLOCK + "." + "rates", "random local clock rates", false);

            	this.localClockRatesStatistic.priorType = PriorType.GAMMA_PRIOR;
            	this.localClockRatesStatistic.gammaAlpha = 0.5;
            	this.localClockRatesStatistic.gammaBeta = 2.0;
            }
            
            this.localClockRateChangesStatistic.setPrefix(getPrefix());
            params.add(this.localClockRateChangesStatistic);
            this.localClockRatesStatistic.setPrefix(getPrefix());
            params.add(this.localClockRatesStatistic);
        }

        if (clockType != ClockType.STRICT_CLOCK) {
            params.add(getParameter("meanRate"));
            params.add(getParameter(RateStatistic.COEFFICIENT_OF_VARIATION));
            params.add(getParameter("covariance"));
        }

    }
    
    /////////////////////////////////////////////////////////////
    

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

        if (name.startsWith(getName())) {
            name = name.substring(getName().length() + 1);
        }
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

        operator.setPrefix(getName());

        return operator;
    }


    public String getPrefix() {
        String prefix = "";
        if (options.getActivePartitionTreeModels().size() > 1) { //|| options.isSpeciesAnalysis()
            // There is more than one active partition model
            prefix += getName() + ".";
        }
        return prefix;
    }


}
