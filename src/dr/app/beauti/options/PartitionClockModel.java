package dr.app.beauti.options;

import java.util.ArrayList;
import java.util.List;

import dr.app.beauti.priorsPanel.PriorType;

/**
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
        
        for (PartitionTreeModel tree : options.getPartitionTreeModels(getAllPartitionData())) { // borrow the method in BeautiOption        	
            createOperator("upDownRateHeights", "Substitution rate and heights",
                    "Scales substitution rates inversely to node heights of the tree", super.getParameter("clock.rate"),
                    tree.getParameter("treeModel.allInternalNodeHeights"), OperatorType.UP_DOWN, 0.75, rateWeights);
            createOperator("upDownUCEDMeanHeights", "UCED mean and heights",
                    "Scales UCED mean inversely to node heights of the tree", super.getParameter(ClockType.UCED_MEAN),
                    tree.getParameter("treeModel.allInternalNodeHeights"), OperatorType.UP_DOWN, 0.75, rateWeights);
            createOperator("upDownUCLDMeanHeights", "UCLD mean and heights",
                    "Scales UCLD mean inversely to node heights of the tree", super.getParameter(ClockType.UCLD_MEAN),
                    tree.getParameter("treeModel.allInternalNodeHeights"), OperatorType.UP_DOWN, 0.75, rateWeights);
        }        
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
                        ops.add(getOperator("upDownRateHeights"));
                        break;

                    case UNCORRELATED_EXPONENTIAL:
                        ops.add(getOperator(ClockType.UCED_MEAN));
                        ops.add(getOperator("upDownUCEDMeanHeights"));
                        for (PartitionTreeModel tree : options.getPartitionTreeModels(getAllPartitionData())) {
	                        ops.add(tree.getOperator("swapBranchRateCategories"));
	                        ops.add(tree.getOperator("randomWalkBranchRateCategories"));
	                        ops.add(tree.getOperator("unformBranchRateCategories"));
                        }
                        break;

                    case UNCORRELATED_LOGNORMAL:
                        ops.add(getOperator(ClockType.UCLD_MEAN));
                        ops.add(getOperator(ClockType.UCLD_STDEV));
                        ops.add(getOperator("upDownUCLDMeanHeights"));
                        for (PartitionTreeModel tree : options.getPartitionTreeModels(getAllPartitionData())) {
	                        ops.add(tree.getOperator("swapBranchRateCategories"));
	                        ops.add(tree.getOperator("randomWalkBranchRateCategories"));
	                        ops.add(tree.getOperator("unformBranchRateCategories"));
                        }
                        break;

                    case AUTOCORRELATED_LOGNORMAL:
                    	for (PartitionTreeModel tree : options.getPartitionTreeModels(getAllPartitionData())) {
                    		ops.add(tree.getOperator("scaleRootRate"));
                            ops.add(tree.getOperator("scaleOneRate"));
                            ops.add(tree.getOperator("scaleAllRates"));
                            ops.add(tree.getOperator("scaleAllRatesIndependently"));
                            ops.add(tree.getOperator("upDownAllRatesHeights"));
                        }                        
                        ops.add(getOperator("branchRates.var"));
                        break;

                    case RANDOM_LOCAL_CLOCK:
                        ops.add(getOperator("clock.rate"));
                        ops.add(getOperator("upDownRateHeights"));
                        for (PartitionTreeModel tree : options.getPartitionTreeModels(getAllPartitionData())) {
	                        ops.add(tree.getOperator(ClockType.LOCAL_CLOCK + "." + "rates"));
	                        ops.add(tree.getOperator(ClockType.LOCAL_CLOCK + "." + "changes"));                        
                        	ops.add(tree.getOperator("treeBitMove"));                            
                        }
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
                    	for (PartitionTreeModel tree : options.getPartitionTreeModels(getAllPartitionData())) {
	                        ops.add(tree.getOperator("swapBranchRateCategories"));
	                        ops.add(tree.getOperator("randomWalkBranchRateCategories"));
	                        ops.add(tree.getOperator("unformBranchRateCategories"));
                        }
                        break;

                    case UNCORRELATED_LOGNORMAL:
                    	for (PartitionTreeModel tree : options.getPartitionTreeModels(getAllPartitionData())) {
	                        ops.add(tree.getOperator("swapBranchRateCategories"));
	                        ops.add(tree.getOperator("randomWalkBranchRateCategories"));
	                        ops.add(tree.getOperator("unformBranchRateCategories"));
                        }
                        break;

                    case AUTOCORRELATED_LOGNORMAL:
                        for (PartitionTreeModel tree : options.getPartitionTreeModels(getAllPartitionData())) {
                    		ops.add(tree.getOperator("scaleOneRate"));
                            ops.add(tree.getOperator("scaleAllRatesIndependently"));                            
                        } 
                        ops.add(getOperator("branchRates.var"));
                        break;

                    case RANDOM_LOCAL_CLOCK:
                    	for (PartitionTreeModel tree : options.getPartitionTreeModels(getAllPartitionData())) {
	                        ops.add(tree.getOperator(ClockType.LOCAL_CLOCK + "." + "rates"));
	                        ops.add(tree.getOperator(ClockType.LOCAL_CLOCK + "." + "changes"));                        
                        	ops.add(tree.getOperator("treeBitMove"));                            
                        }
                        break;

                    default:
                        throw new IllegalArgumentException("Unknown clock model");
                }
            }
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
        if (options.getPartitionTreeModels().size() > 1) { //|| options.isSpeciesAnalysis()
            // There is more than one active partition model
            prefix += getName() + ".";
        }
        return prefix;
    }


}
