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

    private List<PartitionData> allPartitionData = new ArrayList<PartitionData>();

    private ClockType clockType = ClockType.STRICT_CLOCK;
    private boolean isFixedRate = true;

    public PartitionClockModel(BeautiOptions options, PartitionData partition) {
        this.options = options;
        this.name = partition.getName();

        allPartitionData.clear();
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

    private void initClockModelParaAndOpers() {
        createParameter("clock.rate", "substitution rate", SUBSTITUTION_RATE_SCALE, 1.0, 0.0, Double.POSITIVE_INFINITY);
        createParameter(ClockType.UCED_MEAN, "uncorrelated exponential relaxed clock mean", SUBSTITUTION_RATE_SCALE, 1.0, 0.0, Double.POSITIVE_INFINITY);
        createParameter(ClockType.UCLD_MEAN, "uncorrelated lognormal relaxed clock mean", SUBSTITUTION_RATE_SCALE, 1.0, 0.0, Double.POSITIVE_INFINITY);
        createParameter(ClockType.UCLD_STDEV, "uncorrelated lognormal relaxed clock stdev", LOG_STDEV_SCALE, 0.1, 0.0, Double.POSITIVE_INFINITY);
        
        createScaleOperator("clock.rate", demoTuning, rateWeights);
        createScaleOperator(ClockType.UCED_MEAN, demoTuning, rateWeights);
        createScaleOperator(ClockType.UCLD_MEAN, demoTuning, rateWeights);
        createScaleOperator(ClockType.UCLD_STDEV, demoTuning, rateWeights);
    }

    /**
     * return a list of parameters that are required
     *
     * @param params the parameter list
     */
    public void selectParameters(List<Parameter> params) {    	    	
        if (options.hasData()) {
            // if not fixed then do mutation rate move and up/down move
            boolean fixed = isFixedRate;
            Parameter rateParam;

            switch (clockType) {
                case STRICT_CLOCK:
                    rateParam = getParameter("clock.rate");
                    rateParam.isFixed = fixed;
                    if (fixed) rateParam.initial = options.getMeanSubstitutionRate();
                    if (!fixed) params.add(rateParam);
                    break;

                case UNCORRELATED_EXPONENTIAL:
                    rateParam = getParameter(ClockType.UCED_MEAN);
                    rateParam.isFixed = fixed;
                    if (fixed) rateParam.initial = options.getMeanSubstitutionRate();
                    if (!fixed) params.add(rateParam);
                    break;

                case UNCORRELATED_LOGNORMAL:
                    rateParam = getParameter(ClockType.UCLD_MEAN);
                    rateParam.isFixed = fixed;
                    if (fixed) rateParam.initial = options.getMeanSubstitutionRate();
                    if (!fixed) params.add(rateParam);
                    params.add(getParameter(ClockType.UCLD_STDEV));
                    break;

                case AUTOCORRELATED_LOGNORMAL:                    
//                    rateParam = getParameter("treeModel.rootRate");//TODO fix tree?
//                    rateParam.isFixed = fixed;
//                    if (!fixed) params.add(rateParam);
//                    
//                    params.add(getParameter("branchRates.var"));
                    break;

                case RANDOM_LOCAL_CLOCK:
                    rateParam = getParameter("clock.rate");
                    rateParam.isFixed = fixed;
                    if (fixed) rateParam.initial = options.getMeanSubstitutionRate();
                    if (!fixed) params.add(rateParam);
                    break;

                default:
                    throw new IllegalArgumentException("Unknown clock model");
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

            if (isFixedRate) {
            	switch (clockType) {
	                case STRICT_CLOCK:	
	                case UNCORRELATED_EXPONENTIAL:
	                case AUTOCORRELATED_LOGNORMAL: 
	                case RANDOM_LOCAL_CLOCK:
	                	// no parameter to operator on
	                    break;      
	
	                case UNCORRELATED_LOGNORMAL:	
	                    ops.add(getOperator(ClockType.UCLD_STDEV));
	                    break;
	
	                default:
	                    throw new IllegalArgumentException("Unknown clock model");
            	}
            } else {                
                switch (clockType) {
	                case STRICT_CLOCK:
	                    ops.add(getOperator("clock.rate"));	
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
	
	                case RANDOM_LOCAL_CLOCK:
	                    ops.add(getOperator("clock.rate"));
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

	public void setFixedRate(boolean isFixedRate) {
		this.isFixedRate = isFixedRate;
	}

	public boolean isFixedRate() {
		return isFixedRate;
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
