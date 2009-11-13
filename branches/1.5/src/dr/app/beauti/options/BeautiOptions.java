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

import dr.app.beauti.components.ComponentFactory;
import dr.app.beauti.mcmcpanel.MCMCPanel;
import dr.app.beauti.enumTypes.PriorScaleType;
import dr.app.beauti.enumTypes.TreePriorType;
import dr.app.beauti.util.BEAUTiImporter;
import dr.app.beauti.util.BeautiTemplate;
import dr.evolution.datatype.DataType;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxa;
import dr.evolution.util.Units;
import dr.evolution.alignment.Alignment;
import dr.evoxml.DateUnitsType;
import dr.inference.operators.OperatorSchedule;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 */
public class BeautiOptions extends ModelOptions {

    public BeautiOptions() {
        this(new ComponentFactory[]{});       
    }

    public BeautiOptions(ComponentFactory[] components) {    	

        // Install all the component's options from the given list of factories:
        for (ComponentFactory component : components) {
            addComponent(component.getOptions(this));
        }
    }   

    /**
     * resets the options to the initial conditions
     */
    public void reset() {
        // Data options
        allowDifferentTaxa = false;
        dataType = null;
//        dataReset = true;

        taxonList = null;
        taxonSets.clear();
        taxonSetsMono.clear();
        
//        meanDistance = 1.0;
        datesUnits = DateUnitsType.YEARS;
        datesDirection = DateUnitsType.FORWARDS;
        maximumTipHeight = 0.0;
        translation = 0;
        
//        selecetedTraits.clear();
//        traitTypes.clear();
        
        dataPartitions.clear();
//        partitionModels.clear();
//        partitionTreeModels.clear();
//        partitionTreePriors.clear();
        partitionClockTreeLinks.clear();
//        activedSameTreePrior = null;
//        shareSameTreePrior = true;
        userTrees.clear();

//        rateOptionClockModel = FixRateType.FIX_FIRST_PARTITION;
//        meanSubstitutionRate = 1.0;
        unlinkPartitionRates = true;

        units = Units.Type.SUBSTITUTIONS;

        // Operator schedule options
        coolingSchedule = OperatorSchedule.DEFAULT_SCHEDULE;

        // MCMC options
        chainLength = 10000000;
        logEvery = 1000;
        echoEvery = 1000;
        burnIn = 100000;
        fileName = null;
        autoOptimize = true;
        performTraceAnalysis = false;
        generateCSV = true;  // until/if a button
        samplePriorOnly = false;
        
        fileNameStem = MCMCPanel.fileNameStem;
        logFileName = null;
//        mapTreeLog = false;
//        mapTreeFileName = null;
        treeFileName.clear();
        substTreeLog = false;
        substTreeFileName.clear();
        operatorAnalysis = false;
        operatorAnalysisFileName = null;
        
//        substitutionModelOptions = new SubstitutionModelOptions(this);
        clockModelOptions = new ClockModelOptions(this);
        treeModelOptions = new TreeModelOptions(this);
        priorOptions = new PriorOptions(this);
        
        traitOptions = new TraitOptions(this);
        starBEASTOptions = new STARBEASTOptions(this);
        
        beautiImporter = new BEAUTiImporter(this);    
        
        beautiTemplate = new BeautiTemplate(this);
        
        parameters.clear();
        operators.clear();
        statistics.clear();
    }
    
    public void selectTaxonSetsStatistics(List<Parameter> params) {
    	
        if (taxonSets != null) {
            for (Taxa taxa : taxonSets) {
                Parameter statistic = statistics.get(taxa);
                if (statistic == null) {
                    statistic = new Parameter.Builder(taxa.getId(), "").taxa(taxa)
                            .isStatistic(true).isNodeHeight(true).scaleType(PriorScaleType.TIME_SCALE)
                            .lower(0.0).upper(Double.MAX_VALUE).build();
                    statistics.put(taxa, statistic);
                }
                params.add(statistic);
            }
        } else {
            System.err.println("TaxonSets are null");
        }
    }
    
    /**
     * return an list of parameters that are required
     *
     * @return the parameter list
     */
    public ArrayList<Parameter> selectParameters() {

        ArrayList<Parameter> parameters = new ArrayList<Parameter>();
        
        selectTaxonSetsStatistics(parameters); // have to be before clockModelOptions.selectParameters(parameters);       
        
        for (PartitionSubstitutionModel model : getPartitionSubstitutionModels()) {
//          parameters.addAll(model.getParameters(multiplePartitions));
        	model.selectParameters(parameters);
        }
//        substitutionModelOptions.selectParameters(parameters);
        
        for (PartitionClockModel model : getPartitionClockModels()) {
            model.selectParameters(parameters);
        }
        clockModelOptions.selectParameters();

        for (PartitionTreeModel tree : getPartitionTreeModels()) {
            tree.selectParameters(parameters);
        }
        treeModelOptions.selectParameters(parameters);

        for (PartitionTreePrior prior : getPartitionTreePriors()) {
            prior.selectParameters(parameters);
        }

        for (PartitionClockModelTreeModelLink clockTree : getPartitionClockTreeLinks()) {
        	clockTree.selectParameters(parameters);
            clockTree.selectStatistics(parameters);
        }
        
        if (starBEASTOptions.isSpeciesAnalysis()) { // species
        	starBEASTOptions.selectParameters(parameters);
        }

        selectComponentParameters(this, parameters);
        
        selectComponentStatistics(this, parameters);
        
        priorOptions.selectParameters(parameters);

//        boolean multiplePartitions = getTotalActivePartitionSubstitutionModelCount() > 1;
        // add all Parameter (with prefix) into parameters list     
        

//
//        double growthRateMaximum = 1E6;
//        double birthRateMaximum = 1E6;
//        double substitutionRateMaximum = 100;
//        double logStdevMaximum = 10;
//        double substitutionParameterMaximum = 100;
//        double initialRootHeight = 1;
//        double initialRate = 1;
//
//
//        if (clockModelOptions.getRateOptionClockModel() == FixRateType.FIX_MEAN) {
//            double rate = clockModelOptions.getMeanRelativeRate();
//
//            growthRateMaximum = 1E6 * rate;
//            birthRateMaximum = 1E6 * rate;
//
//            if (hasData()) {
//                initialRootHeight = meanDistance / rate;
//
//                initialRootHeight = round(initialRootHeight, 2);
//            }
//
//        } else {
//            if (maximumTipHeight > 0) {
//                initialRootHeight = maximumTipHeight * 10.0;
//            }
//
//            initialRate = round((meanDistance * 0.2) / initialRootHeight, 2);
//        }
//
//        double timeScaleMaximum = round(initialRootHeight * 1000.0, 2);
//
//        for (Parameter param : parameters) {
////            if (dataReset) param.priorEdited = false;
//
//            if (!param.priorEdited) {
//                switch (param.scale) {
//                    case TIME_SCALE:
//                        param.uniformLower = Math.max(0.0, param.lower);
//                        param.uniformUpper = Math.min(timeScaleMaximum, param.upper);
//                        param.initial = initialRootHeight;
//                        break;
//                    case T50_SCALE:
//                        param.uniformLower = Math.max(0.0, param.lower);
//                        param.uniformUpper = Math.min(timeScaleMaximum, param.upper);
//                        param.initial = initialRootHeight / 5.0;
//                        break;
//                    case GROWTH_RATE_SCALE:
//                        param.uniformLower = Math.max(-growthRateMaximum, param.lower);
//                        param.uniformUpper = Math.min(growthRateMaximum, param.upper);
//                        break;
//                    case BIRTH_RATE_SCALE:
//                        param.uniformLower = Math.max(0.0, param.lower);
//                        param.uniformUpper = Math.min(birthRateMaximum, param.upper);
//                        break;
//                    case SUBSTITUTION_RATE_SCALE:
//                        param.uniformLower = Math.max(0.0, param.lower);
//                        param.uniformUpper = Math.min(substitutionRateMaximum, param.upper);
//                        param.initial = initialRate;
//                        break;
//                    case LOG_STDEV_SCALE:
//                        param.uniformLower = Math.max(0.0, param.lower);
//                        param.uniformUpper = Math.min(logStdevMaximum, param.upper);
//                        break;
//                    case SUBSTITUTION_PARAMETER_SCALE:
//                        param.uniformLower = Math.max(0.0, param.lower);
//                        param.uniformUpper = Math.min(substitutionParameterMaximum, param.upper);
//                        break;
//
//                    case UNITY_SCALE:
//                        param.uniformLower = 0.0;
//                        param.uniformUpper = 1.0;
//                        break;
//
//                    case ROOT_RATE_SCALE:
//                        param.initial = initialRate;
//                        param.gammaAlpha = 0.5;
//                        param.gammaBeta = param.initial / 0.5;
//                        break;
//
//                    case LOG_VAR_SCALE:
//                        param.initial = initialRate;
//                        param.gammaAlpha = 2.0;
//                        param.gammaBeta = param.initial / 2.0;
//                        break;
//
//                }
//                if (param.isNodeHeight) {
//                    param.lower = maximumTipHeight;
//                    param.uniformLower = maximumTipHeight;
//                    param.uniformUpper = timeScaleMaximum;
//                    param.initial = initialRootHeight;
//                }
//            }
//        }
//
////        dataReset = false;

        return parameters;
    }

    /**
     * return an list of operators that are required
     *
     * @return the operator list
     */
    public List<Operator> selectOperators() {

        ArrayList<Operator> ops = new ArrayList<Operator>();

        for (PartitionSubstitutionModel model : getPartitionSubstitutionModels()) {
        	model.selectOperators(ops);
        }
//        substitutionModelOptions.selectOperators(ops);
        
        for (PartitionClockModel model : getPartitionClockModels()) {
            model.selectOperators(ops);
        }
        clockModelOptions.selectOperators(ops);

        for (PartitionTreeModel tree : getPartitionTreeModels()) {
            tree.selectOperators(ops);
        }
        treeModelOptions.selectOperators(ops);

        for (PartitionTreePrior prior : getPartitionTreePriors()) {
            prior.selectOperators(ops);
        }

        for (PartitionClockModelTreeModelLink clockTree : getPartitionClockTreeLinks()) {
            clockTree.selectOperators(ops);
        }

        if (starBEASTOptions.isSpeciesAnalysis()) { // species
        	starBEASTOptions.selectOperators(ops);
        }

        selectComponentOperators(this, ops);

//        priorOptions.selectOperators(ops);
        
//        if (multiplePartitions) {
//        if (hasCodon()) {
//            Operator deltaMuOperator = getOperator("deltaMu");
//
//            // update delta mu operator weight
//            deltaMuOperator.weight = 0.0;
//            for (PartitionSubstitutionModel pm : getPartitionSubstitutionModels()) {
//                deltaMuOperator.weight += pm.getCodonPartitionCount();
//            }
//
//            ops.add(deltaMuOperator);
//        }

//        double initialRootHeight = 1;
//
//        if (clockModelOptions.getRateOptionClockModel() == FixRateType.FIX_MEAN) {
//            double rate = clockModelOptions.getMeanRelativeRate();
//
//            if (hasData()) {
//                initialRootHeight = meanDistance / rate;
//                initialRootHeight = round(initialRootHeight, 2);
//            }
//
//        } else {
//            if (maximumTipHeight > 0) {
//                initialRootHeight = maximumTipHeight * 10.0;
//            }
//        }

//        for (PartitionTreeModel tree : getPartitionTreeModels()) {
//            Operator op = tree.getOperator("subtreeSlide");
//            if (!op.tuningEdited) {
//                op.tuning = initialRootHeight / 10.0;
//            }
//        }
        
//        for (Operator op : ops) {
//        	System.out.println(op.prefix + " + " + op.getName());
//        }

        return ops;
    }
    
    public boolean hasData() {
        return dataPartitions.size() > 0;
    }

    public PartitionData getPartitionData(Alignment alignment) {
        for (PartitionData pd : dataPartitions) {
            if (pd.getAlignment() == alignment) {
                return pd;
            }
        }
        return null;
    }
    
//    public boolean isFixedSubstitutionRate() {
//        return fixedSubstitutionRate;
//    }
//    public void updateFixedRateClockModel() {
//    	if (getPartitionClockModels().size() > 0) {
//    		
//	    	if (rateOptionClockModel == FixRateType.FIX_FIRST_PARTITION) {
//	    		// fix rate of 1st partition
//	    		for (PartitionClockModel model : getPartitionClockModels()) {
//	    			if (getPartitionClockModels().indexOf(model) < 1) {
//	    				model.setFixedRate(true);
//	    			} else {
//	    				model.setFixedRate(false);
//	    			}
//	            }
//	    		
//	    	} else if (rateOptionClockModel == FixRateType.FIX_MEAN) {
//	    		// TODO check
//	    		for (PartitionClockModel model : getPartitionClockModels()) {
//	    			model.setFixedRate(true);
//	            }
//	    		
//	    	} else {
//	    		// estimate all rate
//	    		for (PartitionClockModel model : getPartitionClockModels()) {
//	    			model.setFixedRate(false);
//	            }
//	    		
//	    	}
//    	}
//    }
    


    public boolean isEBSPSharingSamePrior() {
        return getPartitionTreePriors().size() >= 1 && 
                (isShareSameTreePrior() && getPartitionTreePriors().get(0).getNodeHeightPrior() == TreePriorType.EXTENDED_SKYLINE);
    }

    // ++++++++++++++ Partition Substitution Model ++++++++++++++ 
//    public void addPartitionSubstitutionModel(PartitionSubstitutionModel model) {
//        if (!partitionModels.contains(model)) {
//            partitionModels.add(model);
//        }
//    }   

//    public List<PartitionSubstitutionModel> getPartitionSubstitutionModels() {
//        return partitionModels;
//    }

    public List<PartitionSubstitutionModel> getPartitionSubstitutionModels(DataType dataType, List<PartitionData> givenDataPartitions) {
        List<PartitionSubstitutionModel> models = new ArrayList<PartitionSubstitutionModel>();
        for (PartitionSubstitutionModel model : getPartitionSubstitutionModels(givenDataPartitions)) {
            if (model.getDataType() == dataType) {
                models.add(model);
            }
        }
        return models;
    }

    public List<PartitionSubstitutionModel> getPartitionSubstitutionModels(List<PartitionData> givenDataPartitions) {

        List<PartitionSubstitutionModel> activeModels = new ArrayList<PartitionSubstitutionModel>();

        for (PartitionData partition : givenDataPartitions) {
            PartitionSubstitutionModel model = partition.getPartitionSubstitutionModel();
            if (model != null && (!activeModels.contains(model))) {
                activeModels.add(model);
            }
        }

        return activeModels;
    }

    public List<PartitionSubstitutionModel> getPartitionSubstitutionModels() {
        return getPartitionSubstitutionModels(dataPartitions);
    }

//    public int getTotalActivePartitionSubstitutionModelCount() {
//        int totalPartitionCount = 0;
//        for (PartitionSubstitutionModel model : getPartitionSubstitutionModels()) {
//            totalPartitionCount += model.getCodonPartitionCount();
//        }
//        return totalPartitionCount;
//    }

 
    // ++++++++++++++ Partition Clock Model ++++++++++++++    
//    public void addPartitionClockModel (PartitionClockModel model) {
//        if (!clockModels.contains(model)) {
//        	clockModels.add(model);
//        }
//    }

//    public List<PartitionClockModel> getPartitionClockModels() {
//        return clockModels;
//    }

    public List<PartitionClockModel> getPartitionClockModels(List<PartitionData> givenDataPartitions) {

        List<PartitionClockModel> activeModels = new ArrayList<PartitionClockModel>();

        for (PartitionData partition : givenDataPartitions) {
            PartitionClockModel model = partition.getPartitionClockModel();
            if (model != null && (!activeModels.contains(model))) {
                activeModels.add(model);
            }
        }

        return activeModels;
    }

    public List<PartitionClockModel> getPartitionClockModels() {
        return getPartitionClockModels(dataPartitions);
    }

    // ++++++++++++++ Partition Tree Model ++++++++++++++ 
//    public void addPartitionTreeModel(PartitionTreeModel tree) {
//
//        if (!partitionTreeModels.contains(tree)) {
//            partitionTreeModels.add(tree);
//        }
//    }

//    public List<PartitionTreeModel> getPartitionTreeModels() {
//        return partitionTreeModels;
//    }

    public List<PartitionTreeModel> getPartitionTreeModels(List<PartitionData> givenDataPartitions) {

        List<PartitionTreeModel> activeTrees = new ArrayList<PartitionTreeModel>();

        for (PartitionData partition : givenDataPartitions) {
            PartitionTreeModel tree = partition.getPartitionTreeModel();
            if (tree != null && (!activeTrees.contains(tree))) {
                activeTrees.add(tree);
            }
        }

        return activeTrees;
    }

    public List<PartitionTreeModel> getPartitionTreeModels() {
        return getPartitionTreeModels(dataPartitions);
    }

    // ++++++++++++++ Partition Tree Prior ++++++++++++++ 
    public List<PartitionTreePrior> getPartitionTreePriors() {

        List<PartitionTreePrior> activeTrees = new ArrayList<PartitionTreePrior>();

        // # tree prior = 1 or # tree model
		for (PartitionTreeModel model : getPartitionTreeModels()) {
			PartitionTreePrior prior = model.getPartitionTreePrior();
			if (prior != null && (!activeTrees.contains(prior))) {
				activeTrees.add(prior);
			}			
		}

        return activeTrees;
    }
    
    public void unLinkTreePriors() {
    	for (PartitionTreeModel model : getPartitionTreeModels()) {
    		PartitionTreePrior prior = model.getPartitionTreePrior();    		
    		if (prior == null || (!prior.getName().equals(model.getName()))) {
    			PartitionTreePrior ptp = new PartitionTreePrior(this, model);
    			model.setPartitionTreePrior(ptp);
    		}    		
		}  
    }
    
    public void linkTreePriors(PartitionTreePrior treePrior) {
    	if (treePrior == null) treePrior = new PartitionTreePrior(this, getPartitionTreeModels().get(0));
    	for (PartitionTreeModel model : getPartitionTreeModels()) {
			model.setPartitionTreePrior(treePrior);		
		}    	
    }
    
    public boolean isShareSameTreePrior() {
    	return getPartitionTreePriors().size() <= 1;
    }
    
    // ++++++++++++++ Partition Clock Model ++++++++++++++    
    public List<PartitionClockModelTreeModelLink> getPartitionClockTreeLinks() {    	
    	return partitionClockTreeLinks;
    }
    
    public PartitionClockModelTreeModelLink getPartitionClockTreeLink(PartitionClockModel model, PartitionTreeModel tree) {    	
    	for (PartitionClockModelTreeModelLink clockTree : getPartitionClockTreeLinks()) {
        	if (clockTree.getPartitionClockModel().equals(model) && clockTree.getPartitionTreeTree().equals(tree)) {
        		return clockTree;        		
        	} 
        }
    	
    	return null;
    }
    
    public void updatePartitionClockTreeLinks() {
    	partitionClockTreeLinks.clear();
    	
    	for (PartitionClockModel model : getPartitionClockModels()) {
            for (PartitionTreeModel tree : getPartitionTreeModels(model.getAllPartitionData())) {
                PartitionClockModelTreeModelLink clockTree = new PartitionClockModelTreeModelLink(this, model, tree);
                
                if (!partitionClockTreeLinks.contains(clockTree)) {
                	partitionClockTreeLinks.add(clockTree);
                }
            }
        }    	
    }

    // update links (e.g List<PartitionData> allPartitionData), after use (e.g partition.setPartitionSubstitutionModel(model))
    public void updateLinksBetweenPDPCMPSMPTMPTPP() {
        for (PartitionSubstitutionModel model : getPartitionSubstitutionModels()) {
            model.clearAllPartitionData();
        }

        for (PartitionClockModel model : getPartitionClockModels()) {
            model.clearAllPartitionData();
        }

        for (PartitionTreeModel tree : getPartitionTreeModels()) {
            tree.clearAllPartitionData();
        }

        //TODO update PartitionTreePrior ?

        for (PartitionData partition : dataPartitions) {
            PartitionSubstitutionModel psm = partition.getPartitionSubstitutionModel();
            if (!psm.getAllPartitionData().contains(partition)) {
                psm.addPartitionData(partition);
            }

            PartitionClockModel pcm = partition.getPartitionClockModel();
            if (!pcm.getAllPartitionData().contains(partition)) {
                pcm.addPartitionData(partition);
            }

            PartitionTreeModel ptm = partition.getPartitionTreeModel();
            if (!ptm.getAllPartitionData().contains(partition)) {
                ptm.addPartitionData(partition);
            }
        }

    }

	public double getAveWeightedMeanDistance(List<PartitionData> partitions) {
		double meanDistance = 0;
		double totalSite = 0;
		for (PartitionData partition : partitions) {
			meanDistance = meanDistance + partition.getMeanDistance() * partition.getSiteCount();
			totalSite = totalSite + partition.getSiteCount();	
		}
		
		if (totalSite == 0) {
			return 0;
		} else {
			return meanDistance / totalSite;
		}
    }
	
	public String statusMessage() {
        String message = "";
        if (hasData()) {
            message += "Data: " + taxonList.getTaxonCount() + " taxa, " +
                    dataPartitions.size() +
                    (dataPartitions.size() > 1 ? " partitions" : " partition");

            if (starBEASTOptions.isSpeciesAnalysis()) {                
                int num = starBEASTOptions.getSpeciesList().size();
                message += ", " + num + " species"; // species is both singular and plural
            }

            if (userTrees.size() > 0) {
                message += ", " + userTrees.size() +
                        (userTrees.size() > 1 ? " trees" : " tree");
            }
            
            if (allowDifferentTaxa) {
                message += " in total";
            }

            if (starBEASTOptions.isSpeciesAnalysis()) {
                message += ";    Species Tree Ancestral Reconstruction (*BEAST)";
            }
            
            message += ";    " + clockModelOptions.statusMessageClockModel();
            
        } else if (userTrees.size() > 0) {
            message += "Trees only : " + userTrees.size() +
                    (userTrees.size() > 1 ? " trees, " : " tree, ") +
                    taxonList.getTaxonCount() + " taxa";
        } else if (taxonList != null && taxonList.getTaxonCount() > 0) {
            message += "Taxa only: " + taxonList.getTaxonCount() + " taxa";
        }

        return message;
    }

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // Data options
    public boolean allowDifferentTaxa = false;
    public DataType dataType = null;
//    public boolean dataReset = true;

    public Taxa taxonList = null;

    public List<Taxa> taxonSets = new ArrayList<Taxa>();
    public Map<Taxa, Boolean> taxonSetsMono = new HashMap<Taxa, Boolean>();

//    public double meanDistance = 1.0;
    public DateUnitsType datesUnits = DateUnitsType.YEARS;
    public DateUnitsType datesDirection = DateUnitsType.FORWARDS;
    public double maximumTipHeight = 0.0;
    public int translation = 0;

    public DateGuesser dateGuesser = new DateGuesser();
//    public TraitGuesser traitGuesser = new TraitGuesser();
//
//    public List<String> selecetedTraits = new ArrayList<String>();
//    public Map<String, TraitGuesser.TraitType> traitTypes = new HashMap<String, TraitGuesser.TraitType>();
   
    // Data 
    public List<PartitionData> dataPartitions = new ArrayList<PartitionData>();
    // ClockModel <=> TreeModel
    private List<PartitionClockModelTreeModelLink> partitionClockTreeLinks = new ArrayList<PartitionClockModelTreeModelLink>();
    
//    public boolean shareSameTreePrior = true;
    // list of starting tree from user import
    public List<Tree> userTrees = new ArrayList<Tree>();

//    public FixRateType rateOptionClockModel = FixRateType.FIX_FIRST_PARTITION; 
//    public boolean fixedSubstitutionRate = true;
//    public double meanSubstitutionRate = 1.0;
    public boolean unlinkPartitionRates = true;

    public Units.Type units = Units.Type.SUBSTITUTIONS;

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

    public String fileNameStem = MCMCPanel.fileNameStem;
    public String logFileName = null;
//    public boolean mapTreeLog = false;
//    public String mapTreeFileName = null;
    public List<String> treeFileName = new ArrayList<String>();
    public boolean substTreeLog = false;
    public List<String> substTreeFileName = new ArrayList<String>();
    public boolean operatorAnalysis = false;
    public String operatorAnalysisFileName = null;
    
//    public SubstitutionModelOptions substitutionModelOptions = new SubstitutionModelOptions(this);
    public ClockModelOptions clockModelOptions = new ClockModelOptions(this);
    public TreeModelOptions treeModelOptions = new TreeModelOptions(this);
    public PriorOptions priorOptions = new PriorOptions(this);

    public TraitOptions traitOptions = new TraitOptions(this);
    public STARBEASTOptions starBEASTOptions = new STARBEASTOptions(this);     
        
    public BEAUTiImporter beautiImporter = new BEAUTiImporter(this);    
    
    public BeautiTemplate beautiTemplate = new BeautiTemplate(this);

}