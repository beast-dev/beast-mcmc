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
import dr.app.beauti.generator.Generator;
import dr.app.beauti.mcmcpanel.MCMCPanel;
import dr.app.beauti.priorsPanel.PriorType;
import dr.app.beauti.util.BeautiTemplate;
import dr.evolution.datatype.DataType;
import dr.evolution.tree.Tree;
import dr.evolution.util.Date;
import dr.evolution.util.Taxa;
import dr.evolution.util.Taxon;
import dr.evolution.util.Units;
import dr.evomodel.coalescent.GMRFFixedGridImportanceSampler;
import dr.evomodel.operators.TreeNodeSlide;
import dr.evomodel.speciation.SpeciesTreeModel;
import dr.evomodelxml.BirthDeathModelParser;
import dr.evomodelxml.YuleModelParser;
import dr.evoxml.TaxaParser;
import dr.evoxml.TaxonParser;
import dr.inference.operators.OperatorSchedule;
import dr.util.NumberFormatter;
import dr.xml.XMLParser;
import org.jdom.Document;
import org.jdom.Element;

import java.text.NumberFormat;
import java.text.ParseException;
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

        initGlobalParaAndOpers();
    }

    public BeautiOptions(ComponentFactory[] components) {
    	initGlobalParaAndOpers();

        // Install all the component's options from the given list of factories:
        for (ComponentFactory component : components) {
            addComponent(component.getOptions(this));
        }
    }   

    private void initGlobalParaAndOpers() {
//        double rateWeights = 3.0;
//
//        // A vector of relative rates across all partitions...
//        createParameter("allMus", "All the relative rates regarding codon positions");
//
//        // This only works if the partitions are of the same size...
////      createOperator("centeredMu", "Relative rates",
////              "Scales codon position rates relative to each other maintaining mean", "allMus",
////              OperatorType.CENTERED_SCALE, 0.75, rateWeights);
//        createOperator("deltaMu", RelativeRatesType.MU_RELATIVE_RATES.toString(),
//        		 "Currently use to scale codon position rates relative to each other maintaining mean", "allMus",
//                OperatorType.DELTA_EXCHANGE, 0.75, rateWeights);
        
        // only available for *BEAST and EBSP
        createUpDownAllOperator("upDownAllRatesHeights", "Up down all rates and heights", "Scales all rates inversely to node heights of the tree", 
        		demoTuning, branchWeights);        

    }

    public void initSpeciesParametersAndOperators() {
        double spWeights = 5.0;
        double spTuning = 0.9;

        createScaleParameter(TraitGuesser.Traits.TRAIT_SPECIES + "." + POP_MEAN, "Species tree: population hyper-parameter operator",
        		PriorScaleType.TIME_SCALE, 1.0, 0.0, Double.POSITIVE_INFINITY);
        // species tree Yule
        createScaleParameter(TraitGuesser.Traits.TRAIT_SPECIES + "." + YuleModelParser.YULE + "." + YuleModelParser.BIRTH_RATE,
                "Speices tree: Yule process birth rate", PriorScaleType.BIRTH_RATE_SCALE, 1.0, 0.0, Double.POSITIVE_INFINITY);

        // species tree Birth Death
        createScaleParameter(TraitGuesser.Traits.TRAIT_SPECIES + "." + BirthDeathModelParser.BIRTHDIFF_RATE_PARAM_NAME,
                "Speices tree: Birth Death model mean growth rate", PriorScaleType.BIRTH_RATE_SCALE, 1.0, 0.0, Double.POSITIVE_INFINITY);
        createParameter(TraitGuesser.Traits.TRAIT_SPECIES + "." + BirthDeathModelParser.RELATIVE_DEATH_RATE_PARAM_NAME,
                "Speices tree: Birth Death model relative death rate", PriorScaleType.BIRTH_RATE_SCALE, 0.5, 0.0, 1.0);

        createScaleParameter(SpeciesTreeModel.SPECIES_TREE + "." + Generator.SPLIT_POPS, "Species tree: population size operator",
        		PriorScaleType.TIME_SCALE, 1.0, 0.0, Double.POSITIVE_INFINITY);

        createParameter(TraitGuesser.Traits.TRAIT_SPECIES + "." + TreeNodeSlide.TREE_NODE_REHEIGHT, "Species tree: tree node operator");

        createScaleOperator(TraitGuesser.Traits.TRAIT_SPECIES + "." + POP_MEAN, spTuning, spWeights);

        createScaleOperator(TraitGuesser.Traits.TRAIT_SPECIES + "." + YuleModelParser.YULE + "." + YuleModelParser.BIRTH_RATE, demoTuning, demoWeights);
        createScaleOperator(TraitGuesser.Traits.TRAIT_SPECIES + "." + BirthDeathModelParser.BIRTHDIFF_RATE_PARAM_NAME, demoTuning, demoWeights);
        createScaleOperator(TraitGuesser.Traits.TRAIT_SPECIES + "." + BirthDeathModelParser.RELATIVE_DEATH_RATE_PARAM_NAME, demoTuning, demoWeights);

        createScaleOperator(SpeciesTreeModel.SPECIES_TREE + "." + Generator.SPLIT_POPS, 0.5, 94);

        createOperator(TraitGuesser.Traits.TRAIT_SPECIES + "." + TreeNodeSlide.TREE_NODE_REHEIGHT, OperatorType.NODE_REHIGHT, demoTuning, 94);
                
        //TODO: more
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
        
        meanDistance = 1.0;
        datesUnits = YEARS;
        datesDirection = FORWARDS;
        maximumTipHeight = 0.0;
        translation = 0;
        
        selecetedTraits.clear();
        traitTypes.clear();
        
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
        
        substitutionModelOptions = new SubstitutionModelOptions(this);
        clockModelOptions = new ClockModelOptions(this);
        
        beautiTemplate = new BeautiTemplate(this);
    }
    
    public void selectTaxonSetsStatistics(List<Parameter> params) {
    	
        if (taxonSets != null) {
            for (Taxa taxonSet : taxonSets) {
                Parameter statistic = statistics.get(taxonSet);
                if (statistic == null) {
                    statistic = new Parameter(taxonSet, "tMRCA for taxon set ");
                    statistics.put(taxonSet, statistic);
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
        
        for (PartitionSubstitutionModel model : getPartitionSubstitutionModels()) {
//          parameters.addAll(model.getParameters(multiplePartitions));
        	model.selectParameters(parameters);
        }
        
        for (PartitionClockModel model : getPartitionClockModels()) {
            model.selectParameters(parameters);
        }
        clockModelOptions.selectParameters(parameters);

        for (PartitionTreeModel tree : getPartitionTreeModels()) {
            tree.selectParameters(parameters);
        }

        for (PartitionTreePrior prior : getPartitionTreePriors()) {
            prior.selectParameters(parameters);
        }

        for (PartitionClockModelTreeModelLink clockTree : getPartitionClockTreeLinks()) {
        	clockTree.selectParameters(parameters);
            clockTree.selectStatistics(parameters);
        }
        
        if (isSpeciesAnalysis()) { // species
            selectParametersForSpecies(parameters);
        }

        selectComponentParameters(this, parameters);

        selectTaxonSetsStatistics(parameters);
        
        selectComponentStatistics(this, parameters);

//        boolean multiplePartitions = getTotalActivePartitionSubstitutionModelCount() > 1;
        // add all Parameter (with prefix) into parameters list     
        


        double growthRateMaximum = 1E6;
        double birthRateMaximum = 1E6;
        double substitutionRateMaximum = 100;
        double logStdevMaximum = 10;
        double substitutionParameterMaximum = 100;
        double initialRootHeight = 1;
        double initialRate = 1;


        if (clockModelOptions.getRateOptionClockModel() == FixRateType.FIX_MEAN) {
            double rate = clockModelOptions.getMeanRelativeRate();

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
//            if (dataReset) param.priorEdited = false;

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

//        dataReset = false;

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
        substitutionModelOptions.selectOperators(ops);
        
        for (PartitionClockModel model : getPartitionClockModels()) {
            model.selectOperators(ops);
        }
        clockModelOptions.selectOperators(ops);

        for (PartitionTreeModel tree : getPartitionTreeModels()) {
            tree.selectOperators(ops);
        }

        for (PartitionTreePrior prior : getPartitionTreePriors()) {
            prior.selectOperators(ops);
        }

        for (PartitionClockModelTreeModelLink clockTree : getPartitionClockTreeLinks()) {
            clockTree.selectOperators(ops);
        }

        if (isSpeciesAnalysis()) { // species
            selectOperatorsForSpecies(ops);
        }

        //up down all rates and trees operator only available for *BEAST and EBSP
        if (clockModelOptions.getRateOptionClockModel() == FixRateType.ESTIMATE && 
        		(isSpeciesAnalysis() || isEBSPSharingSamePrior())) {
        	ops.add(getOperator("upDownAllRatesHeights")); 
        }
        
        selectComponentOperators(this, ops);



        
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

        double initialRootHeight = 1;

        if (clockModelOptions.getRateOptionClockModel() == FixRateType.FIX_MEAN) {
            double rate = clockModelOptions.getMeanRelativeRate();

            if (hasData()) {
                initialRootHeight = meanDistance / rate;
                initialRootHeight = round(initialRootHeight, 2);
            }

        } else {
            if (maximumTipHeight > 0) {
                initialRootHeight = maximumTipHeight * 10.0;
            }
        }

        for (PartitionTreeModel tree : getPartitionTreeModels()) {
            Operator op = tree.getOperator("subtreeSlide");
            if (!op.tuningEdited) {
                op.tuning = initialRootHeight / 10.0;
            }
        }
        
//        for (Operator op : ops) {
//        	System.out.println(op.prefix + " + " + op.getName());
//        }

        return ops;
    }
    
    public boolean hasData() {
        return dataPartitions.size() > 0;
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
    
    private double round(double value, int sf) {
        NumberFormatter formatter = new NumberFormatter(sf);
        try {
            return NumberFormat.getInstance().parse(formatter.format(value)).doubleValue();
        } catch (ParseException e) {
            return value;
        }
    }
    
    // +++++++++++++++++++++++++++ *BEAST ++++++++++++++++++++++++++++++++++++
    private void selectParametersForSpecies(List<Parameter> params) {

        params.add(getParameter(TraitGuesser.Traits.TRAIT_SPECIES + "." + POP_MEAN));

        if (getPartitionTreePriors().get(0).getNodeHeightPrior() == TreePriorType.SPECIES_BIRTH_DEATH) {
            params.add(getParameter(TraitGuesser.Traits.TRAIT_SPECIES + "." + BirthDeathModelParser.BIRTHDIFF_RATE_PARAM_NAME));
            params.add(getParameter(TraitGuesser.Traits.TRAIT_SPECIES + "." + BirthDeathModelParser.RELATIVE_DEATH_RATE_PARAM_NAME));
        } else if (getPartitionTreePriors().get(0).getNodeHeightPrior() == TreePriorType.SPECIES_YULE) {
            params.add(getParameter(TraitGuesser.Traits.TRAIT_SPECIES + "." + YuleModelParser.YULE + "." + YuleModelParser.BIRTH_RATE));
        }

//    	params.add(getParameter(SpeciesTreeModel.SPECIES_TREE + "." + Generator.SPLIT_POPS));

        //TODO: more

    }
 
    private void selectOperatorsForSpecies(List<Operator> ops) {

        ops.add(getOperator(TraitGuesser.Traits.TRAIT_SPECIES + "." + POP_MEAN));

        if (getPartitionTreePriors().get(0).getNodeHeightPrior() == TreePriorType.SPECIES_BIRTH_DEATH) {
            ops.add(getOperator(TraitGuesser.Traits.TRAIT_SPECIES + "." + BirthDeathModelParser.BIRTHDIFF_RATE_PARAM_NAME));
            ops.add(getOperator(TraitGuesser.Traits.TRAIT_SPECIES + "." + BirthDeathModelParser.RELATIVE_DEATH_RATE_PARAM_NAME));
            
//            ops.add(getOperator("upDownBirthDeathSpeciesTree"));
//            ops.add(getOperator("upDownBirthDeathSTPop"));
//            
//            for (PartitionTreeModel tree : getPartitionTreeModels()) {
//            	ops.add(getOperator(tree.getPrefix() + "upDownBirthDeathGeneTree"));
//            }
        } else if (getPartitionTreePriors().get(0).getNodeHeightPrior() == TreePriorType.SPECIES_YULE) {
            ops.add(getOperator(TraitGuesser.Traits.TRAIT_SPECIES + "." + YuleModelParser.YULE + "." + YuleModelParser.BIRTH_RATE));
            
//            ops.add(getOperator("upDownYuleSpeciesTree"));
//            ops.add(getOperator("upDownYuleSTPop"));
//            
//            for (PartitionTreeModel tree : getPartitionTreeModels()) {
//            	ops.add(getOperator(tree.getPrefix() + "upDownYuleGeneTree"));
//            }
        }

        ops.add(getOperator(SpeciesTreeModel.SPECIES_TREE + "." + Generator.SPLIT_POPS));

        ops.add(getOperator(TraitGuesser.Traits.TRAIT_SPECIES + "." + TreeNodeSlide.TREE_NODE_REHEIGHT));
        //TODO: more
    }

    public boolean isSpeciesAnalysis() {
        return selecetedTraits.contains(TraitGuesser.Traits.TRAIT_SPECIES.toString());
    }

    public List<String> getSpeciesList() {
        List<String> species = new ArrayList<String>();
        String sp;

        if (taxonList != null) {
            for (int i = 0; i < taxonList.getTaxonCount(); i++) {
                Taxon taxon = taxonList.getTaxon(i);
                sp = taxon.getAttribute(TraitGuesser.Traits.TRAIT_SPECIES.toString()).toString();

                if (!species.contains(sp)) {
                    species.add(sp);
                }
            }
            return species;
        } else {
            return null;
        }
    } 

    public boolean isEBSPSharingSamePrior() {
    	if (getPartitionTreePriors().size() < 1) {
    		return false;
    	} else {
    		return (isShareSameTreePrior() && getPartitionTreePriors().get(0).getNodeHeightPrior() == TreePriorType.EXTENDED_SKYLINE);
    	}
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

    public List<PartitionSubstitutionModel> getPartitionSubstitutionModels(DataType dataType) {
        List<PartitionSubstitutionModel> models = new ArrayList<PartitionSubstitutionModel>();
        for (PartitionSubstitutionModel model : getPartitionSubstitutionModels()) {
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

    public int getTotalActivePartitionSubstitutionModelCount() {
        int totalPartitionCount = 0;
        for (PartitionSubstitutionModel model : getPartitionSubstitutionModels()) {
            totalPartitionCount += model.getCodonPartitionCount();
        }
        return totalPartitionCount;
    }

    /**
     * This returns an integer vector of the number of sites in each partition (including any codon partitions). These
     * are strictly in the same order as the 'mu' relative rates are listed.
     */
    public int[] getPartitionCodonWeights() {
        int[] weights = new int[getTotalActivePartitionSubstitutionModelCount()];

        int k = 0;
        for (PartitionSubstitutionModel model : getPartitionSubstitutionModels()) {
            for (PartitionData partition : dataPartitions) {
                if (partition.getPartitionSubstitutionModel() == model) {
                    model.addWeightsForPartition(partition, weights, k);
                }
            }
            k += model.getCodonPartitionCount();
        }

        assert (k == weights.length);

        return weights;
    }

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



//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // Data options
    public boolean allowDifferentTaxa = false;
    public DataType dataType = null;
//    public boolean dataReset = true;

    public Taxa taxonList = null;

    public List<Taxa> taxonSets = new ArrayList<Taxa>();
    public Map<Taxa, Boolean> taxonSetsMono = new HashMap<Taxa, Boolean>();

    public double meanDistance = 1.0;
    public int datesUnits = YEARS;
    public int datesDirection = FORWARDS;
    public double maximumTipHeight = 0.0;
    public int translation = 0;

    public DateGuesser dateGuesser = new DateGuesser();
    public TraitGuesser traitGuesser = new TraitGuesser();

    public List<String> selecetedTraits = new ArrayList<String>();
    public Map<String, TraitGuesser.TraitType> traitTypes = new HashMap<String, TraitGuesser.TraitType>();

    public final String SPECIES_TREE_FILE_NAME = TraitGuesser.Traits.TRAIT_SPECIES + "." + GMRFFixedGridImportanceSampler.TREE_FILE_NAME; // species.trees
    public final String POP_MEAN = "popMean";

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
    
    public SubstitutionModelOptions substitutionModelOptions = new SubstitutionModelOptions(this);
    public ClockModelOptions clockModelOptions = new ClockModelOptions(this);
    
    public BeautiTemplate beautiTemplate = new BeautiTemplate(this);
        
    @Override
    public String getPrefix() {
        // TODO Auto-generated method stub
        return null;
    }

}