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
        double rateWeights = 3.0;

        // A vector of relative rates across all partitions...
        createParameter("allMus", "All the relative rates");

        // This only works if the partitions are of the same size...
//      createOperator("centeredMu", "Relative rates",
//              "Scales codon position rates relative to each other maintaining mean", "allMus",
//              OperatorType.CENTERED_SCALE, 0.75, rateWeights);
        createOperator("deltaMu", "Relative rates",
//                "Changes partition relative rates relative to each other maintaining their mean", "allMus",
        		 "Currently use to scale codon position rates relative to each other maintaining mean", "allMus",
                OperatorType.DELTA_EXCHANGE, 0.75, rateWeights);
    }

    public void initSpeciesParametersAndOperators() {
        double spWeights = 5.0;
        double spTuning = 0.9;

        createScaleParameter(TraitGuesser.Traits.TRAIT_SPECIES + "." + POP_MEAN, "Species tree: population hyper-parameter operator",
                TIME_SCALE, 1.0, 0.0, Double.POSITIVE_INFINITY);
        // species tree Yule
        createScaleParameter(TraitGuesser.Traits.TRAIT_SPECIES + "." + YuleModelParser.YULE + "." + YuleModelParser.BIRTH_RATE,
                "Speices tree: Yule process birth rate", BIRTH_RATE_SCALE, 1.0, 0.0, Double.POSITIVE_INFINITY);

        // species tree Birth Death
        createScaleParameter(TraitGuesser.Traits.TRAIT_SPECIES + "." + BirthDeathModelParser.BIRTHDIFF_RATE_PARAM_NAME,
                "Speices tree: Birth Death Model BminusD rate", BIRTH_RATE_SCALE, 1.0, 0.0, Double.POSITIVE_INFINITY);
        createParameter(TraitGuesser.Traits.TRAIT_SPECIES + "." + BirthDeathModelParser.RELATIVE_DEATH_RATE_PARAM_NAME,
                "Speices tree: Birth Death Model DoverB rate", BIRTH_RATE_SCALE, 0.5, 0.0, 1.0);

        createScaleParameter(SpeciesTreeModel.SPECIES_TREE + "." + Generator.SPLIT_POPS, "Species tree: population size operator",
                TIME_SCALE, 1.0, 0.0, Double.POSITIVE_INFINITY);

        createParameter(TraitGuesser.Traits.TRAIT_SPECIES + "." + TreeNodeSlide.TREE_NODE_REHEIGHT, "Species tree: tree node operator");

        createScaleOperator(TraitGuesser.Traits.TRAIT_SPECIES + "." + POP_MEAN, spTuning, spWeights);

        createScaleOperator(TraitGuesser.Traits.TRAIT_SPECIES + "." + YuleModelParser.YULE + "." + YuleModelParser.BIRTH_RATE, demoTuning, demoWeights);
        createScaleOperator(TraitGuesser.Traits.TRAIT_SPECIES + "." + BirthDeathModelParser.BIRTHDIFF_RATE_PARAM_NAME, demoTuning, demoWeights);
        createScaleOperator(TraitGuesser.Traits.TRAIT_SPECIES + "." + BirthDeathModelParser.RELATIVE_DEATH_RATE_PARAM_NAME, demoTuning, demoWeights);

        createScaleOperator(SpeciesTreeModel.SPECIES_TREE + "." + Generator.SPLIT_POPS, 0.5, 94);

        createOperator(TraitGuesser.Traits.TRAIT_SPECIES + "." + TreeNodeSlide.TREE_NODE_REHEIGHT, OperatorType.NODE_REHIGHT, demoTuning, 94);
        
        // species tree Yule
//        createTagOperator("upDownYuleSpeciesTree", "Yule birth rate and species tree", "Scales Yule birth rate inversely to the species tree", 
//        		TraitGuesser.Traits.TRAIT_SPECIES + "." + YuleModelParser.YULE + "." + YuleModelParser.BIRTH_RATE,
//        		SpeciesTreeModel.SPECIES_TREE, Generator.SP_TREE, OperatorType.UP_DOWN, 0.75, branchWeights);
//        createOperator("upDownYuleSTPop", "Yule birth rate and species tree population size", 
//        		"Scales Yule birth rate inversely to the species tree population size", 
//                this.getParameter(TraitGuesser.Traits.TRAIT_SPECIES + "." + YuleModelParser.YULE + "." + YuleModelParser.BIRTH_RATE),
//                this.getParameter(SpeciesTreeModel.SPECIES_TREE + "." + Generator.SPLIT_POPS), OperatorType.UP_DOWN, 0.75, branchWeights);
//        
//        for (PartitionTreeModel tree : getPartitionTreeModels()) {
//	        createOperator(tree.getPrefix() + "upDownYuleGeneTree", "Species tree Yule and heights", 
//	        		"Scales Yule birth rate inversely to the gene tree", 
//	                this.getParameter(TraitGuesser.Traits.TRAIT_SPECIES + "." + YuleModelParser.YULE + "." + YuleModelParser.BIRTH_RATE),
//	                tree.getParameter("treeModel.allInternalNodeHeights"), OperatorType.UP_DOWN, 0.75, branchWeights);
//        }
//        // species tree Birth Death
//        createTagOperator("upDownBirthDeathSpeciesTree", "Birth death and species tree", 
//        		"Scales birth death BminusD rate inversely to the species tree", 
//        		TraitGuesser.Traits.TRAIT_SPECIES + "." + BirthDeathModelParser.BIRTHDIFF_RATE_PARAM_NAME,
//        		SpeciesTreeModel.SPECIES_TREE, Generator.SP_TREE, OperatorType.UP_DOWN, 0.75, branchWeights);
//        createOperator("upDownBirthDeathSTPop", "Yule birth rate and species tree population size", 
//        		"Scales birth death BminusD rate inversely to the species tree population size", 
//                this.getParameter(TraitGuesser.Traits.TRAIT_SPECIES + "." + BirthDeathModelParser.BIRTHDIFF_RATE_PARAM_NAME),
//                this.getParameter(SpeciesTreeModel.SPECIES_TREE + "." + Generator.SPLIT_POPS), OperatorType.UP_DOWN, 0.75, branchWeights);
//        
//        for (PartitionTreeModel tree : getPartitionTreeModels()) {
//	        createOperator(tree.getPrefix() + "upDownBirthDeathGeneTree", "Species tree Yule and heights", 
//	        		"Scales birth death BminusD rate inversely to the gene tree", 
//	                this.getParameter(TraitGuesser.Traits.TRAIT_SPECIES + "." + BirthDeathModelParser.BIRTHDIFF_RATE_PARAM_NAME),
//	                tree.getParameter("treeModel.allInternalNodeHeights"), OperatorType.UP_DOWN, 0.75, branchWeights);
//        }
        
        //TODO: more
    }

    /**
     * resets the options to the initial conditions
     */
    public void reset() {
        // Data options
        allowDifferentTaxa = false;
        dataType = null;
        dataReset = true;

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
        activedSameTreePrior = null;
        shareSameTreePrior = true;
        userTrees.clear();

        rateOptionClockModel = FixRateType.FIX_FIRST_PARTITION;
        meanSubstitutionRate = 1.0;
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
        
        for (PartitionClockModel model : getPartitionClockModels()) {
            model.selectParameters(parameters);
        }

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
        for (PartitionSubstitutionModel model : getPartitionSubstitutionModels()) {
//            parameters.addAll(model.getParameters(multiplePartitions));
        	parameters.addAll(model.getParameters());
        }


        double growthRateMaximum = 1E6;
        double birthRateMaximum = 1E6;
        double substitutionRateMaximum = 100;
        double logStdevMaximum = 10;
        double substitutionParameterMaximum = 100;
        double initialRootHeight = 1;
        double initialRate = 1;


        if (rateOptionClockModel == FixRateType.FIX_FIRST_PARTITION || rateOptionClockModel == FixRateType.FIX_MEAN) {
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

    /**
     * return an list of operators that are required
     *
     * @return the operator list
     */
    public List<Operator> selectOperators() {

        ArrayList<Operator> ops = new ArrayList<Operator>();

        for (PartitionClockModel model : getPartitionClockModels()) {
            model.selectOperators(ops);
        }

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

        selectComponentOperators(this, ops);

//        boolean multiplePartitions = getTotalActivePartitionSubstitutionModelCount() > 1;

        for (PartitionSubstitutionModel model : getPartitionSubstitutionModels()) {
            ops.addAll(model.getOperators());
        }

//        if (multiplePartitions) {
        if (hasCodon()) {
            Operator deltaMuOperator = getOperator("deltaMu");

            // update delta mu operator weight
            deltaMuOperator.weight = 0.0;
            for (PartitionSubstitutionModel pm : getPartitionSubstitutionModels()) {
                deltaMuOperator.weight += pm.getCodonPartitionCount();
            }

            ops.add(deltaMuOperator);
        }

        double initialRootHeight = 1;

        if (rateOptionClockModel == FixRateType.FIX_FIRST_PARTITION || rateOptionClockModel == FixRateType.FIX_MEAN) {
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
    
    /**
     * @return true either if the options have more than one partition or any partition is
     *         broken into codon positions.
     */
    public boolean hasCodon() {
//        final List<PartitionSubstitutionModel> models = options.getPartitionSubstitutionModels();
//        return (models.size() > 1 || models.get(0).getCodonPartitionCount() > 1);
    	for (PartitionSubstitutionModel model : getPartitionSubstitutionModels()) {
            if (model.getCodonPartitionCount() > 1) {
            	return true;
            }
        }
        return false;
    }
//    public boolean isFixedSubstitutionRate() {
//        return fixedSubstitutionRate;
//    }
    public void updateFixedRateClockModel() {
    	if (getPartitionClockModels().size() > 0) {
    		
	    	if (rateOptionClockModel == FixRateType.FIX_FIRST_PARTITION) {
	    		// fix rate of 1st partition
	    		for (PartitionClockModel model : getPartitionClockModels()) {
	    			if (getPartitionClockModels().indexOf(model) < 1) {
	    				model.setFixedRate(true);
	    			} else {
	    				model.setFixedRate(false);
	    			}
	            }
	    		
	    	} else if (rateOptionClockModel == FixRateType.FIX_MEAN) {
	    		// TODO check
	    		for (PartitionClockModel model : getPartitionClockModels()) {
	    			model.setFixedRate(true);
	            }
	    		
	    	} else {
	    		// estimate all rate
	    		for (PartitionClockModel model : getPartitionClockModels()) {
	    			model.setFixedRate(false);
	            }
	    		
	    	}
    	}
    }
    

    public double getMeanSubstitutionRate() {
        return meanSubstitutionRate;
    }

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

        if (activedSameTreePrior.getNodeHeightPrior() == TreePrior.SPECIES_BIRTH_DEATH) {
            params.add(getParameter(TraitGuesser.Traits.TRAIT_SPECIES + "." + BirthDeathModelParser.BIRTHDIFF_RATE_PARAM_NAME));
            params.add(getParameter(TraitGuesser.Traits.TRAIT_SPECIES + "." + BirthDeathModelParser.RELATIVE_DEATH_RATE_PARAM_NAME));
        } else if (activedSameTreePrior.getNodeHeightPrior() == TreePrior.SPECIES_YULE) {
            params.add(getParameter(TraitGuesser.Traits.TRAIT_SPECIES + "." + YuleModelParser.YULE + "." + YuleModelParser.BIRTH_RATE));
        }

//    	params.add(getParameter(SpeciesTreeModel.SPECIES_TREE + "." + Generator.SPLIT_POPS));

        //TODO: more

    }
 
    private void selectOperatorsForSpecies(List<Operator> ops) {

        ops.add(getOperator(TraitGuesser.Traits.TRAIT_SPECIES + "." + POP_MEAN));

        if (activedSameTreePrior.getNodeHeightPrior() == TreePrior.SPECIES_BIRTH_DEATH) {
            ops.add(getOperator(TraitGuesser.Traits.TRAIT_SPECIES + "." + BirthDeathModelParser.BIRTHDIFF_RATE_PARAM_NAME));
            ops.add(getOperator(TraitGuesser.Traits.TRAIT_SPECIES + "." + BirthDeathModelParser.RELATIVE_DEATH_RATE_PARAM_NAME));
            
//            ops.add(getOperator("upDownBirthDeathSpeciesTree"));
//            ops.add(getOperator("upDownBirthDeathSTPop"));
//            
//            for (PartitionTreeModel tree : getPartitionTreeModels()) {
//            	ops.add(getOperator(tree.getPrefix() + "upDownBirthDeathGeneTree"));
//            }
        } else if (activedSameTreePrior.getNodeHeightPrior() == TreePrior.SPECIES_YULE) {
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
    public int[] getPartitionWeights() {
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
        if (shareSameTreePrior) {
            if (activedSameTreePrior == null) {
                return activeTrees;
            } else {
                activeTrees.add(activedSameTreePrior);
            }
        } else {
            for (PartitionTreeModel model : getPartitionTreeModels()) {
                activeTrees.add(model.getPartitionTreePrior());
            }
        }

        return activeTrees;
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

    /**
     * return a list of parameters that are required
     *
     * @param params the parameter list
     */
//    private void selectParameters(List<Parameter> params) {
//
//        if (hasData()) {
//
//            // if not fixed then do mutation rate move and up/down move
//            boolean fixed = isFixedSubstitutionRate();
//            Parameter rateParam;
//
//            switch (clockType) {
//                case STRICT_CLOCK:
//                    rateParam = getParameter("clock.rate");
//                    if (!fixed) params.add(rateParam);
//                    break;
//
//                case UNCORRELATED_EXPONENTIAL:
//                    rateParam = getParameter(ClockType.UCED_MEAN);
//                    if (!fixed) params.add(rateParam);
//                    break;
//
//                case UNCORRELATED_LOGNORMAL:
//                    rateParam = getParameter(ClockType.UCLD_MEAN);
//                    if (!fixed) params.add(rateParam);
//                    params.add(getParameter(ClockType.UCLD_STDEV));
//                    break;
//
//                case AUTOCORRELATED_LOGNORMAL:
//                    rateParam = getParameter("treeModel.rootRate");
//                    if (!fixed) params.add(rateParam);
//                    params.add(getParameter("branchRates.var"));
//                    break;
//
//                case RANDOM_LOCAL_CLOCK:
//                    rateParam = getParameter("clock.rate");
//                    if (!fixed) params.add(rateParam);
//                    break;
//
//                default:
//                    throw new IllegalArgumentException("Unknown clock model");
//            }
//
//            /*if (clockType == ClockType.STRICT_CLOCK || clockType == ClockType.RANDOM_LOCAL_CLOCK) {
//				rateParam = getParameter("clock.rate");
//				if (!fixed) params.add(rateParam);
//			} else {
//				if (clockType == ClockType.UNCORRELATED_EXPONENTIAL) {
//					rateParam = getParameter("uced.mean");
//					if (!fixed) params.add(rateParam);
//				} else if (clockType == ClockType.UNCORRELATED_LOGNORMAL) {
//					rateParam = getParameter("ucld.mean");
//					if (!fixed) params.add(rateParam);
//					params.add(getParameter("ucld.stdev"));
//				} else {
//					throw new IllegalArgumentException("Unknown clock model");
//				}
//			}*/
//
//            rateParam.isFixed = fixed;
//
//        }
//
////        if (nodeHeightPrior == TreePrior.CONSTANT) {
////            params.add(getParameter("constant.popSize"));
////        } else if (nodeHeightPrior == TreePrior.EXPONENTIAL) {
////            params.add(getParameter("exponential.popSize"));
////            if (parameterization == GROWTH_RATE) {
////                params.add(getParameter("exponential.growthRate"));
////            } else {
////                params.add(getParameter("exponential.doublingTime"));
////            }
////        } else if (nodeHeightPrior == TreePrior.LOGISTIC) {
////            params.add(getParameter("logistic.popSize"));
////            if (parameterization == GROWTH_RATE) {
////                params.add(getParameter("logistic.growthRate"));
////            } else {
////                params.add(getParameter("logistic.doublingTime"));
////            }
////            params.add(getParameter("logistic.t50"));
////        } else if (nodeHeightPrior == TreePrior.EXPANSION) {
////            params.add(getParameter("expansion.popSize"));
////            if (parameterization == GROWTH_RATE) {
////                params.add(getParameter("expansion.growthRate"));
////            } else {
////                params.add(getParameter("expansion.doublingTime"));
////            }
////            params.add(getParameter("expansion.ancestralProportion"));
////        } else if (nodeHeightPrior == TreePrior.SKYLINE) {
////            params.add(getParameter("skyline.popSize"));
////        } else if (nodeHeightPrior == TreePrior.EXTENDED_SKYLINE) {
////            params.add(getParameter("demographic.populationSizeChanges"));
////            params.add(getParameter("demographic.populationMean"));
////        } else if (nodeHeightPrior == TreePrior.GMRF_SKYRIDE) {
//////            params.add(getParameter("skyride.popSize"));
////            params.add(getParameter("skyride.precision"));
////        } else if (nodeHeightPrior == TreePrior.YULE) {
////            params.add(getParameter("yule.birthRate"));
////        } else if (nodeHeightPrior == TreePrior.BIRTH_DEATH) {
////            params.add(getParameter(BirthDeathModelParser.BIRTHDIFF_RATE_PARAM_NAME));
////            params.add(getParameter(BirthDeathModelParser.RELATIVE_DEATH_RATE_PARAM_NAME));
////        }
////
////        params.add(getParameter("treeModel.rootHeight"));
//    }

//    private void selectStatistics(List<Parameter> params) {
//
//        if (taxonSets != null) {
//            for (Taxa taxonSet : taxonSets) {
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
//
//        if (clockType == ClockType.RANDOM_LOCAL_CLOCK) {
//            if (localClockRateChangesStatistic == null) {
//                localClockRateChangesStatistic = new Parameter("rateChanges", "number of random local clocks", true);
//                localClockRateChangesStatistic.priorType = PriorType.POISSON_PRIOR;
//                localClockRateChangesStatistic.poissonMean = 1.0;
//                localClockRateChangesStatistic.poissonOffset = 0.0;
//            }
//            if (localClockRatesStatistic == null) {
//                localClockRatesStatistic = new Parameter(ClockType.LOCAL_CLOCK + "." + "rates", "random local clock rates", false);
//
//                localClockRatesStatistic.priorType = PriorType.GAMMA_PRIOR;
//                localClockRatesStatistic.gammaAlpha = 0.5;
//                localClockRatesStatistic.gammaBeta = 2.0;
//            }
//            params.add(localClockRatesStatistic);
//            params.add(localClockRateChangesStatistic);
//        }
//
//        if (clockType != ClockType.STRICT_CLOCK) {
//            params.add(getParameter("meanRate"));
//            params.add(getParameter(RateStatistic.COEFFICIENT_OF_VARIATION));
//            params.add(getParameter("covariance"));
//        }
//
//    }

    /**
     * return a list of operators that are required
     *
     * @param ops the operator list
     */
//    private void selectOperators(List<Operator> ops) {
//
//        if (hasData()) {
//
//            if (!isFixedSubstitutionRate()) {
//                switch (clockType) {
//                    case STRICT_CLOCK:
//                        ops.add(getOperator("clock.rate"));
//                        ops.add(getOperator("upDownRateHeights"));
//                        break;
//
//                    case UNCORRELATED_EXPONENTIAL:
//                        ops.add(getOperator(ClockType.UCED_MEAN));
//                        ops.add(getOperator("upDownUCEDMeanHeights"));
//                        ops.add(getOperator("swapBranchRateCategories"));
//                        ops.add(getOperator("randomWalkBranchRateCategories"));
//                        ops.add(getOperator("unformBranchRateCategories"));
//                        break;
//
//                    case UNCORRELATED_LOGNORMAL:
//                        ops.add(getOperator(ClockType.UCLD_MEAN));
//                        ops.add(getOperator(ClockType.UCLD_STDEV));
//                        ops.add(getOperator("upDownUCLDMeanHeights"));
//                        ops.add(getOperator("swapBranchRateCategories"));
//                        ops.add(getOperator("randomWalkBranchRateCategories"));
//                        ops.add(getOperator("unformBranchRateCategories"));
//                        break;
//
//                    case AUTOCORRELATED_LOGNORMAL:
//                        ops.add(getOperator("scaleRootRate"));
//                        ops.add(getOperator("scaleOneRate"));
//                        ops.add(getOperator("scaleAllRates"));
//                        ops.add(getOperator("scaleAllRatesIndependently"));
//                        ops.add(getOperator("upDownAllRatesHeights"));
//                        ops.add(getOperator("branchRates.var"));
//                        break;
//
//                    case RANDOM_LOCAL_CLOCK:
//                        ops.add(getOperator("clock.rate"));
//                        ops.add(getOperator("upDownRateHeights"));
//                        ops.add(getOperator(ClockType.LOCAL_CLOCK + "." + "rates"));
//                        ops.add(getOperator(ClockType.LOCAL_CLOCK + "." + "changes"));
//                        ops.add(getOperator("treeBitMove"));
//                        break;
//
//                    default:
//                        throw new IllegalArgumentException("Unknown clock model");
//                }
//            } else {
//                switch (clockType) {
//                    case STRICT_CLOCK:
//                        // no parameter to operator on
//                        break;
//
//                    case UNCORRELATED_EXPONENTIAL:
//                        ops.add(getOperator("swapBranchRateCategories"));
//                        ops.add(getOperator("randomWalkBranchRateCategories"));
//                        ops.add(getOperator("unformBranchRateCategories"));
//                        break;
//
//                    case UNCORRELATED_LOGNORMAL:
//                        ops.add(getOperator(ClockType.UCLD_STDEV));
//                        ops.add(getOperator("swapBranchRateCategories"));
//                        ops.add(getOperator("randomWalkBranchRateCategories"));
//                        ops.add(getOperator("unformBranchRateCategories"));
//                        break;
//
//                    case AUTOCORRELATED_LOGNORMAL:
//                        ops.add(getOperator("scaleOneRate"));
//                        ops.add(getOperator("scaleAllRatesIndependently"));
//                        ops.add(getOperator("branchRates.var"));
//                        break;
//
//                    case RANDOM_LOCAL_CLOCK:
//                        ops.add(getOperator(ClockType.LOCAL_CLOCK + "." + "rates"));
//                        ops.add(getOperator(ClockType.LOCAL_CLOCK + "." + "changes"));
//                        ops.add(getOperator("treeBitMove"));
//                        break;
//
//                    default:
//                        throw new IllegalArgumentException("Unknown clock model");
//                }
//            }
//        }
//
////        if (nodeHeightPrior == TreePrior.CONSTANT) {
////            ops.add(getOperator("constant.popSize"));
////        } else if (nodeHeightPrior == TreePrior.EXPONENTIAL) {
////            ops.add(getOperator("exponential.popSize"));
////            if (parameterization == GROWTH_RATE) {
////                ops.add(getOperator("exponential.growthRate"));
////            } else {
////                ops.add(getOperator("exponential.doublingTime"));
////            }
////        } else if (nodeHeightPrior == TreePrior.LOGISTIC) {
////            ops.add(getOperator("logistic.popSize"));
////            if (parameterization == GROWTH_RATE) {
////                ops.add(getOperator("logistic.growthRate"));
////            } else {
////                ops.add(getOperator("logistic.doublingTime"));
////            }
////            ops.add(getOperator("logistic.t50"));
////        } else if (nodeHeightPrior == TreePrior.EXPANSION) {
////            ops.add(getOperator("expansion.popSize"));
////            if (parameterization == GROWTH_RATE) {
////                ops.add(getOperator("expansion.growthRate"));
////            } else {
////                ops.add(getOperator("expansion.doublingTime"));
////            }
////            ops.add(getOperator("expansion.ancestralProportion"));
////        } else if (nodeHeightPrior == TreePrior.SKYLINE) {
////            ops.add(getOperator("skyline.popSize"));
////            ops.add(getOperator("skyline.groupSize"));
////        } else if (nodeHeightPrior == TreePrior.GMRF_SKYRIDE) {
////            ops.add(getOperator("gmrfGibbsOperator"));
////        } else if (nodeHeightPrior == TreePrior.EXTENDED_SKYLINE) {
////            ops.add(getOperator("demographic.populationMean"));
////            ops.add(getOperator("demographic.popSize"));
////            ops.add(getOperator("demographic.indicators"));
////            ops.add(getOperator("demographic.scaleActive"));
////        } else if (nodeHeightPrior == TreePrior.YULE) {
////            ops.add(getOperator("yule.birthRate"));
////        } else if (nodeHeightPrior == TreePrior.BIRTH_DEATH) {
////            ops.add(getOperator(BirthDeathModelParser.BIRTHDIFF_RATE_PARAM_NAME));
////            ops.add(getOperator(BirthDeathModelParser.RELATIVE_DEATH_RATE_PARAM_NAME));
////        }
////
////        ops.add(getOperator("treeModel.rootHeight"));
////        ops.add(getOperator("uniformHeights"));
////
////        // if not a fixed tree then sample tree space
////        if (!fixedTree) {
////            ops.add(getOperator("subtreeSlide"));
////            ops.add(getOperator("narrowExchange"));
////            ops.add(getOperator("wideExchange"));
////            ops.add(getOperator("wilsonBalding"));
////        }
//
//    }

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
        //TODO:
//        dataElement.addContent(createChild("startingTreeType", startingTreeType.name()));

        dataElement.addContent(createChild("guessDates", guessDates));
        dataElement.addContent(createChild("guessType", dateGuesser.guessType.name()));
        dataElement.addContent(createChild("fromLast", dateGuesser.fromLast));
        dataElement.addContent(createChild("order", dateGuesser.order));
        dataElement.addContent(createChild("prefix", dateGuesser.prefix));
        dataElement.addContent(createChild("offset", dateGuesser.offset));
        dataElement.addContent(createChild("unlessLessThan", dateGuesser.unlessLessThan));
        dataElement.addContent(createChild("offset2", dateGuesser.offset2));

        root.addContent(dataElement);

        Element taxaElement = new Element(TaxaParser.TAXA);

        for (Taxa taxonSet : taxonSets) {
            Element taxonSetElement = new Element("taxonSet");
            taxonSetElement.addContent(createChild(XMLParser.ID, taxonSet.getId()));
            taxonSetElement.addContent(createChild("enforceMonophyly",
                    taxonSetsMono.get(taxonSet) ? "true" : "false"));
            for (int j = 0; j < taxonSet.getTaxonCount(); j++) {
                Element taxonElement = new Element(TaxonParser.TAXON);
                taxonElement.addContent(createChild(XMLParser.ID, taxonSet.getTaxon(j).getId()));
                taxonSetElement.addContent(taxonElement);
            }
            taxaElement.addContent(taxonSetElement);
        }

        root.addContent(taxaElement);

//        for (PartitionSubstitutionModel model : partitionModels) {
//
//            Element modelElement = new Element("model");
//
//            /*modelElement.addContent(createChild("nucSubstitutionModel", nucSubstitutionModel));
//                           modelElement.addContent(createChild("aaSubstitutionModel", aaSubstitutionModel));
//                           modelElement.addContent(createChild("binarySubstitutionModel", binarySubstitutionModel));
//                           modelElement.addContent(createChild("frequencyPolicy", frequencyPolicy));
//                           modelElement.addContent(createChild("gammaHetero", gammaHetero));
//                           modelElement.addContent(createChild("gammaCategories", gammaCategories));
//                           modelElement.addContent(createChild("invarHetero", invarHetero));
//                           modelElement.addContent(createChild("codonHeteroPattern", codonHeteroPattern));
//                           modelElement.addContent(createChild("maximumTipHeight", maximumTipHeight));
//                           modelElement.addContent(createChild("hasSetFixedSubstitutionRate", hasSetFixedSubstitutionRate));
//                           modelElement.addContent(createChild("meanSubstitutionRate", meanSubstitutionRate));
//                           modelElement.addContent(createChild("fixedSubstitutionRate", fixedSubstitutionRate));
//                           modelElement.addContent(createChild("unlinkedSubstitutionModel", unlinkedSubstitutionModel));
//                           modelElement.addContent(createChild("unlinkedHeterogeneityModel", unlinkedHeterogeneityModel));
//                           modelElement.addContent(createChild("unlinkedFrequencyModel", unlinkedFrequencyModel));
//                           modelElement.addContent(createChild("clockModel", clockModel));
//                           modelElement.addContent(createChild("nodeHeightPrior", nodeHeightPrior));
//                           modelElement.addContent(createChild("parameterization", parameterization));
//                           modelElement.addContent(createChild("skylineGroupCount", skylineGroupCount));
//                           modelElement.addContent(createChild("skylineModel", skylineModel));
//                           modelElement.addContent(createChild("fixedTree", fixedTree)); */
//
//            root.addContent(modelElement);
//        }

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

        Element taxaElement = root.getChild(TaxaParser.TAXA);
        Element modelElement = root.getChild("model");
        Element priorsElement = root.getChild("priors");
        Element operatorsElement = root.getChild("operators");
        Element mcmcElement = root.getChild("mcmc");
        /*
                  if (taxaElement != null) {
                      for (Object ts : taxaElement.getChildren("taxonSet")) {
                          Element taxonSetElement = (Element) ts;

                          String id = getStringChild(taxonSetElement, XMLParser.ID, "");
                          final Taxa taxonSet = new Taxa(id);

                          Boolean enforceMonophyly = Boolean.valueOf(getStringChild(taxonSetElement, "enforceMonophyly", "false"));
                          for (Object o : taxonSetElement.getChildren("taxon")) {
                              Element taxonElement = (Element) o;
                              String taxonId = getStringChild(taxonElement, XMLParser.ID, "");
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

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // Data options
    public boolean allowDifferentTaxa = false;
    public DataType dataType = null;
    public boolean dataReset = true;

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
    // Substitution Model
//    List<PartitionSubstitutionModel> partitionModels = new ArrayList<PartitionSubstitutionModel>();
    // Clock Model
    //    List<PartitionClockModel> clockModels = new ArrayList<PartitionClockModel>();
    // Tree
    //    List<PartitionTreeModel> partitionTreeModels = new ArrayList<PartitionTreeModel>();
    // PopSize
    public PartitionTreePrior activedSameTreePrior = null;
    //    List<PartitionTreePrior> partitionTreePriors = new ArrayList<PartitionTreePrior>();
    public boolean shareSameTreePrior = true;
    // list of starting tree from user import
    public List<Tree> userTrees = new ArrayList<Tree>();

    public FixRateType rateOptionClockModel = FixRateType.FIX_FIRST_PARTITION; 
//    public boolean fixedSubstitutionRate = true;
    public double meanSubstitutionRate = 1.0;
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
    
    
    @Override
    public String getPrefix() {
        // TODO Auto-generated method stub
        return null;
    }

}