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
import dr.app.beauti.enumTypes.PriorScaleType;
import dr.app.beauti.enumTypes.TreePriorType;
import dr.app.beauti.mcmcpanel.MCMCPanel;
import dr.app.beauti.util.BeautiTemplate;
import dr.evolution.alignment.Alignment;
import dr.evolution.datatype.DataType;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxa;
import dr.evolution.util.Units;
import dr.evoxml.util.DateUnitsType;
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
        
        siteModelOptions = new SiteModelOptions(this);
        clockModelOptions = new ClockModelOptions(this);
        treeModelOptions = new TreeModelOptions(this);
        priorOptions = new PriorOptions(this);
        
//        traitsOptions = new TraitsOptions(this);
        starBEASTOptions = new STARBEASTOptions(this);

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

//        for (TraitData trait : getTraitsList()) { // all traits including locations
//            if (!trait.getName().equalsIgnoreCase(TraitData.Traits.TRAIT_SPECIES.toString()))
//        	   trait.gets.selectParameters(parameters);
//        }

        selectComponentParameters(this, parameters);
        
        selectComponentStatistics(this, parameters);
        
        priorOptions.selectParameters(parameters);

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

//        for (TraitData trait : getTraitsList()) { // all traits including locations
//        	if (!trait.getName().equalsIgnoreCase(TraitData.Traits.TRAIT_SPECIES.toString()))
//                trait.getTraitData().selectOperators(ops);
//        }

        selectComponentOperators(this, ops);

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
            if (model != null && (!activeModels.contains(model))
                    // species excluded
                    && (!partition.getName().equalsIgnoreCase(TraitData.Traits.TRAIT_SPECIES.toString()))) {
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

    public boolean validateDiffTaxa(List<PartitionData> partitionDataList) {
        Alignment ref = null;
        boolean legal = true;
        for (PartitionData partition : partitionDataList) {
            final Alignment a = partition.getAlignment();
            if (ref == null) {
                ref = a;
            } else {
                if (a.getTaxonCount() != ref.getTaxonCount()) {
                    legal = false;
                } else {
                    for (int k = 0; k < a.getTaxonCount(); ++k) {
                        if (ref.getTaxonIndex(a.getTaxonId(k)) == -1) {
                            legal = false;
                        }
                    }
                }
            }
        }        
        return legal;
    }

    // +++++++++++++ Traits +++++++++++++
    public static List<TraitData> getTraitsList() {
        List<TraitData> traits = new ArrayList<TraitData>();
        for (PartitionData partition : dataPartitions) {
            if (partition.getTraitType() != null) {
                traits.add((TraitData) partition);
            }
        }
        return traits;
    }

    public static List<TraitData> getDiscreteIntegerTraits() {
        List<TraitData> traits = new ArrayList<TraitData>();
        for (PartitionData partition : dataPartitions) {
            if (partition.getTraitType() != null && partition.getTraitType() != TraitData.TraitType.CONTINUOUS) {
                traits.add((TraitData) partition);
            }
        }
        return traits;
    }

    public static boolean hasDiscreteIntegerTraitsExcludeSpecies() { // exclude species at moment
        return getDiscreteIntegerTraits().size() > 1
                || (getDiscreteIntegerTraits().size() > 0 && (!containTrait(TraitData.Traits.TRAIT_SPECIES.toString())));
    }

    public static boolean containTrait(String traitName) {
        for (TraitData trait : getTraitsList()) {
            if (trait.getName().equalsIgnoreCase(traitName))
                return true;
        }
        return false;
    }

    public static int addTrait(TraitData newTrait) {
        int selRow;
        String traitName = newTrait.getName();
        if (containTrait(traitName)) {
            clearTraitValues(traitName); // Clear trait values
            selRow = dataPartitions.indexOf(getTrait(traitName));
            dataPartitions.set(selRow, newTrait);
        } else {
            dataPartitions.add(newTrait);
            selRow = getTraitsList().size() - 1; // start 0
        }
        return selRow;
    }

    public static void removeTrait(String traitName) {
        if (containTrait(traitName)) {
            clearTraitValues(traitName); // Clear trait values
            dataPartitions.remove(getTrait(traitName));
        }
    }

    public static void clearTraitValues(String traitName) {
        for (int i = 0; i < taxonList.getTaxonCount(); i++) {
            taxonList.getTaxon(i).setAttribute(traitName, "");
        }
    }

    public static TraitData getTrait(String traitName) {
        for (TraitData trait : getTraitsList()) {
            if (trait.getName().equalsIgnoreCase(traitName))
                return trait;
        }
        return null;
    }



    // ++++++++++++++++++++ message bar +++++++++++++++++
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

            if (TraitData.hasPhylogeographic()) {
                message += ";    Phylogeographic Analysis";
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

    public static Taxa taxonList = null;

    public List<Taxa> taxonSets = new ArrayList<Taxa>();
    public Map<Taxa, Boolean> taxonSetsMono = new HashMap<Taxa, Boolean>();

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
    public static final List<PartitionData> dataPartitions = new ArrayList<PartitionData>();
    
    // ClockModel <=> TreeModel
    private List<PartitionClockModelTreeModelLink> partitionClockTreeLinks = new ArrayList<PartitionClockModelTreeModelLink>();
    
    // list of starting tree from user import
    public List<Tree> userTrees = new ArrayList<Tree>();

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
    
    public SiteModelOptions siteModelOptions = new SiteModelOptions(this);
    public ClockModelOptions clockModelOptions = new ClockModelOptions(this);
    public TreeModelOptions treeModelOptions = new TreeModelOptions(this);
    public PriorOptions priorOptions = new PriorOptions(this);

//    public TraitsOptions traitsOptions = new TraitsOptions(this);
    public STARBEASTOptions starBEASTOptions = new STARBEASTOptions(this);

    public BeautiTemplate beautiTemplate = new BeautiTemplate(this);

    public static ArrayList<TraitData> getDiscreteTraitsExcludeSpecies() {
        return new ArrayList<TraitData>();  //Todo remove after
    }
}