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
import dr.app.beauti.types.PriorScaleType;
import dr.app.beauti.types.TreePriorType;
import dr.app.beauti.util.BeautiTemplate;
import dr.evolution.alignment.Alignment;
import dr.evolution.datatype.DataType;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxa;
import dr.evolution.util.TaxonList;
import dr.evolution.util.Units;
import dr.evoxml.util.DateUnitsType;
import dr.inference.operators.OperatorSchedule;

import java.util.*;

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
//        dataType = null;
//        dataReset = true;

        taxonList = null;
        taxonSets.clear();
        taxonSetsMono.clear();
        taxonSetsIncludeStem.clear();
        taxonSetsTreeModel.clear();

//        meanDistance = 1.0;
        datesUnits = DateUnitsType.YEARS;
        datesDirection = DateUnitsType.FORWARDS;
        maximumTipHeight = 0.0;
        translation = 0;

//        selecetedTraits.clear();
//        traitTypes.clear();

        dataPartitions.clear();
        traits.clear();
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
        allowOverwriteLog = false;
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
        useStarBEAST = false;
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
                    PartitionTreeModel treeModel = taxonSetsTreeModel.get(taxa);

                    statistic = new Parameter.Builder(taxa.getId(), "")
                            .taxaId(treeModel.getPrefix() + taxa.getId())
                            .isStatistic(true).isNodeHeight(true).scaleType(PriorScaleType.NONE)
                            .initial(Double.NaN).lower(0.0).upper(Double.POSITIVE_INFINITY).build();
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

        for (PartitionClockModelSubstModelLink clockSubst : getTraitClockSubstLinks()) {
            clockSubst.selectParameters(parameters);
        }

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

        if (useStarBEAST) { // species
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

        for (PartitionClockModelSubstModelLink clockSubst : getTraitClockSubstLinks()) {
            clockSubst.selectOperators(ops);
        }

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

        if (useStarBEAST) { // species
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

    public boolean contains(DataType dataType) {
        for (AbstractPartitionData pd : dataPartitions) {
            if (pd.getDataType().getType() == dataType.getType()) {
                return true;
            }
        }
        return false;
    }

    public boolean hasPartitionData(String name) {
        for (AbstractPartitionData pd : dataPartitions) {
            if (name.equalsIgnoreCase(pd.getName())) {
                return true;
            }
        }
        return false;
    }

    public PartitionData getPartitionData(Alignment alignment) {
        for (AbstractPartitionData partition : dataPartitions) {
            if (partition instanceof PartitionData) {
                if (((PartitionData) partition).getAlignment() == alignment)
                    return (PartitionData) partition;
            }
        }
        return null;
    }

    public List<AbstractPartitionData> getAllPartitionData(TraitData trait) {
        List<AbstractPartitionData> pdList = new ArrayList<AbstractPartitionData>();
        for (AbstractPartitionData pd : dataPartitions) {
            if (pd.getTrait() == trait) {
                pdList.add(pd);
            }
        }
        return pdList;
    }

    public List<AbstractPartitionData> getAllPartitionData(PartitionOptions model) {
        if (model instanceof PartitionSubstitutionModel) {
            return getAllPartitionData((PartitionSubstitutionModel) model);
        } else if (model instanceof PartitionClockModel) {
            return getAllPartitionData((PartitionClockModel) model);
        } else if (model instanceof PartitionTreeModel) {
            return getAllPartitionData((PartitionTreeModel) model);
        } else {
            return null;
        }
    }

    public List<AbstractPartitionData> getAllPartitionData(PartitionSubstitutionModel model) {
        List<AbstractPartitionData> pdList = new ArrayList<AbstractPartitionData>();
        for (AbstractPartitionData pd : dataPartitions) {
            if (pd.getPartitionSubstitutionModel() == model) {
                pdList.add(pd);
            }
        }
        return pdList;
    }

    public List<AbstractPartitionData> getAllPartitionData(PartitionClockModel model) {
        List<AbstractPartitionData> pdList = new ArrayList<AbstractPartitionData>();
        for (AbstractPartitionData pd : dataPartitions) {
            if (pd.getPartitionClockModel() == model) {
                pdList.add(pd);
            }
        }
        return pdList;
    }

    public List<AbstractPartitionData> getAllPartitionData(PartitionTreeModel model) {
        List<AbstractPartitionData> pdList = new ArrayList<AbstractPartitionData>();
        for (AbstractPartitionData pd : dataPartitions) {
            if (pd.getPartitionTreeModel() == model) {
                pdList.add(pd);
            }
        }
        return pdList;
    }

    public List<AbstractPartitionData> getAllPartitionData(ClockModelGroup clockModelGroup) {
        List<AbstractPartitionData> pdList = new ArrayList<AbstractPartitionData>();
        for (AbstractPartitionData pd : dataPartitions) {
            if (pd.getPartitionClockModel().getClockModelGroup() == clockModelGroup) {
                pdList.add(pd);
            }
        }
        return pdList;
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

    public List<PartitionSubstitutionModel> getPartitionSubstitutionModels(DataType dataType) {
        List<PartitionSubstitutionModel> models = new ArrayList<PartitionSubstitutionModel>();
        for (PartitionSubstitutionModel model : getPartitionSubstitutionModels(dataPartitions)) {
            if (model.getDataType() == dataType) {
                models.add(model);
            }
        }
        return models;
    }

    public List<PartitionSubstitutionModel> getPartitionSubstitutionModels(List<? extends AbstractPartitionData> givenDataPartitions) {

        List<PartitionSubstitutionModel> activeModels = new ArrayList<PartitionSubstitutionModel>();

        for (AbstractPartitionData partition : givenDataPartitions) {
            PartitionSubstitutionModel model = partition.getPartitionSubstitutionModel();
            if (model != null && (!activeModels.contains(model))
                    // species excluded
                    && (!partition.getName().equalsIgnoreCase(TraitData.TRAIT_SPECIES))) {
                activeModels.add(model);
            }
        }

        return activeModels;
    }

    public List<PartitionSubstitutionModel> getPartitionSubstitutionModels() {
        return getPartitionSubstitutionModels(dataPartitions);
    }

//    public List<PartitionSubstitutionModel> getPartitionTraitsSubstitutionModels() {
//        return getPartitionSubstitutionModels(getTraitsList());
//    }
//
//    public List<PartitionSubstitutionModel> getPartitionNonTraitsSubstitutionModels() {
//        return getPartitionSubstitutionModels(getNonTraitsDataList());
//    }


    // ++++++++++++++ Partition Clock Model ++++++++++++++

    public List<PartitionClockModel> getPartitionClockModels(List<? extends AbstractPartitionData> givenDataPartitions) {

        List<PartitionClockModel> activeModels = new ArrayList<PartitionClockModel>();

        for (AbstractPartitionData partition : givenDataPartitions) {
            PartitionClockModel model = partition.getPartitionClockModel();
            if (model != null && (!activeModels.contains(model))
                    // species excluded
                    && (!partition.getName().equalsIgnoreCase(TraitData.TRAIT_SPECIES))) {
                activeModels.add(model);
            }
        }

        return activeModels;
    }

    public List<PartitionClockModel> getPartitionClockModels(DataType dataType) {
        List<PartitionClockModel> models = new ArrayList<PartitionClockModel>();
        for (PartitionClockModel model : getPartitionClockModels()) {
            if (model.getDataType().getType() == dataType.getType()) {
                models.add(model);
            }
        }
        return models;
    }

    public List<PartitionClockModel> getPartitionClockModels(ClockModelGroup group) {
        List<PartitionClockModel> models = new ArrayList<PartitionClockModel>();
        for (PartitionClockModel model : getPartitionClockModels()) {
            if (model.getClockModelGroup() == group) {
                models.add(model);
            }
        }
        return models;
    }

//    public List<PartitionClockModel> getPartitionNonTraitsClockModels() {
//        return getPartitionClockModels(getNonTraitsDataList());
//    }
//
//    public List<PartitionClockModel> getPartitionTraitsClockModels() {
//        return getPartitionClockModels(getTraitsList());
//    }

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

    public List<PartitionTreeModel> getPartitionTreeModels(List<? extends AbstractPartitionData> givenDataPartitions) {

        List<PartitionTreeModel> activeTrees = new ArrayList<PartitionTreeModel>();

        for (AbstractPartitionData partition : givenDataPartitions) {
            PartitionTreeModel tree = partition.getPartitionTreeModel();
            if (tree != null && (!activeTrees.contains(tree))
                    // species excluded
                    && (!partition.getName().equalsIgnoreCase(TraitData.TRAIT_SPECIES))) {
                activeTrees.add(tree);
            }
        }

        return activeTrees;
    }

    public List<PartitionTreeModel> getPartitionTreeModels() {
        return getPartitionTreeModels(dataPartitions);
    }

//    public List<PartitionTreeModel> getNonTraitPartitionTreeModels() {
//        return getPartitionTreeModels(getNonTraitsDataList());
//    }

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

    public List<PartitionClockModelSubstModelLink> getTraitClockSubstLinks() {
        return partitionClockSubstLinks;
    }

    public PartitionClockModelTreeModelLink getPartitionClockTreeLink(PartitionClockModel model, PartitionTreeModel tree) {
        for (PartitionClockModelTreeModelLink clockTree : getPartitionClockTreeLinks()) {
            if (clockTree.getPartitionClockModel().equals(model) && clockTree.getPartitionTreeTree().equals(tree)) {
                return clockTree;
            }
        }

        return null;
    }

    public void updatePartitionAllLinks() {
        partitionClockTreeLinks.clear();
        partitionClockSubstLinks.clear();

        for (PartitionClockModel model : getPartitionClockModels()) {
            for (PartitionTreeModel tree : getPartitionTreeModels(getAllPartitionData(model))) {
                PartitionClockModelTreeModelLink clockTree = new PartitionClockModelTreeModelLink(this, model, tree);

                if (!partitionClockTreeLinks.contains(clockTree)) {
                    partitionClockTreeLinks.add(clockTree);
                }
            }
        }

    }

    public void updateAll() {
        updatePartitionAllLinks();
        for (ClockModelGroup clockModelGroup : clockModelOptions.getClockModelGroups()) {
            clockModelOptions.fixRateOfFirstClockPartition(clockModelGroup);
        }
    }

    // update links (e.g List<PartitionData> allPartitionData), after use (e.g partition.setPartitionSubstitutionModel(model))

//    public void updateLinksBetweenPDPCMPSMPTMPTPP() {
//        for (PartitionSubstitutionModel model : getPartitionSubstitutionModels()) {
//            model.clearAllPartitionData();
//        }
//
//        for (PartitionClockModel model : getPartitionClockModels()) {
//            model.clearAllPartitionData();
//        }
//
//        for (PartitionTreeModel tree : getPartitionTreeModels()) {
//            tree.clearAllPartitionData();
//        }
//
//        //TODO update PartitionTreePrior ?
//
//        for (PartitionData partition : dataPartitions) {
//            PartitionSubstitutionModel psm = partition.getPartitionSubstitutionModel();
//            if (!psm.getAllPartitionData().contains(partition)) {
//                psm.addPartitionData(partition);
//            }
//
//            PartitionClockModel pcm = partition.getPartitionClockModel();
//            if (!pcm.getAllPartitionData().contains(partition)) {
//                pcm.addPartitionData(partition);
//            }
//
//            PartitionTreeModel ptm = partition.getPartitionTreeModel();
//            if (!ptm.getAllPartitionData().contains(partition)) {
//                ptm.addPartitionData(partition);
//            }
//        }
//
//    }

    public double getAveWeightedMeanDistance(List<AbstractPartitionData> partitions) {
        double meanDistance = 0;
        double totalSite = 0;
        for (AbstractPartitionData partition : partitions) {
            meanDistance = meanDistance + partition.getMeanDistance() * partition.getSiteCount();
            totalSite = totalSite + partition.getSiteCount();
        }

        if (totalSite == 0) {
            return 0;
        } else {
            return meanDistance / totalSite;
        }
    }

    public boolean validateDiffTaxa(List<AbstractPartitionData> partitionDataList) {
        TaxonList ref = null;
        boolean legal = true;
        for (AbstractPartitionData partition : partitionDataList) {
            final TaxonList a = partition.getTaxonList();
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

//    public List<PartitionData> getNonTraitsDataList() {
//        List<PartitionData> nonTraitsData = new ArrayList<PartitionData>();
//        for (PartitionData partition : dataPartitions) {
//            if (partition.getTraitType() == null) {
//                nonTraitsData.add(partition);
//            }
//        }
//        return nonTraitsData;
//    }
//
//    public List<TraitData> getTraitsList() {
//        List<TraitData> traits = new ArrayList<TraitData>();
//        for (PartitionData partition : dataPartitions) {
//            if (partition.getTraitType() != null) {
//                traits.add((TraitData) partition);
//            }
//        }
//        return traits;
//    }
//
//    public List<TraitData> getDiscreteIntegerTraits() {
//        List<TraitData> traits = new ArrayList<TraitData>();
//        for (PartitionData partition : dataPartitions) {
//            if (partition.getTraitType() != null && partition.getTraitType() != TraitData.TraitType.CONTINUOUS) {
//                traits.add((TraitData) partition);
//            }
//        }
//        return traits;
//    }
//
//    public boolean hasDiscreteIntegerTraitsExcludeSpecies() { // exclude species at moment
//        return getDiscreteIntegerTraits().size() > 1
//                || (getDiscreteIntegerTraits().size() > 0 && (!traitExists(TraitData.TRAIT_SPECIES)));
//    }

    public int getIndexOfTrait(String traitName) {
        int i = 0;
        for (TraitData trait : traits) {
            if (trait.getName().equalsIgnoreCase(traitName)) {
                return i;
            }
            i++;
        }
        return -1;
    }

    public boolean traitExists(String traitName) {
        return getIndexOfTrait(traitName) != -1;
    }

    public int addTrait(TraitData newTrait) {
        int selRow = getIndexOfTrait(newTrait.getName());

        if (selRow == -1) {
            traits.add(newTrait);
            selRow = traits.size() - 1; // start 0
        }

        return selRow; // only for trait panel
    }

    public int createPartitionForTrait(String name, TraitData trait) {
        int selRow = -1;

        PartitionData partition = new PartitionData(this, name, trait);
        dataPartitions.add(partition);
        selRow = dataPartitions.size() - 1;

        if (partition.getPartitionSubstitutionModel() == null) {
            PartitionSubstitutionModel substModel = new PartitionSubstitutionModel(this, partition);
            partition.setPartitionSubstitutionModel(substModel);
        }

        if (partition.getPartitionClockModel() == null) {
            // PartitionClockModel based on PartitionData
            PartitionClockModel pcm = new PartitionClockModel(this, partition);
            partition.setPartitionClockModel(pcm);

            clockModelOptions.addClockModelGroup(pcm);
        }

        if (partition.getPartitionTreeModel() == null) {
            partition.setPartitionTreeModel(getPartitionTreeModels().get(0));// always use 1st tree
//            getPartitionTreeModels().get(0).addPartitionData(newTrait);
        }

        return selRow; // only for trait panel
    }

    public void removeTrait(String traitName) {
        if (traitExists(traitName)) {
            clearTraitValues(traitName); // Clear trait values
            dataPartitions.remove(getTrait(traitName));
        }
    }

    public void clearTraitValues(String traitName) {
        for (int i = 0; i < taxonList.getTaxonCount(); i++) {
            taxonList.getTaxon(i).setAttribute(traitName, "");
        }
    }

    public TraitData getTrait(String traitName) {
        for (TraitData trait : traits) {
            if (trait.getName().equalsIgnoreCase(traitName))
                return trait;
        }
        return null;
    }

//    public boolean hasAlignmentPartition() {
//        for (AbstractPartitionData partition : dataPartitions) {
//            if (partition instanceof PartitionData) {
//                if (((PartitionData) partition).getAlignment() != null) {
//                    return true;
//                }
//            }
//        }
//        return false;
//    }


    public boolean hasDiscreteTrait() {
        for (TraitData traitData : traits) {
            if (traitData.getTraitType() == TraitData.TraitType.DISCRETE) {
                return true;
            }
        }
        return false;
    }

    public boolean hasDiscreteTraitPartition() {
        for (AbstractPartitionData partition : dataPartitions) {
            if (partition.getTrait() != null && partition.getTrait().getTraitType() == TraitData.TraitType.DISCRETE) {
                return true;
            }
        }
        return false;
    }

    public boolean hasContinuousTrait() {
        for (TraitData traitData : traits) {
            if (traitData.getTraitType() == TraitData.TraitType.CONTINUOUS) {
                return true;
            }
        }
        return false;
    }

    public boolean hasContinuousTraitPartition() {
        for (AbstractPartitionData partition : dataPartitions) {
            if (partition.getTrait() != null && partition.getTrait().getTraitType() == TraitData.TraitType.CONTINUOUS) {
                return true;
            }
        }
        return false;
    }

    public Set<String> getStatesForDiscreteModel(PartitionSubstitutionModel model) {
        Set<String> states = new TreeSet<String>();
        for (AbstractPartitionData partition : getAllPartitionData(model)) {
            Set<String> newStates = partition.getTrait().getStatesOfTrait(taxonList);

            if (states.size() > 0) {
                Set<String> shared = new HashSet<String>(states);
                shared.retainAll(newStates);
                if (shared.size() == 0) {
                    throw new IllegalArgumentException("For discrete trait partitions to have a linked model they must share states");
                }
            }

            states.addAll(newStates);
        }

        if (states.size() < 1) throw new IllegalArgumentException("The number of states must be greater than 1");

        return states;
    }

    public static TraitData.TraitType guessTraitType(TaxonList taxa, String name) {
        TraitData.TraitType type = TraitData.TraitType.DISCRETE;
        return type;
    }

    // ++++++++++++++++++++ message bar +++++++++++++++++

    public String statusMessage() {
//        String message = "<html><p>";
        String message = "";
        if (hasData()) {
            if (allowDifferentTaxa) {
                message += "Data contains different taxa: " + taxonList.getTaxonCount() + " taxa in total, ";
            } else {
                message += "Data: " + taxonList.getTaxonCount() + " taxa, ";
            }

            message += dataPartitions.size() + (dataPartitions.size() > 1 ? " partitions" : " partition");

            if (starBEASTOptions.getSpeciesList() != null && useStarBEAST) {
                int num = starBEASTOptions.getSpeciesList().size();
                message += ", " + num + " species"; // species is both singular and plural
            }

            if (userTrees.size() > 0) {
                message += ", " + userTrees.size() + " user" +
                        (userTrees.size() > 1 ? " trees" : " tree");
            }

            if (useStarBEAST) {
                message += "; Species Tree Ancestral Reconstruction (*BEAST)";
            }

//            if (hasPhylogeographic()) {
//                message += ";    Phylogeographic Analysis";
//            }

            message += "; " + clockModelOptions.statusMessageClockModel();

        } else if (userTrees.size() > 0) { // TODO
            message += "Trees only : " + userTrees.size() +
                    (userTrees.size() > 1 ? " trees, " : " tree, ") +
                    taxonList.getTaxonCount() + " taxa";
        } else if (taxonList != null && taxonList.getTaxonCount() > 0) {  // TODO
            message += "Taxa only: " + taxonList.getTaxonCount() + " taxa";
        } else {
            message += "No data loaded";
        }
//        message += "</p></html>";
        return message;
    }

    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // Data options
    public boolean allowDifferentTaxa = false;
//    public DataType dataType = null;

    public Taxa taxonList = null; // union set of all taxa in all partitions. todo change to List<Taxa> regarding data type?

    public List<Taxa> taxonSets = new ArrayList<Taxa>();
    public Map<Taxa, Boolean> taxonSetsMono = new HashMap<Taxa, Boolean>();
    public Map<Taxa, Boolean> taxonSetsIncludeStem = new HashMap<Taxa, Boolean>();
    public Map<Taxa, PartitionTreeModel> taxonSetsTreeModel = new HashMap<Taxa, PartitionTreeModel>();

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
    public List<AbstractPartitionData> dataPartitions = new ArrayList<AbstractPartitionData>();
    public List<TraitData> traits = new ArrayList<TraitData>();

    // ClockModel <=> TreeModel
    private List<PartitionClockModelTreeModelLink> partitionClockTreeLinks = new ArrayList<PartitionClockModelTreeModelLink>();
    // ClockModel <=> SubstModel
    private List<PartitionClockModelSubstModelLink> partitionClockSubstLinks = new ArrayList<PartitionClockModelSubstModelLink>();

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
    public boolean allowOverwriteLog = false;
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
    public boolean useStarBEAST = false;
    public STARBEASTOptions starBEASTOptions = new STARBEASTOptions(this);

    public BeautiTemplate beautiTemplate = new BeautiTemplate(this);

//    public static ArrayList<TraitData> getDiscreteTraitsExcludeSpecies() {
//        return new ArrayList<TraitData>();  //Todo remove after
//    }

}