/*
 * BeastGenerator.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
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
 *  BEAST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 *
 */

package dr.app.treelime.generator;

import dr.app.beast.BeastVersion;
import dr.app.treelime.TreeLimeApp;
import dr.app.treelime.TreeLimeFrame;
import dr.app.treelime.components.ComponentFactory;
import dr.app.treelime.options.*;
import dr.app.treelime.types.*;
import dr.app.beauti.util.XMLWriter;
import dr.app.util.Arguments;
import dr.evolution.alignment.Alignment;
import dr.evolution.datatype.DataType;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxa;
import dr.evolution.util.Taxon;
import dr.evolution.util.TaxonList;
import dr.evolution.util.Units;
import dr.evoxml.AlignmentParser;
import dr.evoxml.DateParser;
import dr.evoxml.TaxaParser;
import dr.evoxml.TaxonParser;
import dr.evomodel.tree.DefaultTreeModel;
import dr.evomodelxml.tree.CTMCScalePriorParser;
import dr.evomodelxml.tree.TimeTreeExportParser;
import dr.evomodelxml.tree.TreeLoggerParser;
import dr.inferencexml.model.CompoundLikelihoodParser;
import dr.inferencexml.model.CompoundParameterParser;
import dr.inference.model.ParameterParser;
import dr.util.TransformParsers;
import dr.xml.Report;
import dr.xml.XMLParseException;
import dr.evomodelxml.continuous.hmc.MaximizeWrtParameterParser;
import dr.inferencexml.operators.SimpleOperatorScheduleParser;
import dr.util.Attribute;
import dr.util.Pair;
import dr.util.Version;
import dr.xml.AttributeParser;
import dr.xml.XMLParser;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;

/**
 * This class holds all the data for the current BEAUti Document
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @author Walter Xie
 */
public class BeastGenerator extends Generator {

    private final static Version VERSION = new BeastVersion();
    private static final String MESSAGE_CAL_YULE = "Calibrated Yule requires 1 calibrated internal node \n" +
            "with a proper prior and monophyly enforced for each tree.";
    private final String MESSAGE_CAL = "\nas another element (taxon, sequence, taxon set, species, etc.):\nAll ids should be unique.";

    private final AlignmentGenerator alignmentGenerator;
    private final PatternListGenerator patternListGenerator;
    private final TreePriorGenerator treePriorGenerator;
    private final TreeLikelihoodGenerator treeLikelihoodGenerator;
    private final SubstitutionModelGenerator substitutionModelGenerator;
    private final InitialTreeGenerator initialTreeGenerator;
    private final TreeModelGenerator treeModelGenerator;
    private final ClockModelGenerator clockModelGenerator;
    private final OperatorsGenerator operatorsGenerator;
    private final ParameterPriorGenerator parameterPriorGenerator;
    private final LogGenerator logGenerator;
    private final TMRCAStatisticsGenerator tmrcaStatisticsGenerator;

    public BeastGenerator(TreeLimeOptions options, ComponentFactory[] components) {
        super(options, components);

        alignmentGenerator = new AlignmentGenerator(options, components);
        patternListGenerator = new PatternListGenerator(options, components);
        tmrcaStatisticsGenerator = new TMRCAStatisticsGenerator(options, components);
        substitutionModelGenerator = new SubstitutionModelGenerator(options, components);
        treePriorGenerator = new TreePriorGenerator(options, components);
        treeLikelihoodGenerator = new TreeLikelihoodGenerator(options, components);

        initialTreeGenerator = new InitialTreeGenerator(options, components);
        treeModelGenerator = new TreeModelGenerator(options, components);
        clockModelGenerator = new ClockModelGenerator(options, components);

        operatorsGenerator = new OperatorsGenerator(options, components);
        parameterPriorGenerator = new ParameterPriorGenerator(options, components);
        logGenerator = new LogGenerator(options, components);

    }

    /**
     * Checks various options to check they are valid. Throws IllegalArgumentExceptions with
     * descriptions of the problems.
     *
     * @throws IllegalArgumentException if there is a problem with the current settings
     */
    public void checkOptions() throws GeneratorException {
        try{

            //++++++++++++++++ Taxon List ++++++++++++++++++
            TaxonList taxonList = options.taxonList;
            Set<String> ids = new HashSet<String>();

            ids.add(TaxaParser.TAXA);
            ids.add(AlignmentParser.ALIGNMENT);

            if (taxonList != null) {
                if (taxonList.getTaxonCount() < 2) {
                    throw new GeneratorException("BEAST requires at least two taxa to run.");
                }

                for (int i = 0; i < taxonList.getTaxonCount(); i++) {
                    Taxon taxon = taxonList.getTaxon(i);
                    if (ids.contains(taxon.getId())) {
                        throw new GeneratorException("A taxon has the same id," + taxon.getId() +  MESSAGE_CAL);
                    }
                    ids.add(taxon.getId());
                }
            }

            //++++++++++++++++ Taxon Sets ++++++++++++++++++
            for (PartitionTreeModel model : options.getPartitionTreeModels()) {
                // should be only 1 calibrated internal node with a proper prior and monophyletic for each tree at moment
                if (model.getPartitionTreePrior().getNodeHeightPrior() == TreePriorType.YULE_CALIBRATION) {
                    if (options.treeModelOptions.isNodeCalibrated(model) < 0) // invalid node calibration
                        throw new GeneratorException(MESSAGE_CAL_YULE);

                    if (options.treeModelOptions.isNodeCalibrated(model) > 0) { // internal node calibration
                        List taxonSetsList = options.getKeysFromValue(options.taxonSetsTreeModel, model);
                        if (taxonSetsList.size() != 1 || !options.taxonSetsMono.get(taxonSetsList.get(0))) { // 1 tmrca per tree && monophyletic
                            throw new GeneratorException(MESSAGE_CAL_YULE, TreeLimeFrame.TAXON_SETS);
                        }
                    }
                }
            }

            for (Taxa taxa : options.taxonSets) {
                // AR - we should allow single taxon taxon sets...
                if (taxa.getTaxonCount() < 1 // && !options.taxonSetsIncludeStem.get(taxa)
                        ) {
                    throw new GeneratorException(
                            "Taxon set, " + taxa.getId() + ", should contain \n" +
                                    "at least one taxa. Please go back to Taxon Sets \n" +
                                    "panel to correct this.", TreeLimeFrame.TAXON_SETS);
                }
                if (ids.contains(taxa.getId())) {
                    throw new GeneratorException("A taxon set has the same id," + taxa.getId() +
                            MESSAGE_CAL, TreeLimeFrame.TAXON_SETS);
                }
                ids.add(taxa.getId());
            }

            //++++++++++++++++ Traits ++++++++++++++++++
            // missing data is not necessarily an issue...
//        for (TraitData trait : options.traits) {
//            for (int i = 0; i < trait.getTaxaCount(); i++) {
////                System.out.println("Taxon " + trait.getTaxon(i).getId() + " : [" + trait.getTaxon(i).getAttribute(trait.getName()) + "]");
//                if (!trait.hasValue(i))
//                    throw new IllegalArgumentException("Taxon " + trait.getTaxon(i).getId() +
//                            " has no value for Trait " + trait.getName());
//            }
//        }

            //++++++++++++++++ Tree Prior ++++++++++++++++++
//        if (options.isShareSameTreePrior()) {
            if (options.getPartitionTreeModels().size() > 1) { //TODO not allowed multi-prior yet
                for (PartitionTreePrior prior : options.getPartitionTreePriors()) {
                    if (prior.getNodeHeightPrior() == TreePriorType.GMRF_SKYRIDE) {
                        throw new GeneratorException("For the Skyride, tree model/tree prior combination not implemented by BEAST." +
                                "\nThe Skyride is only available for a single tree model partition in this release.", TreeLimeFrame.TREES);
                    }
                }
            }

            for (PartitionTreePrior prior : options.getPartitionTreePriors()) {
                if ((prior.getNodeHeightPrior() == TreePriorType.SKYGRID || prior.getNodeHeightPrior() == TreePriorType.SKYGRID_HMC)
                        && Double.isNaN(prior.getSkyGridInterval())) {
                    throw new GeneratorException("The Skygrid cut-off time must be set and greater than 0.0.", TreeLimeFrame.TREES);
                }
            }

            //+++++++++++++++ Starting tree ++++++++++++++++
            for (PartitionTreeModel model : options.getPartitionTreeModels()) {
                if (model.getStartingTreeType() == StartingTreeType.USER) {
                    if (model.getUserStartingTree() == null) {
                        throw new GeneratorException("Please select a starting tree in " + TreeLimeFrame.TREES + " panel, " +
                                "\nwhen choosing user specified starting tree option.", TreeLimeFrame.TREES);
                    }
                }
            }

            //++++++++++++++++ Random local clock model validation ++++++++++++++++++
            for (PartitionClockModel model : options.getPartitionClockModels()) {
                // 1 random local clock CANNOT have different tree models
                if (model.getClockType() == ClockType.RANDOM_LOCAL_CLOCK) { // || AUTOCORRELATED_LOGNORMAL
                    PartitionTreeModel treeModel = null;
                    for (AbstractPartitionData pd : options.getDataPartitions(model)) { // only the PDs linked to this tree model
                        if (treeModel != null && treeModel != pd.getPartitionTreeModel()) {
                            throw new GeneratorException("A single random local clock cannot be applied to multiple trees.", TreeLimeFrame.CLOCK_MODELS);
                        }
                        treeModel = pd.getPartitionTreeModel();
                    }
                }
            }

            //++++++++++++++++ Tree Model ++++++++++++++++++
            for (PartitionTreeModel model : options.getPartitionTreeModels()) {
                int numOfTaxa = -1;
                for (AbstractPartitionData pd : options.getDataPartitions(model)) {
                    if (pd.getTaxonCount() > 0) {
                        if (numOfTaxa > 0) {
                            if (numOfTaxa != pd.getTaxonCount()) {
                                throw new GeneratorException("Partitions with different taxa cannot share the same tree.", TreeLimeFrame.DATA_PARTITIONS);
                            }
                        } else {
                            numOfTaxa = pd.getTaxonCount();
                        }
                    }
                }
            }

            //++++++++++++++++ Prior Bounds ++++++++++++++++++
            for (Parameter param : options.selectParameters()) {
                if (!Double.isNaN(param.getInitial()) ) {
                    if (param.isTruncated && (param.getInitial() < param.truncationLower || param.getInitial() > param.truncationUpper)) {
                        throw new GeneratorException("Parameter \"" + param.getName() + "\":" +
                                "\ninitial value " + param.getInitial() + " is NOT in the range [" + param.truncationLower + ", " + param.truncationUpper + "]," +
                                "\nor this range is wrong. Please check the Prior panel.", TreeLimeFrame.PRIORS);
                    } else if (param.priorType == PriorType.UNIFORM_PRIOR && (param.getInitial() < param.uniformLower || param.getInitial() > param.uniformUpper)) {
                        throw new GeneratorException("Parameter \"" + param.getName() + "\":" +
                                "\ninitial value " + param.getInitial() + " is NOT in the range [" + param.uniformLower + ", " + param.uniformUpper + "]," +
                                "\nor this range is wrong. Please check the Prior panel.", TreeLimeFrame.PRIORS);
                    }
                    if (param.isNonNegative && param.getInitial() < 0.0) {
                        throw new GeneratorException("Parameter \"" + param.getName() + "\":" +
                                "\ninitial value " + param.getInitial() + " should be non-negative. Please check the Prior panel.", TreeLimeFrame.PRIORS);
                    }

                    if (param.isZeroOne && (param.getInitial() < 0.0 || param.getInitial() > 1.0)) {
                        throw new GeneratorException("Parameter \"" + param.getName() + "\":" +
                                "\ninitial value " + param.getInitial() + " should lie in the interval [0, 1]. Please check the Prior panel.", TreeLimeFrame.PRIORS);
                    }
                }
            }

            checkComponentOptions();

            // add other tests and warnings here
            // Speciation model with dated tips
            // Sampling rates without dated tips or priors on rate or nodes

        } catch (Exception e) {
            // catch any other exceptions here and rethrow to generate messages
            throw new GeneratorException(e.getMessage());
        }

        // outside the try above, which would drop the panel to switch to
        for (PartitionClockModel clock : options.getPartitionClockModels()) {
            if (!clock.getClockType().isMaximizable()) {
                throw new GeneratorException("TreeLime does not support the " + clock.getClockType().toString().toLowerCase() +
                        " (clock model " + clock.getName() + "). Use the shrinkage local clock for local rate shifts.",
                        TreeLimeFrame.CLOCK_MODELS);
            }
        }
        for (AbstractPartitionData data : options.dataPartitions) {
            double[] counts = data instanceof PartitionData ? ((PartitionData) data).getConstantSiteCounts() : null;
            if (counts == null) {
                continue;
            }
            if (data.getDataType().getType() != DataType.NUCLEOTIDES || counts.length != 4) {
                throw new GeneratorException("Constant-site counts (partition " + data.getName() + ") need " +
                        "nucleotide data and four counts, for A, C, G and T.", TreeLimeFrame.DATA_PARTITIONS);
            }
            for (double count : counts) {
                if (!(count >= 0) || Double.isInfinite(count) || count != Math.rint(count)) {
                    throw new GeneratorException("Constant-site counts (partition " + data.getName() + ") must be " +
                            "whole numbers, not negative.", TreeLimeFrame.DATA_PARTITIONS);
                }
            }
            if (data.getPartitionSubstitutionModel().getCodonPartitionCount() > 1) {
                throw new GeneratorException("Constant-site counts (partition " + data.getName() + ") cannot be " +
                        "split by codon position: turn off codon partitioning for this partition.",
                        TreeLimeFrame.DATA_PARTITIONS);
            }
        }
        String maxRunTime = options.maximizationOptions.maxRunTime;
        if (!maxRunTime.isEmpty()) {
            try {
                MaximizeWrtParameterParser.parseRunTime(maxRunTime);
            } catch (XMLParseException e) {
                throw new GeneratorException("Maximum run time: " + e.getMessage(), TreeLimeFrame.MAXIMIZATION);
            }
        }
    }
    public boolean checkUserTreeIsBifurcating(){
        for (PartitionTreeModel model : options.getPartitionTreeModels()) {
            if(model.getStartingTreeType()==StartingTreeType.USER){
                if(!(isBifurcatingTree(model.getUserStartingTree(),model.getUserStartingTree().getRoot()))){
                    return false;
                }
            }
        }
        return true;
    }
    public boolean isBifurcatingTree(Tree tree, NodeRef node) {
        if (tree.getChildCount(node) > 2) return false;
        for (int i = 0; i < tree.getChildCount(node); i++) {
            if (!isBifurcatingTree(tree, tree.getChild(node, i))) return false;
        }
        return true;
    }

    /**
     * Generate a beast xml file from these beast options
     *
     * @param file File
     * @throws java.io.IOException IOException
     * @throws dr.app.util.Arguments.ArgumentException
     *                             ArgumentException
     */
    public void generateXML(File file) throws GeneratorException, IOException, Arguments.ArgumentException {

        XMLWriter writer = new XMLWriter(new BufferedWriter(new FileWriter(file)));

        writer.writeText("<?xml version=\"1.0\" standalone=\"yes\"?>");
        writer.writeComment("Generated by TreeLime " + VERSION.getVersionString(),
                "      by Alexei J. Drummond, Andrew Rambaut and Marc A. Suchard",
                "      Department of Computer Science, University of Auckland and",
                "      Institute of Evolutionary Biology, University of Edinburgh",
                "      David Geffen School of Medicine, University of California, Los Angeles",
                "      http://beast.community/");
        writer.writeOpenTag("beast", new Attribute.Default<String>("version", TreeLimeApp.VERSION.getVersion()) );
        writer.writeText("");

        // this gives any added implementations of the 'Component' interface a
        // chance to generate XML at this point in the BEAST file.
        generateInsertionPoint(ComponentGenerator.InsertionPoint.BEFORE_TAXA, writer);

        if (options.originDate != null) {
            // Create a dummy taxon whose job is to specify the origin date
            Taxon originTaxon = new Taxon("originTaxon");
            options.originDate.setUnits(options.units);
            originTaxon.setDate(options.originDate);
            writeTaxon(originTaxon, true, false, writer);
        }


        //++++++++++++++++ Taxon List ++++++++++++++++++
        try {
            // write complete taxon list
            writeTaxa(options.taxonList, writer);

            writer.writeText("");

        } catch (Exception e) {
            e.printStackTrace(System.err);
            throw new GeneratorException("Taxon list generation has failed:\n" + e.getMessage());
        }

        //++++++++++++++++ Taxon Sets ++++++++++++++++++
        List<Taxa> taxonSets = options.taxonSets;
        try {
            if (taxonSets != null && !taxonSets.isEmpty()) {
                tmrcaStatisticsGenerator.writeTaxonSets(writer, taxonSets);
            }
        } catch (Exception e) {
            e.printStackTrace();
            throw new GeneratorException("Taxon sets generation has failed:\n" + e.getMessage());
        }

        generateInsertionPoint(ComponentGenerator.InsertionPoint.AFTER_TAXA, writer);

        //++++++++++++++++ Alignments ++++++++++++++++++
        List<Alignment> alignments = new ArrayList<Alignment>();
        try {
            for (AbstractPartitionData partition : options.dataPartitions) {
                Alignment alignment = null;
                if (partition instanceof PartitionData) {
                    alignment = ((PartitionData) partition).getAlignment();
                }
                if (alignment != null && !alignments.contains(alignment)) {
                    alignments.add(alignment);
                }
            }
            if (alignments.size() > 0) {
                alignmentGenerator.writeAlignments(alignments, writer);
                // a bootstrap replicate: each partition's resampled columns; the originals above still give the
                // starting tree, so the topology is the original fit's
                if (options.bootstrapReplicate != null) {
                    for (AbstractPartitionData partition : options.dataPartitions) {
                        Alignment replicate = partition instanceof PartitionData ?
                                options.bootstrapReplicate.getAlignment((PartitionData) partition) : null;
                        if (replicate != null) {
                            alignmentGenerator.writeBootstrapAlignment(replicate, writer);
                        }
                    }
                }
                generateInsertionPoint(ComponentGenerator.InsertionPoint.AFTER_SEQUENCES, writer);
            }
        } catch (Exception e) {
            e.printStackTrace();
            throw new GeneratorException("Alignments generation has failed:\n" + e.getMessage());
        }

        //++++++++++++++++ Pattern Lists ++++++++++++++++++
        try {
            // Construct pattern lists even if sampling from a null alignment
            //if (!options.samplePriorOnly) {
                for (AbstractPartitionData partition : options.dataPartitions) { // Each PD has one TreeLikelihood
                    if (partition.getTaxonList() != null) {
                        switch (partition.getDataType().getType()) {
                            case DataType.NUCLEOTIDES:
                            case DataType.AMINO_ACIDS:
                            case DataType.CODONS:
                            case DataType.COVARION:
                            case DataType.TWO_STATES:
                                patternListGenerator.writePatternList((PartitionData) partition, writer);
                                break;

                            case DataType.GENERAL:
                            case DataType.CONTINUOUS:
                                // no patternlist for trait data - discrete (general) data type uses an
                                // attribute patterns which is generated next bit of this method.
                                break;

                            case DataType.TREE:
                            case DataType.DUMMY:
                                //Do nothing
                                break;

                            default:
                                throw new IllegalArgumentException("Unsupported data type");
                        }
                        writer.writeText("");
                    }
                }
            //}
        } catch (Exception e) {
            e.printStackTrace();
            throw new GeneratorException("Pattern lists generation has failed:\n" + e.getMessage());
        }

        generateInsertionPoint(ComponentGenerator.InsertionPoint.AFTER_PATTERNS, writer);

        //++++++++++++++++ Tree Prior Model ++++++++++++++++++
        try {
            for (PartitionTreePrior prior : options.getPartitionTreePriors()) {
                treePriorGenerator.writeTreePriorModel(prior, writer);
                writer.writeText("");

                if (prior.getSubtreeTaxonSet() != null) {
                    treePriorGenerator.writeSubtreePriorModel(prior, writer);
                    writer.writeText("");
                }
            }
        } catch (Exception e) {
            e.printStackTrace();
            throw new GeneratorException("Tree prior model generation has failed:\n" + e.getMessage());
        }

        //++++++++++++++++ Starting Tree ++++++++++++++++++
        try {
            for (PartitionTreeModel model : options.getPartitionTreeModels()) {
                initialTreeGenerator.writeStartingTree(model, writer);
                writer.writeText("");
            }
        } catch (Exception e) {
            e.printStackTrace();
            throw new GeneratorException("Starting tree generation has failed:\n" + e.getMessage());
        }

        //++++++++++++++++ Tree Model +++++++++++++++++++
        try {
            for (PartitionTreeModel model : options.getPartitionTreeModels()) {
                treeModelGenerator.writeTreeModel(model, writer);
                writer.writeText("");
            }

            generateInsertionPoint(ComponentGenerator.InsertionPoint.AFTER_TREE_MODEL, writer);
        } catch (Exception e) {
            e.printStackTrace();
            throw new GeneratorException("Tree model generation has failed:\n" + e.getMessage());
        }

        //++++++++++++++++ Statistics ++++++++++++++++++
        try {
            if (taxonSets != null && taxonSets.size() > 0) {
                tmrcaStatisticsGenerator.writeTMRCAStatistics(writer);
            }
        } catch (Exception e) {
            e.printStackTrace();
            throw new GeneratorException("TMRCA statistics generation has failed:\n" + e.getMessage());
        }

        //++++++++++++++++ Tree Prior Likelihood ++++++++++++++++++
        try {
            for (PartitionTreeModel model : options.getPartitionTreeModels()) {
                treePriorGenerator.writePriorLikelihood(model, writer);
                writer.writeText("");
                if (model.getPartitionTreePrior().getSubtreeTaxonSet() != null) {
                    treePriorGenerator.writeSubtreePriorLikelihood(model, writer);
                    writer.writeText("");
                }
            }

            for (PartitionTreePrior prior : options.getPartitionTreePriors()) {
                treePriorGenerator.writeMultiLociTreePriors(prior, writer);
            }

            generateInsertionPoint(ComponentGenerator.InsertionPoint.AFTER_TREE_PRIOR, writer);
        } catch (Exception e) {
            e.printStackTrace();
            throw new GeneratorException("Tree prior likelihood generation has failed:\n" + e.getMessage());
        }

        //++++++++++++++++ Branch Rates Model ++++++++++++++++++
        try {
            for (PartitionClockModel model : options.getPartitionClockModels()) {
                clockModelGenerator.writeBranchRatesModel(model, writer);
                writer.writeText("");
            }
        } catch (Exception e) {
            e.printStackTrace();
            throw new GeneratorException("Branch rates model generation has failed:\n" + e.getMessage());
        }

        //++++++++++++++++ Substitution Model & Site Model ++++++++++++++++++
        try {
            for (PartitionSubstitutionModel model : options.getPartitionSubstitutionModels()) {
                substitutionModelGenerator.writeSubstitutionSiteModel(model, writer);
                writer.writeText("");
            }

            generateInsertionPoint(ComponentGenerator.InsertionPoint.AFTER_SUBSTITUTION_MODEL, writer);
        } catch (Exception e) {
            e.printStackTrace();
            throw new GeneratorException("Substitution model or site model generation has failed:\n" + e.getMessage());
        }

        //++++++++++++++++ AllMus parameter ++++++++++++++++++
        try {
            for (PartitionClockModel model : options.getPartitionClockModels()) {
                clockModelGenerator.writeAllMus(model, writer);
            }
        } catch (Exception e) {
            e.printStackTrace();
            throw new GeneratorException("Clock model generation has failed:\n" + e.getMessage());
        }

        //++++++++++++++++ Site Model ++++++++++++++++++
//        for (PartitionSubstitutionModel model : options.getPartitionSubstitutionModels()) {
//            substitutionModelGenerator.writeSiteModel(model, writer); // site model
//            substitutionModelGenerator.writeAllMus(model, writer); // allMus
//            writer.writeText("");
//        }

        generateInsertionPoint(ComponentGenerator.InsertionPoint.AFTER_SITE_MODEL, writer);

        //++++++++++++++++ Tree Likelihood ++++++++++++++++++
        try {
            Map<Pair<Pair<PartitionTreeModel, PartitionClockModel>, DataType>, List<PartitionData>> partitionLists = new HashMap<Pair<Pair<PartitionTreeModel, PartitionClockModel>, DataType>, List<PartitionData>>();
            options.multiPartitionLists.clear();
            options.otherPartitions.clear();

            for (AbstractPartitionData partition : options.dataPartitions) {
                // generate tree likelihoods for alignment data partitions
                if (partition.getTaxonList() != null) {

                    if (treeLikelihoodGenerator.canUseMultiPartition(partition)) {
                        // all sequence partitions of the same type as the first into the list for use in a
                        // MultipartitionTreeDataLikelihood. Must also share the same tree, clock model and not be doing
                        // ancestral reconstruction or counting
                        Pair<Pair<PartitionTreeModel, PartitionClockModel>, DataType> key = new Pair(new Pair(partition.getPartitionTreeModel(),partition.getPartitionClockModel()), partition.getDataType());
                        List<PartitionData> partitions = partitionLists.get(key);

                        if (partitions == null) {
                            partitions = new ArrayList<PartitionData>();
                            options.multiPartitionLists.add(partitions);
                        }

                        partitions.add((PartitionData) partition);
                        partitionLists.put(key, partitions);

                    } else {
                        options.otherPartitions.add(partition);
                    }


                }
            }

            treeLikelihoodGenerator.writeAllTreeLikelihoods(writer);

            generateInsertionPoint(ComponentGenerator.InsertionPoint.AFTER_TREE_LIKELIHOOD, writer);
        } catch (Exception e) {
            e.printStackTrace();
            throw new GeneratorException("Tree likelihood generation has failed:\n" + e.getMessage());
        }

        generateInsertionPoint(ComponentGenerator.InsertionPoint.AFTER_TRAITS, writer);

        //++++++++++++++++ Maximization ++++++++++++++++++
        try {
            writeObjective(writer);
            writer.writeText("");

            writeMaximizer(writer);
            writer.writeText("");
        } catch (Exception e) {
            e.printStackTrace();
            throw new GeneratorException("Maximizer generation has failed:\n" + e.getMessage());
        }

        writer.writeCloseTag("beast");
        writer.flush();
        writer.close();
    }

    /**
     * Generate a taxa block from these beast options
     *
     * @param writer    the writer
     * @param taxonList the taxon list to write
     * @throws dr.app.util.Arguments.ArgumentException
     *          ArgumentException
     */
    private void writeTaxa(TaxonList taxonList, XMLWriter writer) throws Arguments.ArgumentException {
        // -1 (single taxa), 0 (1st gene of multi-taxa)

        writer.writeComment("The list of taxa to be analysed (can also include dates/ages).",
                "ntax=" + taxonList.getTaxonCount());
        writer.writeOpenTag(TaxaParser.TAXA, new Attribute[]{new Attribute.Default<String>(XMLParser.ID, TaxaParser.TAXA)});

        boolean hasAttr = options.traits.size() > 0;

        boolean firstDate = true;
        for (int i = 0; i < taxonList.getTaxonCount(); i++) {
            Taxon taxon = taxonList.getTaxon(i);

            boolean hasDate = false;
            if (options.useTipDates) {
                hasDate = TaxonList.Utils.hasAttribute(taxonList, i, dr.evolution.util.Date.DATE);
            }

            if (hasDate) {
                dr.evolution.util.Date date = (dr.evolution.util.Date) taxon.getAttribute(dr.evolution.util.Date.DATE);

                if (firstDate) {
                    options.units = date.getUnits();
                    firstDate = false;
                } else {
                    if (options.units != date.getUnits()) {
                        System.err.println("Error: Units in dates do not match.");
                    }
                }
            }

            writeTaxon(taxon, hasDate, hasAttr, writer);

        }

        writer.writeCloseTag(TaxaParser.TAXA);
    }

    /**
     * Generate a taxa block from these beast options
     *
     * @param writer    the writer
     * @param taxon the taxon to write
     * @throws dr.app.util.Arguments.ArgumentException
     *          ArgumentException
     */
    private void writeTaxon(Taxon taxon, boolean hasDate, boolean hasAttr, XMLWriter writer) throws Arguments.ArgumentException {

        writer.writeTag(TaxonParser.TAXON, new Attribute[]{
                        new Attribute.Default<String>(XMLParser.ID, taxon.getId())},
                !(hasDate || hasAttr)); // false if any of hasDate or hasAttr is true


        if (hasDate) {
            dr.evolution.util.Date date = (dr.evolution.util.Date) taxon.getAttribute(dr.evolution.util.Date.DATE);

            Attribute[] attributes;
            if (date.getUncertainty() > 0.0) {
                attributes = new Attribute[] {
                        new Attribute.Default<Double>(DateParser.VALUE, date.getTimeValue()),
                        new Attribute.Default<String>(DateParser.DIRECTION, date.isBackwards() ? DateParser.BACKWARDS : DateParser.FORWARDS),
                        new Attribute.Default<String>(DateParser.UNITS, Units.Utils.getDefaultUnitName(options.units)),
                        new Attribute.Default<Double>(DateParser.UNCERTAINTY, date.getUncertainty())
                };
            } else {
                attributes = new Attribute[] {
                        new Attribute.Default<Double>(DateParser.VALUE, date.getTimeValue()),
                        new Attribute.Default<String>(DateParser.DIRECTION, date.isBackwards() ? DateParser.BACKWARDS : DateParser.FORWARDS),
                        new Attribute.Default<String>(DateParser.UNITS, Units.Utils.getDefaultUnitName(options.units))
                        //new Attribute.Default("origin", date.getOrigin()+"")
                };
            }

            writer.writeTag(dr.evolution.util.Date.DATE, attributes, true);
        }

        for (TraitData trait : options.traits) {
            // there is no harm in allowing the species trait to be listed in the taxa
//                if (!trait.getName().equalsIgnoreCase(TraitData.TRAIT_SPECIES)) {
            writer.writeOpenTag(AttributeParser.ATTRIBUTE, new Attribute[]{
                    new Attribute.Default<String>(Attribute.NAME, trait.getName())});

            // denotes missing data using '?'
            writer.writeText(taxon.containsAttribute(trait.getName()) ? taxon.getAttribute(trait.getName()).toString() : "?");
            writer.writeCloseTag(AttributeParser.ATTRIBUTE);
//                }
        }

        generateInsertionPoint(ComponentGenerator.InsertionPoint.IN_TAXON, taxon, writer);

        if (hasDate || hasAttr) writer.writeCloseTag(TaxonParser.TAXON);

    }

    /**
     * Write the timer report block.
     *
     * @param writer the writer
     */
    /**
     * Write the densities the maximizer can maximize, at the top level: prior, likelihood and their joint, as BEAUti
     * writes them inside the MCMC element.
     *
     * @param writer XMLWriter
     */
    public void writeObjective(XMLWriter writer) {
        writer.writeComment("The densities to maximize: priors, tree likelihoods and their joint");

        if (options.hasData()) {
            writer.writeOpenTag(CompoundLikelihoodParser.JOINT, new Attribute.Default<String>(XMLParser.ID, "joint"));
        }

        // write prior block
        writer.writeOpenTag(CompoundLikelihoodParser.PRIOR, new Attribute.Default<String>(XMLParser.ID, "prior"));

        parameterPriorGenerator.writeParameterPriors(writer);

        for (PartitionTreeModel model : options.getPartitionTreeModels()) {
            PartitionTreePrior prior = model.getPartitionTreePrior();
            treePriorGenerator.writePriorLikelihoodReference(prior, model, writer);
            writer.writeText("");
            if (prior.getSubtreeTaxonSet() != null) {
                treePriorGenerator.writeSubtreePriorLikelihoodReference(prior, model, writer);
                writer.writeText("");
            }
        }

        for (PartitionTreePrior prior : options.getPartitionTreePriors()) {
            treePriorGenerator.writeMultiLociLikelihoodReference(prior, writer);
            writer.writeText("");
        }

        clockModelGenerator.writeClockLikelihoodReferences(writer);

        generateInsertionPoint(ComponentGenerator.InsertionPoint.IN_MCMC_PRIOR, writer);

        writer.writeCloseTag(CompoundLikelihoodParser.PRIOR);

        if (options.hasData()) {
            // write likelihood block
            writer.writeOpenTag(CompoundLikelihoodParser.LIKELIHOOD, new Attribute.Default<String>(XMLParser.ID, "likelihood"));

            treeLikelihoodGenerator.writeTreeLikelihoodReferences(writer);

            generateInsertionPoint(ComponentGenerator.InsertionPoint.IN_MCMC_LIKELIHOOD, writer);

            writer.writeCloseTag(CompoundLikelihoodParser.LIKELIHOOD);

            writer.writeCloseTag(CompoundLikelihoodParser.JOINT);
        }
    }

    /**
     * Write the maximizer: L-BFGS over the selected parameters, each mapped to the real line by its transform, then a
     * report of the optimum. With analytic gradients, each parameter gets a gradient block: the tree-likelihood term
     * (analytic where BEAST has one, numeric otherwise) and a numeric term from the cheap prior compound, which covers
     * every prior on the parameter. The maximizer's value is the chosen density itself (<objective>), not the union
     * of the blocks' terms.
     *
     * @param writer XMLWriter
     */
    public void writeMaximizer(XMLWriter writer) throws GeneratorException {
        MaximizationOptions maximization = options.maximizationOptions;
        MaximizationOptions.Objective objective = maximization.objective;

        if (objective != MaximizationOptions.Objective.PRIOR && !options.hasData()) {
            throw new GeneratorException("There are no data to compute a likelihood: choose the prior-only objective.");
        }

        List<Parameter> parameters = maximization.selectOptimized(options.selectParameters());
        Map<Parameter, GradientRegistry.LikelihoodTerm> terms = new LinkedHashMap<Parameter, GradientRegistry.LikelihoodTerm>();
        if (maximization.isAnalytic()) {
            for (Parameter parameter : parameters) {
                terms.put(parameter, likelihoodTerm(parameter));
            }
            if (objective == MaximizationOptions.Objective.ML) {
                // the likelihood is flat in these: they are identified only by their priors
                List<Parameter> identified = new ArrayList<Parameter>();
                for (Parameter parameter : parameters) {
                    if (terms.get(parameter) != GradientRegistry.LikelihoodTerm.NONE) {
                        identified.add(parameter);
                    }
                }
                parameters = identified;
            }
        }

        // trees whose internal node heights are optimized, each a block of its own
        Map<PartitionTreeModel, GradientRegistry.NodeHeightPlan> heightPlans =
                new LinkedHashMap<PartitionTreeModel, GradientRegistry.NodeHeightPlan>();
        for (Parameter parameter : options.selectParameters()) {
            if (maximization.classify(parameter, false) == MaximizationOptions.Treatment.NODE_HEIGHTS) {
                PartitionOptions owner = GradientRegistry.findOwner(options, parameter);
                if (owner instanceof PartitionTreeModel) {
                    GradientRegistry.NodeHeightPlan plan =
                            GradientRegistry.nodeHeightPlan(options, (PartitionTreeModel) owner);
                    if (plan.isOptimized()) {
                        heightPlans.put((PartitionTreeModel) owner, plan);
                    }
                }
            }
        }

        // estimated nucleotide frequencies: blocks of their own on the simplex, with the analytic likelihood term only
        // (their prior is flat on the simplex); not under the prior-only objective, which does not depend on them
        List<Parameter> frequencies = new ArrayList<Parameter>();
        if (maximization.isAnalytic() && objective != MaximizationOptions.Objective.PRIOR) {
            for (Parameter parameter : maximization.selectFrequencies(options.selectParameters())) {
                if (likelihoodTerm(parameter).isAnalytic()) {
                    frequencies.add(parameter);
                }
            }
        }

        if (parameters.isEmpty() && frequencies.isEmpty() && heightPlans.isEmpty()) {
            throw new GeneratorException("No parameter is selected for optimization in the Maximization tab.");
        }

        int analyticCount = 0;
        for (GradientRegistry.LikelihoodTerm term : terms.values()) {
            if (term.isAnalytic()) {
                analyticCount++;
            }
        }
        writer.writeComment("Maximize the " + objective.getDensityId() + " density by L-BFGS, with " +
                        (maximization.isAnalytic() ? "analytic tree-likelihood gradients for " + analyticCount + " of " +
                                parameters.size() + " parameters" : "a numeric gradient") +
                                (frequencies.isEmpty() ? "" : "; " + frequencies.size() +
                                        " set(s) of frequencies on the simplex, with analytic gradients"),
                heightPlans.isEmpty() ? "The tree topology and node heights are fixed" :
                        "The tree topology is fixed; the internal node heights of " + heightPlans.size() +
                                " tree(s) are optimized as ratios on the real line");
        for (String warning : MaximizationWarnings.get(options)) {
            writer.writeComment("WARNING: " + warning);
        }

        if (!parameters.isEmpty()) {
            writer.writeOpenTag(CompoundParameterParser.COMPOUND_PARAMETER,
                    new Attribute.Default<String>(XMLParser.ID, OPTIMIZED_PARAMETERS));
            for (Parameter parameter : parameters) {
                writer.writeIDref(ParameterParser.PARAMETER, parameter.getName());
            }
            writer.writeCloseTag(CompoundParameterParser.COMPOUND_PARAMETER);
            writer.writeText("");
        }

        if (maximization.isAnalytic()) {
            writeGradientBlocks(parameters, terms, frequencies, heightPlans, objective, writer);
        }

        List<Attribute> attributes = new ArrayList<Attribute>();
        attributes.add(new Attribute.Default<String>(XMLParser.ID, MAXIMIZER));
        attributes.add(new Attribute.Default<Integer>("nIterations", Math.max(1, maximization.nIterations)));
        if (!maximization.maxRunTime.isEmpty()) {
            attributes.add(new Attribute.Default<String>("maxRunTime", maximization.maxRunTime));
        }
        if (maximization.maxEvaluations > 0) {
            attributes.add(new Attribute.Default<Long>("maxEvaluations", maximization.maxEvaluations));
        }
        attributes.add(new Attribute.Default<Boolean>("startAtCurrentState", true));
        attributes.add(new Attribute.Default<Boolean>("printToScreen", maximization.printToScreen));
        attributes.add(new Attribute.Default<Boolean>("includeJacobian", false));
        // a line search that fails far from the optimum (e.g. from a poor starting tree) restarts from the best point
        attributes.add(new Attribute.Default<Integer>("lineSearchRestarts", 5));
        // also stop when the density has stopped improving: zero-length branches (a node ratio running to 1) leave a
        // plateau where the gradient norm test may never pass
        attributes.add(new Attribute.Default<Integer>("past", 10));
        attributes.add(new Attribute.Default<Double>("delta", 1E-8));
        if (!maximization.isAnalytic()) {
            attributes.add(new Attribute.Default<String>("numericStepSize", maximization.stepSize.getName()));
        }
        writer.writeOpenTag("maximizeWrtParameter", attributes);

        if (maximization.isAnalytic()) {
            writer.writeOpenTag("objective");
            writer.writeIDref(objective.getDensityId(), objective.getDensityId());
            writer.writeCloseTag("objective");
            writer.writeIDref("compoundGradient", MAXIMIZER_GRADIENT);
        } else {
            writer.writeOpenTag("densityWrtParameter");
            writer.writeIDref(objective.getDensityId(), objective.getDensityId());
            writer.writeIDref(CompoundParameterParser.COMPOUND_PARAMETER, OPTIMIZED_PARAMETERS);
            writer.writeCloseTag("densityWrtParameter");
        }

        // the transforms in the order of the gradient blocks: the other parameters, the frequencies (stick-breaking),
        // then each tree's node heights
        boolean multivariate = !heightPlans.isEmpty() || !frequencies.isEmpty();
        if (multivariate) {
            writer.writeOpenTag("multivariateCompoundTransform");
        }
        if (!parameters.isEmpty()) {
            writeCompoundTransform(parameters, writer);
        }
        for (Parameter parameter : frequencies) {
            writer.writeTag(UNIT_SIMPLEX_TRANSFORM, new Attribute.Default<Integer>("dimension",
                    parameter.getOptions() instanceof PartitionSubstitutionModel ?
                            ((PartitionSubstitutionModel) parameter.getOptions()).getDataType().getStateCount() : 4), true);
        }
        for (PartitionTreeModel tree : heightPlans.keySet()) {
            writer.writeIDref("nodeHeightTransform", tree.getPrefix() + NODE_HEIGHT_TRANSFORM);
        }
        if (multivariate) {
            writer.writeCloseTag("multivariateCompoundTransform");
        }

        // e.g. a gamma shape whose optimum lies at 0: held at its bound instead of chased down a plateau
        for (Parameter parameter : parameters) {
            Double bound = maximization.getLowerBound(parameter);
            if (bound != null) {
                writer.writeOpenTag("holdAtLowerBound", new Attribute.Default<Double>("value", bound));
                writer.writeIDref(ParameterParser.PARAMETER, parameter.getName());
                writer.writeCloseTag("holdAtLowerBound");
            }
        }

        writer.writeCloseTag("maximizeWrtParameter");
        writer.writeText("");

        writer.writeComment("The optimum, the gradient there, and why the maximizer stopped");
        writer.writeOpenTag(Report.REPORT, new Attribute.Default<String>(Report.FILENAME,
                maximization.getReportFileName(options.fileNameStem)));
        writer.writeIDref("maximizeWrtParameter", MAXIMIZER);
        writer.writeCloseTag(Report.REPORT);
        writer.writeOpenTag(Report.REPORT);
        writer.writeIDref("maximizeWrtParameter", MAXIMIZER);
        writer.writeCloseTag(Report.REPORT);
        writer.writeText("");

        List<Parameter> optimized = new ArrayList<Parameter>(parameters);
        optimized.addAll(frequencies);
        writeTimeTreeExport(optimized, writer);
    }

    /**
     * Write the trees at the optimum through BEAST's tree loggers, logged once by the export: the time tree with the
     * clock's rates (<stem>.timetree.nex) and the same tree in substitutions (<stem>.subst.nex); then the export's
     * tables of dates, rates and parameters and the run record (<stem>.dates.tsv, .rates.tsv, .params.tsv,
     * .summary.txt).
     *
     * @param optimized the parameters the maximizer moves, node heights aside
     */
    private void writeTimeTreeExport(List<Parameter> optimized, XMLWriter writer) {
        MaximizationOptions maximization = options.maximizationOptions;

        // the clock of each tree; a tree with several clocks has no single rate per branch
        Map<PartitionTreeModel, PartitionClockModel> treeClocks = new LinkedHashMap<PartitionTreeModel, PartitionClockModel>();
        for (PartitionTreeModel tree : options.getPartitionTreeModels()) {
            List<PartitionClockModel> clocks = new ArrayList<PartitionClockModel>();
            for (PartitionClockModel clock : options.getPartitionClockModels()) {
                if (clock.getPartitionTreeModel() == tree) {
                    clocks.add(clock);
                }
            }
            treeClocks.put(tree, clocks.size() == 1 ? clocks.get(0) : null);
        }

        writer.writeComment("The trees at the optimum, logged once by timeTreeExport");
        List<String> treeLogs = new ArrayList<String>();
        for (PartitionTreeModel tree : options.getPartitionTreeModels()) {
            PartitionClockModel clock = treeClocks.get(tree);
            treeLogs.add(writeOptimumTreeLog(tree, clock, false, writer));
            if (clock != null) {
                treeLogs.add(writeOptimumTreeLog(tree, clock, true, writer));
            }
        }
        writer.writeText("");

        writer.writeComment("Log the trees; the dates and rates, the parameters and a record of the run");
        writer.writeOpenTag(TimeTreeExportParser.TIME_TREE_EXPORT, new Attribute[]{
                new Attribute.Default<String>(XMLParser.ID, TIME_TREE_EXPORT),
                new Attribute.Default<String>(TimeTreeExportParser.FILE_NAME_STEM, options.fileNameStem)});
        for (String treeLog : treeLogs) {
            writer.writeIDref(TreeLoggerParser.LOG_TREE, treeLog);
        }
        writer.writeIDref("maximizeWrtParameter", MAXIMIZER);
        if (options.hasData()) {
            writer.writeOpenTag(TimeTreeExportParser.TREE_LIKELIHOODS);
            writer.writeIDref(CompoundLikelihoodParser.LIKELIHOOD, "likelihood");
            writer.writeCloseTag(TimeTreeExportParser.TREE_LIKELIHOODS);
        }
        writer.writeOpenTag(TimeTreeExportParser.PRIORS);
        writer.writeIDref(CompoundLikelihoodParser.PRIOR, "prior");
        writer.writeCloseTag(TimeTreeExportParser.PRIORS);

        for (PartitionTreeModel tree : options.getPartitionTreeModels()) {
            writer.writeOpenTag(TimeTreeExportParser.TIME_TREE);
            writer.writeIDref(DefaultTreeModel.TREE_MODEL, tree.getPrefix() + DefaultTreeModel.TREE_MODEL);
            if (treeClocks.get(tree) != null) {
                ClockModelGenerator.writeBranchRatesModelRef(treeClocks.get(tree), writer);
            }
            for (Taxa taxa : options.taxonSets) {
                if (options.taxonSetsTreeModel.get(taxa) == tree) {
                    writer.writeIDref(TaxaParser.TAXA, taxa.getId());
                }
            }
            writer.writeCloseTag(TimeTreeExportParser.TIME_TREE);
        }

        // the parameters the maximizer does not move, by reason; node heights are reported with the trees
        Map<String, List<Parameter>> fixed = new LinkedHashMap<String, List<Parameter>>();
        for (Parameter parameter : options.selectParameters()) {
            if (parameter.isStatistic || parameter.isDiscrete || parameter.isNodeHeight || optimized.contains(parameter)) {
                continue;
            }
            MaximizationOptions.Treatment treatment = maximization.classify(parameter, false);
            String reason = treatment.isOptimized() ? "not identified by the objective" : treatment.getReason();
            if (!fixed.containsKey(reason)) {
                fixed.put(reason, new ArrayList<Parameter>());
            }
            fixed.get(reason).add(parameter);
        }
        for (Map.Entry<String, List<Parameter>> group : fixed.entrySet()) {
            writer.writeOpenTag(TimeTreeExportParser.NOT_OPTIMIZED,
                    new Attribute.Default<String>(TimeTreeExportParser.REASON, group.getKey()));
            for (Parameter parameter : group.getValue()) {
                writer.writeIDref(ParameterParser.PARAMETER, parameter.getName());
            }
            writer.writeCloseTag(TimeTreeExportParser.NOT_OPTIMIZED);
        }
        writer.writeCloseTag(TimeTreeExportParser.TIME_TREE_EXPORT);
    }

    /**
     * A tree logger for the tree at the optimum: annotated with the clock's rates and, as the tree's comment, the
     * maximized log density (lnP).
     *
     * @return its id
     */
    private String writeOptimumTreeLog(PartitionTreeModel tree, PartitionClockModel clock, boolean substitutions,
                                       XMLWriter writer) {
        String id = tree.getPrefix() + (substitutions ? SUBST_TREE_LOG : TIME_TREE_LOG);
        List<Attribute> attributes = new ArrayList<Attribute>();
        attributes.add(new Attribute.Default<String>(XMLParser.ID, id));
        attributes.add(new Attribute.Default<Integer>(TreeLoggerParser.LOG_EVERY, 1));
        attributes.add(new Attribute.Default<Boolean>(TreeLoggerParser.NEXUS_FORMAT, true));
        attributes.add(new Attribute.Default<Boolean>(TreeLoggerParser.SORT_TRANSLATION_TABLE, true));
        attributes.add(new Attribute.Default<String>(TreeLoggerParser.FILE_NAME, options.fileNameStem + "." +
                tree.getPrefix() + (substitutions ? "subst.nex" : "timetree.nex")));
        if (substitutions) {
            attributes.add(new Attribute.Default<String>(TreeLoggerParser.BRANCH_LENGTHS, TreeLoggerParser.SUBSTITUTIONS));
        }
        writer.writeOpenTag(TreeLoggerParser.LOG_TREE, attributes);
        writer.writeIDref(DefaultTreeModel.TREE_MODEL, tree.getPrefix() + DefaultTreeModel.TREE_MODEL);
        if (clock != null) {
            ClockModelGenerator.writeBranchRatesModelRef(clock, writer);
        }
        String density = options.maximizationOptions.objective.getDensityId();
        writer.writeIDref(density, density);
        writer.writeCloseTag(TreeLoggerParser.LOG_TREE);
        return id;
    }

    // one segment per run of parameters with the same transform, in the order of the compound parameter
    private void writeCompoundTransform(List<Parameter> parameters, XMLWriter writer) {
        MaximizationOptions maximization = options.maximizationOptions;
        writer.writeOpenTag("compoundTransform");
        int start = 0;
        while (start < parameters.size()) {
            String type = maximization.classify(parameters.get(start), false).getTransform();
            writer.writeOpenTag(TransformParsers.TRANSFORM, new Attribute.Default<String>(TransformParsers.TYPE, type));
            int end = start;
            while (end < parameters.size() && maximization.classify(parameters.get(end), false).getTransform().equals(type)) {
                writer.writeIDref(ParameterParser.PARAMETER, parameters.get(end).getName());
                end++;
            }
            writer.writeCloseTag(TransformParsers.TRANSFORM);
            start = end;
        }
        writer.writeIDref(CompoundParameterParser.COMPOUND_PARAMETER, OPTIMIZED_PARAMETERS);
        writer.writeCloseTag("compoundTransform");
    }

    /**
     * The tree-likelihood term for a parameter; analytic terms need every partition in a tree data likelihood
     * (BEAUti's older tree likelihoods and trait likelihoods have no gradients).
     */
    private GradientRegistry.LikelihoodTerm likelihoodTerm(Parameter parameter) {
        GradientRegistry.LikelihoodTerm term = GradientRegistry.likelihoodTerm(options, parameter);
        if (term.isAnalytic() && (!options.otherPartitions.isEmpty() || !options.traits.isEmpty())) {
            return GradientRegistry.LikelihoodTerm.NUMERIC;
        }
        return term;
    }

    private void writeGradientBlocks(List<Parameter> parameters, Map<Parameter, GradientRegistry.LikelihoodTerm> terms,
                                     List<Parameter> frequencies,
                                     Map<PartitionTreeModel, GradientRegistry.NodeHeightPlan> heightPlans,
                                     MaximizationOptions.Objective objective, XMLWriter writer) {
        writer.writeComment("Gradient blocks, one per optimized parameter, in the order of " + OPTIMIZED_PARAMETERS);
        Map<Parameter, String> tags = new LinkedHashMap<Parameter, String>();
        for (Parameter parameter : parameters) {
            GradientRegistry.LikelihoodTerm term = terms.get(parameter);
            boolean likelihoodTerm = objective != MaximizationOptions.Objective.PRIOR &&
                    term != GradientRegistry.LikelihoodTerm.NONE;
            boolean priorTerm = objective != MaximizationOptions.Objective.ML;

            List<String> treeLikelihoods = likelihoodTerm && term.isAnalytic() ? treeLikelihoodsOf(parameter) :
                    new ArrayList<String>();
            int members = (priorTerm ? 1 : 0) + (likelihoodTerm ? Math.max(1, treeLikelihoods.size()) : 0);

            String id = "gradient." + parameter.getName();
            if (members > 1) {
                tags.put(parameter, "jointGradient");
                writer.writeOpenTag("jointGradient", new Attribute.Default<String>(XMLParser.ID, id));
            }
            // the prior term first, so the block's parameter is the parameter itself
            String tag = null;
            if (priorTerm) {
                tag = writeNumericalGradient(members > 1 ? null : id, "prior", parameter, writer);
            }
            if (likelihoodTerm) {
                if (!term.isAnalytic()) {
                    tag = writeNumericalGradient(members > 1 ? null : id, "likelihood", parameter, writer);
                } else {
                    for (String treeLikelihood : treeLikelihoods) {
                        tag = writeLikelihoodGradient(members > 1 ? null : id, term, treeLikelihood, parameter, writer);
                    }
                }
            }
            if (members == 1) {
                tags.put(parameter, tag);
            }
            if (members > 1) {
                writer.writeCloseTag("jointGradient");
            }
        }
        writer.writeText("");

        // frequencies: the analytic likelihood term only (SubstitutionModelGradient behind
        // branchSubstitutionParameterGradient), summed over the tree likelihoods that use the model
        for (Parameter parameter : frequencies) {
            String id = "gradient." + parameter.getName();
            List<String> treeLikelihoods = treeLikelihoodsOf(parameter);
            boolean joint = treeLikelihoods.size() > 1;
            if (joint) {
                writer.writeOpenTag("jointGradient", new Attribute.Default<String>(XMLParser.ID, id));
            }
            String tag = null;
            for (String treeLikelihood : treeLikelihoods) {
                tag = writeLikelihoodGradient(joint ? null : id, GradientRegistry.LikelihoodTerm.SUBSTITUTION,
                        treeLikelihood, parameter, writer);
            }
            if (joint) {
                writer.writeCloseTag("jointGradient");
                tag = "jointGradient";
            }
            tags.put(parameter, tag);
        }

        List<String> heightTags = new ArrayList<String>();
        for (Map.Entry<PartitionTreeModel, GradientRegistry.NodeHeightPlan> entry : heightPlans.entrySet()) {
            heightTags.add(writeNodeHeightBlock(entry.getKey(), entry.getValue(), objective, writer));
        }

        writer.writeOpenTag("compoundGradient", new Attribute.Default<String>(XMLParser.ID, MAXIMIZER_GRADIENT));
        for (Parameter parameter : parameters) {
            writer.writeIDref(tags.get(parameter), "gradient." + parameter.getName());
        }
        for (Parameter parameter : frequencies) {
            writer.writeIDref(tags.get(parameter), "gradient." + parameter.getName());
        }
        int index = 0;
        for (PartitionTreeModel tree : heightPlans.keySet()) {
            writer.writeIDref(heightTags.get(index++), "gradient." + tree.getPrefix() + NODE_HEIGHTS);
        }
        writer.writeCloseTag("compoundGradient");
        writer.writeText("");
    }

    // the tree likelihoods whose partitions use the parameter's model
    private List<String> treeLikelihoodsOf(Parameter parameter) {
        List<String> treeLikelihoods = new ArrayList<String>();
        for (List<PartitionData> partitions : options.multiPartitionLists) {
            if (GradientRegistry.usesOwner(options, partitions, parameter)) {
                treeLikelihoods.add(partitions.get(0).getPartitionTreeModel().getPrefix() +
                        partitions.get(0).getPartitionClockModel().getPrefix() + "treeLikelihood");
            }
        }
        return treeLikelihoods;
    }

    // returns the tag written
    private String writeNumericalGradient(String id, String density, Parameter parameter, XMLWriter writer) {
        return writeNumericalGradient(id, density, parameter.getName(), ParameterParser.PARAMETER, writer);
    }

    private String writeNumericalGradient(String id, String density, String parameterId, String parameterTag,
                                          XMLWriter writer) {
        writer.writeOpenTag("numericalGradient", id == null ? new Attribute[0] :
                new Attribute[]{new Attribute.Default<String>(XMLParser.ID, id)});
        writer.writeIDref(density, density);
        writer.writeIDref(parameterTag, parameterId);
        writer.writeCloseTag("numericalGradient");
        return "numericalGradient";
    }

    // returns the tag written
    private String writeLikelihoodGradient(String id, GradientRegistry.LikelihoodTerm term, String treeLikelihood,
                                         Parameter parameter, XMLWriter writer) {
        List<Attribute> attributes = new ArrayList<Attribute>();
        if (id != null) {
            attributes.add(new Attribute.Default<String>(XMLParser.ID, id));
        }
        switch (term) {
            case SUBSTITUTION:
                attributes.add(new Attribute.Default<Boolean>("homogeneous", true));
                writer.writeOpenTag("branchSubstitutionParameterGradient", attributes);
                writer.writeIDref("treeDataLikelihood", treeLikelihood);
                writer.writeIDref(ParameterParser.PARAMETER, parameter.getName());
                writer.writeCloseTag("branchSubstitutionParameterGradient");
                return "branchSubstitutionParameterGradient";
            case SITE_RATE:
                writer.writeOpenTag("siteRateModelGradient", attributes);
                writer.writeIDref("treeDataLikelihood", treeLikelihood);
                writer.writeIDref(ParameterParser.PARAMETER, parameter.getName());
                writer.writeCloseTag("siteRateModelGradient");
                return "siteRateModelGradient";
            case BRANCH_RATE:
                attributes.add(new Attribute.Default<String>("traitName", "Sequence"));
                writer.writeOpenTag("branchRateGradient", attributes);
                writer.writeIDref("treeDataLikelihood", treeLikelihood);
                writer.writeCloseTag("branchRateGradient");
                return "branchRateGradient";
            case LOCATION_SCALE:
                // as the hmc_clock XMLs differentiate arbitraryBranchRates' location and scale
                attributes.add(new Attribute.Default<String>("traitName", "Sequence"));
                attributes.add(new Attribute.Default<Boolean>("useHessian", false));
                writer.writeOpenTag("locationScaleGradient", attributes);
                writer.writeIDref("treeDataLikelihood", treeLikelihood);
                String role = parameter.getBaseName().equals(ClockType.HMCLN_SCALE) ? "scale" : "location";
                writer.writeOpenTag(role);
                writer.writeIDref(ParameterParser.PARAMETER, parameter.getName());
                writer.writeCloseTag(role);
                writer.writeCloseTag("locationScaleGradient");
                return "locationScaleGradient";
            default:
                throw new IllegalArgumentException("Not an analytic term: " + term);
        }
    }

    /**
     * The node heights of one tree, as in the hmc_nodeHeight set-up: a proxy parameter of the root and internal
     * heights, the ratio transform mapped to the real line (log root height, logit ratios), and a gradient block of
     * the tree likelihoods' node-height gradient with the prior side (the tree prior's and CTMC rate reference
     * priors' analytic gradients, or finite differences of the prior compound).
     *
     * @return the tag of the gradient block
     */
    private String writeNodeHeightBlock(PartitionTreeModel tree, GradientRegistry.NodeHeightPlan plan,
                                        MaximizationOptions.Objective objective, XMLWriter writer) {
        String prefix = tree.getPrefix();
        String treeModelId = prefix + DefaultTreeModel.TREE_MODEL;
        String proxyId = prefix + NODE_HEIGHTS;

        writer.writeComment("The internal node heights of " + treeModelId + ": root and internal nodes, mapped to " +
                "the real line through ratios");
        writer.writeOpenTag("nodeHeightProxyParameter", new Attribute[]{
                new Attribute.Default<String>(XMLParser.ID, proxyId),
                new Attribute.Default<Boolean>("rootNode", true)});
        writer.writeIDref(DefaultTreeModel.TREE_MODEL, treeModelId);
        writer.writeCloseTag("nodeHeightProxyParameter");

        writer.writeOpenTag("nodeHeightTransform", new Attribute[]{
                new Attribute.Default<String>(XMLParser.ID, prefix + NODE_HEIGHT_TRANSFORM),
                new Attribute.Default<Boolean>("withRoot", true)});
        writer.writeIDref(DefaultTreeModel.TREE_MODEL, treeModelId);
        ClockModelGenerator.writeBranchRatesModelRef(clockOf(tree), writer);
        writer.writeOpenTag("ratios", new Attribute.Default<Boolean>("realLine", true));
        writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default<String>(XMLParser.ID, prefix + "nodeHeight.ratios"), true);
        writer.writeCloseTag("ratios");
        writer.writeCloseTag("nodeHeightTransform");

        List<String> treeLikelihoods = new ArrayList<String>();
        for (List<PartitionData> partitions : options.multiPartitionLists) {
            if (partitions.get(0).getPartitionTreeModel() == tree) {
                treeLikelihoods.add(prefix + partitions.get(0).getPartitionClockModel().getPrefix() + "treeLikelihood");
            }
        }
        boolean otherLikelihoods = false;
        for (AbstractPartitionData partition : options.otherPartitions) {
            if (partition.getPartitionTreeModel() == tree) {
                otherLikelihoods = true;
            }
        }

        // the first member's parameter is the block's: the proxy for the numeric terms, the tree's heights otherwise
        String id = "gradient." + proxyId;
        writer.writeOpenTag("jointGradient", new Attribute.Default<String>(XMLParser.ID, id));
        if (objective != MaximizationOptions.Objective.ML) {
            if (plan.numericPrior) {
                writeNumericalGradient(null, "prior", proxyId, "nodeHeightProxyParameter", writer);
            } else {
                if (plan.treePriorGradient.equals("coalescentGradient") ||
                        plan.treePriorGradient.equals("speciationLikelihoodGradient")) {
                    writer.writeOpenTag(plan.treePriorGradient, plan.treePriorGradient.equals("coalescentGradient") ?
                            new Attribute[0] : new Attribute[]{new Attribute.Default<String>("wrtParameter", "nodeHeight")});
                    writer.writeIDref(plan.treePriorLikelihood, plan.treePriorId);
                    writer.writeIDref(DefaultTreeModel.TREE_MODEL, treeModelId);
                } else {
                    writer.writeOpenTag(plan.treePriorGradient, new Attribute.Default<String>("wrtParameter", "nodeHeight"));
                    writer.writeIDref(plan.treePriorLikelihood, plan.treePriorId);
                }
                writer.writeCloseTag(plan.treePriorGradient);

                // a copy of each CTMC rate reference prior on the tree's length, as ParameterPriorGenerator writes it
                for (Parameter rate : plan.ctmcPriors) {
                    writer.writeOpenTag("ctmcScalePriorNodeHeightGradient");
                    writer.writeOpenTag(CTMCScalePriorParser.MODEL_NAME,
                            new Attribute.Default<String>(XMLParser.ID, rate.getName() + ".nodeHeightPrior"));
                    writer.writeOpenTag(CTMCScalePriorParser.SCALEPARAMETER);
                    writer.writeIDref(ParameterParser.PARAMETER, rate.getName());
                    writer.writeCloseTag(CTMCScalePriorParser.SCALEPARAMETER);
                    writer.writeIDref(DefaultTreeModel.TREE_MODEL, treeModelId);
                    writer.writeCloseTag(CTMCScalePriorParser.MODEL_NAME);
                    writer.writeIDref(DefaultTreeModel.TREE_MODEL, treeModelId);
                    writer.writeCloseTag("ctmcScalePriorNodeHeightGradient");
                }
            }
        }
        if (objective != MaximizationOptions.Objective.PRIOR) {
            for (String treeLikelihood : treeLikelihoods) {
                writer.writeOpenTag("nodeHeightGradient", new Attribute.Default<String>("traitName", "Sequence"));
                writer.writeIDref("treeDataLikelihood", treeLikelihood);
                writer.writeCloseTag("nodeHeightGradient");
            }
            if (otherLikelihoods) {
                writeNumericalGradient(null, "likelihood", proxyId, "nodeHeightProxyParameter", writer);
            }
        }
        writer.writeCloseTag("jointGradient");
        writer.writeText("");
        return "jointGradient";
    }

    private PartitionClockModel clockOf(PartitionTreeModel tree) {
        for (AbstractPartitionData partition : options.dataPartitions) {
            if (partition.getPartitionTreeModel() == tree && partition.getPartitionClockModel() != null) {
                return partition.getPartitionClockModel();
            }
        }
        throw new IllegalArgumentException("No clock model for tree " + tree.getName());
    }

    public static final String NODE_HEIGHTS = "nodeHeights";
    public static final String NODE_HEIGHT_TRANSFORM = "nodeHeightTransform";
    public static final String MAXIMIZER_GRADIENT = "maximizerGradient";
    private static final String UNIT_SIMPLEX_TRANSFORM = "UnitSimplexTransform";
    public static final String MAXIMIZER = "maximizer";
    public static final String TIME_TREE_EXPORT = "timeTreeExport";
    public static final String TIME_TREE_LOG = "timeTreeLog";
    public static final String SUBST_TREE_LOG = "substTreeLog";
    public static final String OPTIMIZED_PARAMETERS = "optimizedParameters";

}
