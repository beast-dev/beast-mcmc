/*
 * BeastGenerator.java
 *
 * Copyright (c) 2002-2011 Alexei Drummond, Andrew Rambaut and Marc Suchard
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
 */

package dr.app.beauti.generator;

import dr.app.beast.BeastVersion;
import dr.app.beauti.BeautiFrame;
import dr.app.beauti.components.ComponentFactory;
import dr.app.beauti.options.*;
import dr.app.beauti.types.*;
import dr.app.beauti.util.XMLWriter;
import dr.app.util.Arguments;
import dr.evolution.alignment.Alignment;
import dr.evolution.datatype.ContinuousDataType;
import dr.evolution.datatype.DataType;
import dr.evolution.datatype.GeneralDataType;
import dr.evolution.datatype.Microsatellite;
import dr.evolution.util.Taxa;
import dr.evolution.util.Taxon;
import dr.evolution.util.TaxonList;
import dr.evolution.util.Units;
import dr.evomodelxml.speciation.MultiSpeciesCoalescentParser;
import dr.evomodelxml.speciation.SpeciationLikelihoodParser;
import dr.evoxml.AlignmentParser;
import dr.evoxml.DateParser;
import dr.evoxml.TaxaParser;
import dr.evoxml.TaxonParser;
import dr.inferencexml.distribution.MixedDistributionLikelihoodParser;
import dr.inferencexml.model.CompoundLikelihoodParser;
import dr.inferencexml.model.CompoundParameterParser;
import dr.inferencexml.operators.SimpleOperatorScheduleParser;
import dr.util.Attribute;
import dr.util.Version;
import dr.xml.AttributeParser;
import dr.xml.XMLParser;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * This class holds all the data for the current BEAUti Document
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @author Walter Xie
 * @version $Id: BeastGenerator.java,v 1.4 2006/09/05 13:29:34 rambaut Exp $
 */
public class BeastGenerator extends Generator {

    private final static Version version = new BeastVersion();
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
    private final BranchRatesModelGenerator branchRatesModelGenerator;
    private final OperatorsGenerator operatorsGenerator;
    private final ParameterPriorGenerator parameterPriorGenerator;
    private final LogGenerator logGenerator;
    //    private final DiscreteTraitGenerator discreteTraitGenerator;
    private final STARBEASTGenerator starBeastGenerator;
    private final TMRCAStatisticsGenerator tmrcaStatisticsGenerator;

    public BeastGenerator(BeautiOptions options, ComponentFactory[] components) {
        super(options, components);

        alignmentGenerator = new AlignmentGenerator(options, components);
        patternListGenerator = new PatternListGenerator(options, components);
        tmrcaStatisticsGenerator = new TMRCAStatisticsGenerator(options, components);
        substitutionModelGenerator = new SubstitutionModelGenerator(options, components);
        treePriorGenerator = new TreePriorGenerator(options, components);
        treeLikelihoodGenerator = new TreeLikelihoodGenerator(options, components);

        initialTreeGenerator = new InitialTreeGenerator(options, components);
        treeModelGenerator = new TreeModelGenerator(options, components);
        branchRatesModelGenerator = new BranchRatesModelGenerator(options, components);

        operatorsGenerator = new OperatorsGenerator(options, components);
        parameterPriorGenerator = new ParameterPriorGenerator(options, components);
        logGenerator = new LogGenerator(options, components);

        // this has moved into the component system...
//        discreteTraitGenerator = new DiscreteTraitGenerator(options, components);
        starBeastGenerator = new STARBEASTGenerator(options, components);
    }

    /**
     * Checks various options to check they are valid. Throws IllegalArgumentExceptions with
     * descriptions of the problems.
     *
     * @throws IllegalArgumentException if there is a problem with the current settings
     */
    public void checkOptions() throws GeneratorException {
        //++++++++++++++++ Taxon List ++++++++++++++++++
        TaxonList taxonList = options.taxonList;
        Set<String> ids = new HashSet<String>();

        ids.add(TaxaParser.TAXA);
        ids.add(AlignmentParser.ALIGNMENT);
        ids.add(TraitData.TRAIT_SPECIES);

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
                        throw new GeneratorException(MESSAGE_CAL_YULE, BeautiFrame.TAXON_SETS);
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
                                "panel to correct this.", BeautiFrame.TAXON_SETS);
            }
            if (ids.contains(taxa.getId())) {
                throw new GeneratorException("A taxon set has the same id," + taxa.getId() +
                        MESSAGE_CAL, BeautiFrame.TAXON_SETS);
            }
            ids.add(taxa.getId());
        }

        //++++++++++++++++ *BEAST ++++++++++++++++++
        if (options.useStarBEAST) {
            if (!options.traitExists(TraitData.TRAIT_SPECIES))
                throw new GeneratorException("A trait labelled \"species\" is required for *BEAST species designations." +
                        "\nPlease create or import the species designations in the Traits table.", BeautiFrame.TRAITS);

            //++++++++++++++++ Species Sets ++++++++++++++++++
            // should be only 1 calibrated internal node with monophyletic at moment
            if (options.getPartitionTreePriors().get(0).getNodeHeightPrior() == TreePriorType.SPECIES_YULE_CALIBRATION) {
                if (options.speciesSets.size() != 1 || !options.speciesSetsMono.get(options.speciesSets.get(0))) {
                    throw new GeneratorException(MESSAGE_CAL_YULE, BeautiFrame.TAXON_SETS);
                }
            }

            for (Taxa species : options.speciesSets) {
                if (species.getTaxonCount() < 2) {
                    throw new GeneratorException("Species set, " + species.getId() + ",\n should contain" +
                            "at least two species. \nPlease go back to Species Sets panel to select included species.", BeautiFrame.TAXON_SETS);
                }
                if (ids.contains(species.getId())) {
                    throw new GeneratorException("A species set has the same id," + species.getId() +
                            MESSAGE_CAL, BeautiFrame.TAXON_SETS);
                }
                ids.add(species.getId());
            }

            int tId = options.starBEASTOptions.getEmptySpeciesIndex();
            if (tId >= 0) {
                throw new GeneratorException("The taxon " + options.taxonList.getTaxonId(tId) +
                        " has NULL value for \"species\" trait", BeautiFrame.TRAITS);
            }
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
                    throw new GeneratorException("For GMRF, tree model/tree prior combination not implemented by BEAST yet" +
                            "\nIt is only available for single tree model partition for this release.", BeautiFrame.TREES);
                }
            }
        }

        //+++++++++++++++ Starting tree ++++++++++++++++
        for (PartitionTreeModel model : options.getPartitionTreeModels()) {
            if (model.getStartingTreeType() == StartingTreeType.USER) {
                if (model.getUserStartingTree() == null) {
                    throw new GeneratorException("Please select a starting tree in " + BeautiFrame.TREES + " panel, " +
                            "\nwhen choosing user specified starting tree option.", BeautiFrame.TREES);
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
                        throw new GeneratorException("A single random local clock cannot be applied to multiple trees.", BeautiFrame.CLOCK_MODELS);
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
                            throw new GeneratorException("Partitions with different taxa cannot share the same tree.", BeautiFrame.DATA_PARTITIONS);
                        }
                    } else {
                        numOfTaxa = pd.getTaxonCount();
                    }
                }
            }
        }

        //++++++++++++++++ Prior Bounds ++++++++++++++++++
        for (Parameter param : options.selectParameters()) {
            if (param.initial != Double.NaN) {
                if (param.isTruncated && (param.initial < param.truncationLower || param.initial > param.truncationUpper)) {
                    throw new GeneratorException("Parameter \"" + param.getName() + "\":" +
                            "\ninitial value " + param.initial + " is NOT in the range [" + param.truncationLower + ", " + param.truncationUpper + "]," +
                            "\nor this range is wrong. Please check the Prior panel.", BeautiFrame.PRIORS);
                } else if (param.priorType == PriorType.UNIFORM_PRIOR && (param.initial < param.uniformLower || param.initial > param.uniformUpper)) {
                    throw new GeneratorException("Parameter \"" + param.getName() + "\":" +
                            "\ninitial value " + param.initial + " is NOT in the range [" + param.uniformLower + ", " + param.uniformUpper + "]," +
                            "\nor this range is wrong. Please check the Prior panel.", BeautiFrame.PRIORS);
                }
                if (param.isNonNegative && param.initial < 0.0) {
                    throw new GeneratorException("Parameter \"" + param.getName() + "\":" +
                            "\ninitial value " + param.initial + " should be non-negative. Please check the Prior panel.", BeautiFrame.PRIORS);
                }

                if (param.isZeroOne && (param.initial < 0.0 || param.initial > 1.0)) {
                    throw new GeneratorException("Parameter \"" + param.getName() + "\":" +
                            "\ninitial value " + param.initial + " should lie in the interval [0, 1]. Please check the Prior panel.", BeautiFrame.PRIORS);
                }
            }
        }

        // add other tests and warnings here
        // Speciation model with dated tips
        // Sampling rates without dated tips or priors on rate or nodes

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
        writer.writeComment("Generated by BEAUTi " + version.getVersionString(),
                "      by Alexei J. Drummond, Andrew Rambaut and Marc A. Suchard",
                "      Department of Computer Science, University of Auckland and",
                "      Institute of Evolutionary Biology, University of Edinburgh",
                "      David Geffen School of Medicine, University of California, Los Angeles",
                "      http://beast.bio.ed.ac.uk/");
        writer.writeOpenTag("beast");
        writer.writeText("");

        // this gives any added implementations of the 'Component' interface a
        // chance to generate XML at this point in the BEAST file.
        generateInsertionPoint(ComponentGenerator.InsertionPoint.BEFORE_TAXA, writer);

        //++++++++++++++++ Taxon List ++++++++++++++++++
        try {
            // write complete taxon list
            writeTaxa(options.taxonList, writer);

            writer.writeText("");

            if (!options.hasIdenticalTaxa()) {
                // write all taxa in each gene tree regarding each data partition,
                for (AbstractPartitionData partition : options.dataPartitions) {
                    if (partition.getTaxonList() != null) {
                        writeDifferentTaxa(partition, writer);
                    }
                }
            }
        } catch (Exception e) {
            System.err.println(e);
            throw new GeneratorException("Taxon list generation has failed:\n" + e.getMessage());
        }

        //++++++++++++++++ Taxon Sets ++++++++++++++++++
        List<Taxa> taxonSets = options.taxonSets;
        try {
            if (taxonSets != null && taxonSets.size() > 0 && !options.useStarBEAST) {
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
                if (partition instanceof PartitionData) { // microsat has no alignment
                    alignment = ((PartitionData) partition).getAlignment();
                }
                if (alignment != null && !alignments.contains(alignment)) {
                    alignments.add(alignment);
                }
            }
            if (alignments.size() > 0) {
                alignmentGenerator.writeAlignments(alignments, writer);
                generateInsertionPoint(ComponentGenerator.InsertionPoint.AFTER_SEQUENCES, writer);
            }
        } catch (Exception e) {
            e.printStackTrace();
            throw new GeneratorException("Alignments generation has failed:\n" + e.getMessage());
        }

        //++++++++++++++++ Pattern Lists ++++++++++++++++++
        try {
            if (!options.samplePriorOnly) {
                List<Microsatellite> microsatList = new ArrayList<Microsatellite>();
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

                            case DataType.MICRO_SAT:
                                // microsat does not have alignment
                                patternListGenerator.writePatternList((PartitionPattern) partition, microsatList, writer);
                                break;
                            default:
                                throw new IllegalArgumentException("Unsupported data type");
                        }
                        writer.writeText("");
                    }
                }
            }
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
            if (taxonSets != null && taxonSets.size() > 0 && !options.useStarBEAST) {
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
            }

            for (PartitionTreePrior prior : options.getPartitionTreePriors()) {
                treePriorGenerator.writeEBSPVariableDemographic(prior, writer);
            }

            generateInsertionPoint(ComponentGenerator.InsertionPoint.AFTER_TREE_PRIOR, writer);
        } catch (Exception e) {
            e.printStackTrace();
            throw new GeneratorException("Tree prior likelihood generation has failed:\n" + e.getMessage());
        }

        //++++++++++++++++ Branch Rates Model ++++++++++++++++++
        try {
            for (PartitionClockModel model : options.getPartitionClockModels()) {
                branchRatesModelGenerator.writeBranchRatesModel(model, writer);
                writer.writeText("");
            }

            // write allClockRate for fix mean option in clock model panel
            for (ClockModelGroup clockModelGroup : options.clockModelOptions.getClockModelGroups()) {
                if (clockModelGroup.getRateTypeOption() == FixRateType.FIX_MEAN) {
                    writer.writeOpenTag(CompoundParameterParser.COMPOUND_PARAMETER,
                            new Attribute[]{new Attribute.Default<String>(XMLParser.ID, clockModelGroup.getName())});
                    for (PartitionClockModel model : options.getPartitionClockModels(clockModelGroup)) {
                        branchRatesModelGenerator.writeAllClockRateRefs(model, writer);
                    }
                    writer.writeCloseTag(CompoundParameterParser.COMPOUND_PARAMETER);
                    writer.writeText("");
                }
            }
        } catch (Exception e) {
            e.printStackTrace();
            throw new GeneratorException("Branch rates model generation is failed:\n" + e.getMessage());
        }

        //++++++++++++++++ Substitution Model & Site Model ++++++++++++++++++
        try {
            for (PartitionSubstitutionModel model : options.getPartitionSubstitutionModels()) {
                substitutionModelGenerator.writeSubstitutionSiteModel(model, writer);
                substitutionModelGenerator.writeAllMus(model, writer); // allMus
                writer.writeText("");
            }

            generateInsertionPoint(ComponentGenerator.InsertionPoint.AFTER_SUBSTITUTION_MODEL, writer);
        } catch (Exception e) {
            e.printStackTrace();
            throw new GeneratorException("Substitution model or site model generation has failed:\n" + e.getMessage());
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
            for (AbstractPartitionData partition : options.dataPartitions) {
                // generate tree likelihoods for alignment data partitions
                if (partition.getTaxonList() != null) {
                    if (partition instanceof PartitionData) {
                        if (partition.getDataType() != GeneralDataType.INSTANCE &&
                                partition.getDataType() != ContinuousDataType.INSTANCE) {
                            treeLikelihoodGenerator.writeTreeLikelihood((PartitionData) partition, writer);
                            writer.writeText("");
                        }
                    } else if (partition instanceof PartitionPattern) { // microsat
                        treeLikelihoodGenerator.writeTreeLikelihood((PartitionPattern) partition, writer);
                        writer.writeText("");
                    } else {
                        throw new GeneratorException("Find unrecognized partition:\n" + partition.getName());
                    }
                }
            }

            generateInsertionPoint(ComponentGenerator.InsertionPoint.AFTER_TREE_LIKELIHOOD, writer);
        } catch (Exception e) {
            e.printStackTrace();
            throw new GeneratorException("Tree likelihood generation has failed:\n" + e.getMessage());
        }

        //++++++++++++++++ *BEAST ++++++++++++++++++
        if (options.useStarBEAST) {
            //++++++++++++++++ species ++++++++++++++++++
            try {
                starBeastGenerator.writeSpecies(writer);

            } catch (Exception e) {
                e.printStackTrace();
                throw new GeneratorException("*BEAST species section generation has failed:\n" + e.getMessage());
            }

            //++++++++++++++++ Species Sets ++++++++++++++++++
            List<Taxa> speciesSets = options.speciesSets;
            try {
                if (speciesSets != null && speciesSets.size() > 0) {
                    tmrcaStatisticsGenerator.writeTaxonSets(writer, speciesSets);
                }
            } catch (Exception e) {
                e.printStackTrace();
                throw new GeneratorException("Species sets generation has failed:\n" + e.getMessage());
            }

            //++++++++++++++++ trees ++++++++++++++++++
            try {
                if (speciesSets != null && speciesSets.size() > 0) {
                    starBeastGenerator.writeStartingTreeForCalibration(writer);
                }

                starBeastGenerator.writeSpeciesTree(writer, speciesSets != null && speciesSets.size() > 0);
            } catch (Exception e) {
                e.printStackTrace();
                throw new GeneratorException("*BEAST trees generation has failed:\n" + e.getMessage());
            }

            //++++++++++++++++ Statistics ++++++++++++++++++
            try {
                if (speciesSets != null && speciesSets.size() > 0) {
                    tmrcaStatisticsGenerator.writeTMRCAStatistics(writer);
                }
            } catch (Exception e) {
                e.printStackTrace();
                throw new GeneratorException("*BEAST TMRCA statistics generation has failed:\n" + e.getMessage());
            }

            //++++++++++++++++ prior and likelihood ++++++++++++++++++
            try {
                starBeastGenerator.writeSTARBEAST(writer);
            } catch (Exception e) {
                e.printStackTrace();
                throw new GeneratorException("*BEAST trees section generation has failed:\n" + e.getMessage());
            }

        }
        generateInsertionPoint(ComponentGenerator.InsertionPoint.AFTER_TRAITS, writer);

        //++++++++++++++++ Operators ++++++++++++++++++
        try {
            List<Operator> operators = options.selectOperators();
            operatorsGenerator.writeOperatorSchedule(operators, writer);
            writer.writeText("");

            generateInsertionPoint(ComponentGenerator.InsertionPoint.AFTER_OPERATORS, writer);
        } catch (Exception e) {
            e.printStackTrace();
            throw new GeneratorException("Operators generation has failed:\n" + e.getMessage());
        }

        //++++++++++++++++ MCMC ++++++++++++++++++
        try {
            // XMLWriter writer, List<PartitionSubstitutionModel> models,
            writeMCMC(writer);
            writer.writeText("");

            generateInsertionPoint(ComponentGenerator.InsertionPoint.AFTER_MCMC, writer);
        } catch (Exception e) {
            e.printStackTrace();
            throw new GeneratorException("MCMC or log generation has failed:\n" + e.getMessage());
        }

        //++++++++++++++++  ++++++++++++++++++
        try {
            writeTimerReport(writer);
            writer.writeText("");
            if (options.performTraceAnalysis) {
                writeTraceAnalysis(writer);
            }
            if (options.generateCSV) {
                for (PartitionTreePrior prior : options.getPartitionTreePriors()) {
                    treePriorGenerator.writeEBSPAnalysisToCSVfile(prior, writer);
                }
            }
        } catch (Exception e) {
            e.printStackTrace();
            throw new GeneratorException("The last part of XML generation has failed:\n" + e.getMessage());
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

        // why was this here?
//        boolean hasAttr = (options.traits.size() > 0 && !options.useStarBEAST)
//                || (options.traits.size() > 1 && options.useStarBEAST);

        boolean firstDate = true;
        for (int i = 0; i < taxonList.getTaxonCount(); i++) {
            Taxon taxon = taxonList.getTaxon(i);

            boolean hasDate = false;
            if (options.clockModelOptions.isTipCalibrated()) {
                hasDate = TaxonList.Utils.hasAttribute(taxonList, i, dr.evolution.util.Date.DATE);
            }

            writer.writeTag(TaxonParser.TAXON, new Attribute[]{
                    new Attribute.Default<String>(XMLParser.ID, taxon.getId())},
                    !(hasDate || hasAttr)); // false if any of hasDate or hasAttr is true


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

                Attribute[] attributes = {
                        new Attribute.Default<Double>(DateParser.VALUE, date.getTimeValue()),
                        new Attribute.Default<String>(DateParser.DIRECTION, date.isBackwards() ? DateParser.BACKWARDS : DateParser.FORWARDS),
                        new Attribute.Default<String>(DateParser.UNITS, Units.Utils.getDefaultUnitName(options.units))
                        //new Attribute.Default("origin", date.getOrigin()+"")
                };

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

        writer.writeCloseTag(TaxaParser.TAXA);
    }

    public void writeDifferentTaxa(AbstractPartitionData dataPartition, XMLWriter writer) {
        TaxonList taxonList = dataPartition.getTaxonList();

        String name = dataPartition.getName();

        writer.writeComment("gene name = " + name + ", ntax= " + taxonList.getTaxonCount());
        writer.writeOpenTag(TaxaParser.TAXA, new Attribute[]{new Attribute.Default<String>(XMLParser.ID, name + "." + TaxaParser.TAXA)});

        for (int i = 0; i < taxonList.getTaxonCount(); i++) {
            final Taxon taxon = taxonList.getTaxon(i);
            writer.writeIDref(TaxonParser.TAXON, taxon.getId());
        }

        writer.writeCloseTag(TaxaParser.TAXA);
    }

    /**
     * Write the timer report block.
     *
     * @param writer the writer
     */
    public void writeTimerReport(XMLWriter writer) {
        writer.writeOpenTag("report");
        writer.writeOpenTag("property", new Attribute.Default<String>("name", "timer"));
        writer.writeIDref("mcmc", "mcmc");
        writer.writeCloseTag("property");
        writer.writeCloseTag("report");
    }

    /**
     * Write the trace analysis block.
     *
     * @param writer the writer
     */
    public void writeTraceAnalysis(XMLWriter writer) {
        writer.writeTag(
                "traceAnalysis",
                new Attribute[]{
                        new Attribute.Default<String>("fileName", options.logFileName)
                },
                true
        );
    }

    /**
     * Write the MCMC block.
     *
     * @param writer XMLWriter
     */
    public void writeMCMC(XMLWriter writer) {
        writer.writeComment("Define MCMC");

        List<Attribute> attributes = new ArrayList<Attribute>();
        attributes.add(new Attribute.Default<String>(XMLParser.ID, "mcmc"));
        attributes.add(new Attribute.Default<Integer>("chainLength", options.chainLength));
        attributes.add(new Attribute.Default<String>("autoOptimize", options.autoOptimize ? "true" : "false"));

        if (options.operatorAnalysis) {
            attributes.add(new Attribute.Default<String>("operatorAnalysis", options.operatorAnalysisFileName));
        }

        writer.writeOpenTag("mcmc", attributes);

        if (options.hasData()) {
            writer.writeOpenTag(CompoundLikelihoodParser.POSTERIOR, new Attribute.Default<String>(XMLParser.ID, "posterior"));
        }

        // write prior block
        writer.writeOpenTag(CompoundLikelihoodParser.PRIOR, new Attribute.Default<String>(XMLParser.ID, "prior"));

        if (options.useStarBEAST) { // species
            // coalescent prior
            writer.writeIDref(MultiSpeciesCoalescentParser.SPECIES_COALESCENT, TraitData.TRAIT_SPECIES + "." + COALESCENT);
            // prior on population sizes
//            if (options.speciesTreePrior == TreePriorType.SPECIES_YULE) {
            writer.writeIDref(MixedDistributionLikelihoodParser.DISTRIBUTION_LIKELIHOOD, SPOPS);
//            } else {
//                writer.writeIDref(SpeciesTreeBMPrior.STPRIOR, STP);
//            }
            // prior on species tree
            writer.writeIDref(SpeciationLikelihoodParser.SPECIATION_LIKELIHOOD, SPECIATION_LIKE);
        }

        parameterPriorGenerator.writeParameterPriors(writer, options.useStarBEAST);

        for (PartitionTreeModel model : options.getPartitionTreeModels()) {
            PartitionTreePrior prior = model.getPartitionTreePrior();
            treePriorGenerator.writePriorLikelihoodReference(prior, model, writer);
            writer.writeText("");
        }

        for (PartitionTreePrior prior : options.getPartitionTreePriors()) {
            treePriorGenerator.writeEBSPVariableDemographicReference(prior, writer);
            writer.writeText("");
        }

        generateInsertionPoint(ComponentGenerator.InsertionPoint.IN_MCMC_PRIOR, writer);

        writer.writeCloseTag(CompoundLikelihoodParser.PRIOR);

        if (options.hasData()) {
            // write likelihood block
            writer.writeOpenTag(CompoundLikelihoodParser.LIKELIHOOD, new Attribute.Default<String>(XMLParser.ID, "likelihood"));

            treeLikelihoodGenerator.writeTreeLikelihoodReferences(writer);
            branchRatesModelGenerator.writeClockLikelihoodReferences(writer);

            generateInsertionPoint(ComponentGenerator.InsertionPoint.IN_MCMC_LIKELIHOOD, writer);

            writer.writeCloseTag(CompoundLikelihoodParser.LIKELIHOOD);

            writer.writeCloseTag(CompoundLikelihoodParser.POSTERIOR);
        }

        writer.writeIDref(SimpleOperatorScheduleParser.OPERATOR_SCHEDULE, "operators");

        // write log to screen
        logGenerator.writeLogToScreen(writer, branchRatesModelGenerator, substitutionModelGenerator);

        // write log to file
        logGenerator.writeLogToFile(writer, treePriorGenerator, branchRatesModelGenerator,
                substitutionModelGenerator, treeLikelihoodGenerator);

        // write tree log to file
        logGenerator.writeTreeLogToFile(writer);

        writer.writeCloseTag("mcmc");
    }

}
