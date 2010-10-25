/*
 * BeastGenerator.java
 *
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

package dr.app.beauti.generator;

import dr.app.beast.BeastVersion;
import dr.app.beauti.components.ComponentFactory;
import dr.app.beauti.types.ClockType;
import dr.app.beauti.types.FixRateType;
import dr.app.beauti.types.StartingTreeType;
import dr.app.beauti.types.TreePriorType;
import dr.app.beauti.options.*;
import dr.app.beauti.util.XMLWriter;
import dr.app.util.Arguments;
import dr.evolution.alignment.Alignment;
import dr.evolution.alignment.SitePatterns;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.util.Taxa;
import dr.evolution.util.Taxon;
import dr.evolution.util.TaxonList;
import dr.evolution.util.Units;
import dr.evomodel.substmodel.AbstractSubstitutionModel;
import dr.evomodelxml.speciation.MultiSpeciesCoalescentParser;
import dr.evomodelxml.speciation.SpeciationLikelihoodParser;
import dr.evomodelxml.substmodel.GeneralSubstitutionModelParser;
import dr.evoxml.*;
import dr.inferencexml.distribution.MixedDistributionLikelihoodParser;
import dr.inferencexml.model.CompoundLikelihoodParser;
import dr.inferencexml.model.CompoundParameterParser;
import dr.inferencexml.operators.SimpleOperatorScheduleParser;
import dr.util.Attribute;
import dr.util.Version;
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

    private final AlignmentGenerator alignmentGenerator;
    private final TreePriorGenerator treePriorGenerator;
    private final TreeLikelihoodGenerator treeLikelihoodGenerator;
    private final SubstitutionModelGenerator substitutionModelGenerator;
    private final InitialTreeGenerator initialTreeGenerator;
    private final TreeModelGenerator treeModelGenerator;
    private final BranchRatesModelGenerator branchRatesModelGenerator;
    private final OperatorsGenerator operatorsGenerator;
    private final ParameterPriorGenerator parameterPriorGenerator;
    private final LogGenerator logGenerator;
    private final DiscreteTraitGenerator discreteTraitGenerator;
    private final STARBEASTGenerator starEASTGeneratorGenerator;
    private final TMRCAStatisticsGenerator tmrcaStatisticsGenerator;

    public BeastGenerator(BeautiOptions options, ComponentFactory[] components) {
        super(options, components);

        alignmentGenerator = new AlignmentGenerator(options, components);
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

        discreteTraitGenerator = new DiscreteTraitGenerator(options, components);
        starEASTGeneratorGenerator = new STARBEASTGenerator(options, components);
    }

    /**
     * Checks various options to check they are valid. Throws IllegalArgumentExceptions with
     * descriptions of the problems.
     *
     * @throws IllegalArgumentException if there is a problem with the current settings
     */
    public void checkOptions() throws IllegalArgumentException {
        //++++++++++++++++ Taxon List ++++++++++++++++++
        TaxonList taxonList = options.taxonList;
        Set<String> ids = new HashSet<String>();

        ids.add(TaxaParser.TAXA);
        ids.add(AlignmentParser.ALIGNMENT);

        if (taxonList != null) {
            if (taxonList.getTaxonCount() < 2) {
                throw new IllegalArgumentException("BEAST requires at least two taxa to run.");
            }

            for (int i = 0; i < taxonList.getTaxonCount(); i++) {
                Taxon taxon = taxonList.getTaxon(i);
                if (ids.contains(taxon.getId())) {
                    throw new IllegalArgumentException("A taxon has the same id," + taxon.getId() +
                            "\nas another element (taxon, sequence, taxon set etc.):\nAll ids should be unique.");
                }
                ids.add(taxon.getId());
            }
        }

        //++++++++++++++++ Taxon Sets ++++++++++++++++++
        for (Taxa taxa : options.taxonSets) {
            if (taxa.getTaxonCount() < 2) {
                throw new IllegalArgumentException("Taxon set, " + taxa.getId() + ", should contain\n" +
                        "at least two taxa.");
            }
            if (ids.contains(taxa.getId())) {
                throw new IllegalArgumentException("A taxon sets has the same id," + taxa.getId() +
                        "\nas another element (taxon, sequence, taxon set etc.):\nAll ids should be unique.");
            }
            ids.add(taxa.getId());
        }

        //++++++++++++++++ Tree Prior ++++++++++++++++++
//        if (options.isShareSameTreePrior()) {
        if (options.getPartitionTreeModels().size() > 1) { //TODO not allowed multi-prior yet
            for (PartitionTreePrior prior : options.getPartitionTreePriors()) {
                if (prior.getNodeHeightPrior() == TreePriorType.GMRF_SKYRIDE) {
                    throw new IllegalArgumentException("For GMRF, tree model/tree prior combination not implemented by BEAST yet" +
                            "\nIt is only available for single tree model partition for this release.");
                }
            }
        }

        //+++++++++++++++ Starting tree ++++++++++++++++
        for (PartitionTreeModel model : options.getPartitionTreeModels()) {
            if (model.getStartingTreeType() == StartingTreeType.USER) {
                if (model.getUserStartingTree() == null) {
                    throw new IllegalArgumentException("Please selected a starting tree in Trees panel");
                }
            }
        }

        //++++++++++++++++ Random local clock model validation ++++++++++++++++++
        for (PartitionClockModel model : options.getPartitionClockModels()) {
            // 1 random local clock CANNOT have different tree models
            if (model.getClockType() == ClockType.RANDOM_LOCAL_CLOCK) { // || AUTOCORRELATED_LOGNORMAL
                PartitionTreeModel treeModel = null;
                for (PartitionData pd : options.getAllPartitionData(model)) { // only the PDs linked to this tree model
                    if (treeModel != null && treeModel != pd.getPartitionTreeModel()) {
                        throw new IllegalArgumentException("One random local clock CANNOT have different tree models !");
                    }
                    treeModel = pd.getPartitionTreeModel();
                }
            }
        }

        //++++++++++++++++ Tree Model ++++++++++++++++++
        if (options.allowDifferentTaxa) {
            for (PartitionTreeModel model : options.getPartitionTreeModels()) {
                int numOfTaxa = -1;
                for (PartitionData pd : options.getAllPartitionData(model)) {
                    if (pd.getAlignment() != null) {
                        if (numOfTaxa > 0) {
                            if (numOfTaxa != pd.getTaxaCount()) {
                                throw new IllegalArgumentException("Partitions with different taxa cannot share the same tree");
                            }
                        } else {
                            numOfTaxa = pd.getTaxaCount();
                        }
                    }
                }
            }
        }

        //++++++++++++++++ User Modified Prior ++++++++++++++++++
        for (Parameter param : options.selectParameters()) {
            if (param.isPriorEdited() && param.initial != Double.NaN) {
                if (param.initial < param.lower || param.initial > param.upper) {
                    throw new IllegalArgumentException("Parameter \"" + param.getName() + "\":" +
                            "\ninitial value " + param.initial + " is NOT in the range [" + param.lower + ", " + param.upper + "]," +
                            "\nor this range is wrong. Please check the Prior panel.");
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
                "      by Alexei J. Drummond and Andrew Rambaut",
                "      Department of Computer Science, University of Auckland and",
                "      Institute of Evolutionary Biology, University of Edinburgh",
                "      http://beast.bio.ed.ac.uk/");
        writer.writeOpenTag("beast");
        writer.writeText("");

        // this gives any added implementations of the 'Component' interface a
        // chance to generate XML at this point in the BEAST file.
        generateInsertionPoint(ComponentGenerator.InsertionPoint.BEFORE_TAXA, writer);

        //++++++++++++++++ Taxon List ++++++++++++++++++
        try {
            writeTaxa(writer, options.taxonList);

            if (options.allowDifferentTaxa) { // allow diff taxa for multi-gene
                writer.writeText("");
                writer.writeComment("List all taxons regarding each gene (file) for Multispecies Coalescent function");
                // write all taxa in each gene tree regarding each data partition,
                for (PartitionData partition : options.dataPartitions) {
                    // do I need if (!alignments.contains(alignment)) {alignments.add(alignment);} ?
                    if (partition.getAlignment() != null) {
                        writeDifferentTaxaForMultiGene(partition, writer);
                    }
                }
            }
        } catch (Exception e) {
            System.err.println(e);
            throw new GeneratorException("Taxon list generation has failed");
        }

        //++++++++++++++++ Taxon Sets ++++++++++++++++++
        List<Taxa> taxonSets = options.taxonSets;
        try {
            if (taxonSets != null && taxonSets.size() > 0) {
                tmrcaStatisticsGenerator.writeTaxonSets(writer, taxonSets);
            }
            generateInsertionPoint(ComponentGenerator.InsertionPoint.AFTER_TAXA, writer);
        } catch (Exception e) {
            e.printStackTrace();
            throw new GeneratorException("Taxon sets generation has failed");
        }

        //++++++++++++++++ Alignments ++++++++++++++++++
        try {
            alignmentGenerator.writeAlignments(writer);
            generateInsertionPoint(ComponentGenerator.InsertionPoint.AFTER_SEQUENCES, writer);
        } catch (Exception e) {
            e.printStackTrace();
            throw new GeneratorException("Alignments generation has failed");
        }

        //++++++++++++++++ Pattern Lists ++++++++++++++++++
        try {
            if (!options.samplePriorOnly) {
                for (PartitionData partition : options.dataPartitions) { // Each PD has one TreeLikelihood
                    if (partition.getAlignment() != null) {
                        writePatternList(partition, writer);
                        writer.writeText("");
                    }
                }
            }
        } catch (Exception e) {
            e.printStackTrace();
            throw new GeneratorException("Pattern lists generation has failed");
        }

        //++++++++++++++++ General Data of Traits ++++++++++++++++++
        try {
            for (PartitionData partition : options.dataPartitions) { // Each PD has one TreeLikelihood
                TraitData trait = partition.getTrait();
                if (trait != null) {
                    discreteTraitGenerator.writeGeneralDataType(trait, writer);
                    writer.writeText("");
                }
            }

            generateInsertionPoint(ComponentGenerator.InsertionPoint.AFTER_PATTERNS, writer);
        } catch (Exception e) {
            e.printStackTrace();
            throw new GeneratorException("General data of traits generation has failed");
        }

        //++++++++++++++++ Tree Prior Model ++++++++++++++++++
        try {
            for (PartitionTreePrior prior : options.getPartitionTreePriors()) {
                treePriorGenerator.writeTreePriorModel(prior, writer);
                writer.writeText("");
            }
        } catch (Exception e) {
            e.printStackTrace();
            throw new GeneratorException("Tree prior model generation has failed");
        }

        //++++++++++++++++ Starting Tree ++++++++++++++++++
        try {
            for (PartitionTreeModel model : options.getPartitionTreeModels()) {
                initialTreeGenerator.writeStartingTree(model, writer);
                writer.writeText("");
            }
        } catch (Exception e) {
            e.printStackTrace();
            throw new GeneratorException("Starting tree generation has failed");
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
            throw new GeneratorException("Tree model generation has failed");
        }

        //++++++++++++++++ Tree Prior Likelihood ++++++++++++++++++
        try {
            for (PartitionTreeModel model : options.getPartitionTreeModels()) {
                PartitionTreePrior prior = model.getPartitionTreePrior();
                treePriorGenerator.writePriorLikelihood(prior, model, writer);
                writer.writeText("");
            }

            for (PartitionTreePrior prior : options.getPartitionTreePriors()) {
                treePriorGenerator.writeEBSPVariableDemographic(prior, writer);
            }

            generateInsertionPoint(ComponentGenerator.InsertionPoint.AFTER_TREE_PRIOR, writer);
        } catch (Exception e) {
            e.printStackTrace();
            throw new GeneratorException("Tree prior likelihood generation has failed");
        }

        //++++++++++++++++ Branch Rates Model ++++++++++++++++++
        try {
            for (PartitionClockModel model : options.getPartitionClockModels()) {
                branchRatesModelGenerator.writeBranchRatesModel(model, writer);
                writer.writeText("");
            }

            // write allClockRate for fix mean option in clock model panel
            if (options.clockModelOptions.getRateOptionClockModel() == FixRateType.FIX_MEAN) {
                writer.writeOpenTag(CompoundParameterParser.COMPOUND_PARAMETER, new Attribute[]{new Attribute.Default<String>(XMLParser.ID, "allClockRates")});
                for (PartitionClockModel model : options.getPartitionClockModels()) {
                    branchRatesModelGenerator.writeAllClockRateRefs(model, writer);
                }
                writer.writeCloseTag(CompoundParameterParser.COMPOUND_PARAMETER);
                writer.writeText("");
            }
        } catch (Exception e) {
            e.printStackTrace();
            throw new GeneratorException("Branch rates model generation is failed");
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
            throw new GeneratorException("Substitution model or site model generation has failed");
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
            // generate tree likelihoods for alignment data partitions
            if (options.hasAlignmentPartition()) {
                writer.writeComment("Likelihood for tree given sequence data");
            }
            for (PartitionData partition : options.dataPartitions) {
                if (partition.getAlignment() != null) {
                    treeLikelihoodGenerator.writeTreeLikelihood(partition, writer);
                    writer.writeText("");
                }
            }

            generateInsertionPoint(ComponentGenerator.InsertionPoint.AFTER_TREE_LIKELIHOOD, writer);
        } catch (Exception e) {
            e.printStackTrace();
            throw new GeneratorException("Tree likelihood generation has failed");
        }

        //++++++++++++++++ Ancestral Tree Likelihood ++++++++++++++++++
        try {
            // generate tree likelihoods for discrete trait partitions
            if (options.hasDiscreteTraitPartition()) {
                writer.writeComment("Likelihood for tree given discrete trait data");
            }
            for (PartitionData partition : options.dataPartitions) {
                TraitData trait = partition.getTrait();
                if (trait != null && trait.getTraitType() == TraitData.TraitType.DISCRETE) {
                    discreteTraitGenerator.writeAncestralTreeLikelihood(partition, writer);
                }
            }
        } catch (Exception e) {
            e.printStackTrace();
            throw new GeneratorException("Ancestral Tree likelihood generation has failed");
        }

        //++++++++++++++++ Multivariate Diffusion Tree Likelihood ++++++++++++++++++
        try {
            // generate tree likelihoods for continuous trait partitions
            if (options.hasContinuousTraitPartition()) {
                writer.writeComment("Likelihood for tree given continuous multivariate trait data");
            }
            for (PartitionData partition : options.dataPartitions) {
                TraitData trait = partition.getTrait();
                if (trait != null && trait.getTraitType() == TraitData.TraitType.CONTINUOUS) {
                    throw new UnsupportedOperationException("Not implemented yet: writeMultivariateTreeLikelihood");
//                    generalTraitGenerator.writeMultivariateTreeLikelihood(trait, writer);
                }
            }
        } catch (Exception e) {
            e.printStackTrace();
            throw new GeneratorException("Ancestral Tree likelihood generation has failed");
        }

        //++++++++++++++++ *BEAST ++++++++++++++++++
        try {
            if (options.useStarBEAST) { // species
                writeStarBEAST(writer);
            }
        } catch (Exception e) {
            e.printStackTrace();
            throw new GeneratorException("*BEAST special part generation has failed");
        }

        generateInsertionPoint(ComponentGenerator.InsertionPoint.AFTER_TRAITS, writer);

        //++++++++++++++++ Statistics ++++++++++++++++++
        try {
            if (taxonSets != null && taxonSets.size() > 0) {
                tmrcaStatisticsGenerator.writeTMRCAStatistics(writer);
            }
        } catch (Exception e) {
            e.printStackTrace();
            throw new GeneratorException("TMRCA statistics generation has failed");
        }

        //++++++++++++++++ Operators ++++++++++++++++++
        try {
            List<Operator> operators = options.selectOperators();
            operatorsGenerator.writeOperatorSchedule(operators, writer);
            writer.writeText("");

            generateInsertionPoint(ComponentGenerator.InsertionPoint.AFTER_OPERATORS, writer);
        } catch (Exception e) {
            e.printStackTrace();
            throw new GeneratorException("Operators generation has failed");
        }

        //++++++++++++++++ MCMC ++++++++++++++++++
        try {
            // XMLWriter writer, List<PartitionSubstitutionModel> models,
            writeMCMC(writer);
            writer.writeText("");

            generateInsertionPoint(ComponentGenerator.InsertionPoint.AFTER_MCMC, writer);
        } catch (Exception e) {
            e.printStackTrace();
            throw new GeneratorException("MCMC or log generation is failed");
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
            throw new GeneratorException("The last part of XML generation has failed");
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
    private void writeTaxa(XMLWriter writer, TaxonList taxonList) throws Arguments.ArgumentException {
        // -1 (single taxa), 0 (1st gene of multi-taxa)

        writer.writeComment("The list of taxa analyse (can also include dates/ages).",
                "ntax=" + taxonList.getTaxonCount());
        writer.writeOpenTag(TaxaParser.TAXA, new Attribute[]{new Attribute.Default<String>(XMLParser.ID, TaxaParser.TAXA)});

        boolean hasAttr = options.traits.size() > 0;

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

            if (hasAttr) {
                discreteTraitGenerator.writeTaxonTraits(taxon, writer);
            }

            if (hasDate || hasAttr) writer.writeCloseTag(TaxonParser.TAXON);


            generateInsertionPoint(ComponentGenerator.InsertionPoint.IN_TAXON, taxon, writer);
        }

        writer.writeCloseTag(TaxaParser.TAXA);
    }

    public void writeDifferentTaxaForMultiGene(PartitionData dataPartition, XMLWriter writer) {
        Alignment alignment = dataPartition.getAlignment();

        String data = dataPartition.getName();

        writer.writeComment("gene name = " + data + ", ntax= " + alignment.getTaxonCount());
        writer.writeOpenTag(TaxaParser.TAXA, new Attribute[]{new Attribute.Default<String>(XMLParser.ID, data + "." + TaxaParser.TAXA)});

        for (int i = 0; i < alignment.getTaxonCount(); i++) {
            final Taxon taxon = alignment.getTaxon(i);
            writer.writeIDref(TaxonParser.TAXON, taxon.getId());
        }

        writer.writeCloseTag(TaxaParser.TAXA);
    }

    /**
     * *BEAST block
     *
     * @param writer XMLWriter
     */
    private void writeStarBEAST(XMLWriter writer) {
        String traitName = TraitData.TRAIT_SPECIES;
        writer.writeText("");
        writer.writeComment(options.starBEASTOptions.getDescription());

        writer.writeOpenTag(traitName, new Attribute[]{
                new Attribute.Default<String>(XMLParser.ID, traitName)});
        //new Attribute.Default<String>("traitType", traitType)});
        starEASTGeneratorGenerator.writeMultiSpecies(options.taxonList, writer);
        writer.writeCloseTag(traitName);

        starEASTGeneratorGenerator.writeSTARBEAST(writer);
    }

    /**
     * Writes the pattern lists
     *
     * @param partition the partition data to write the pattern lists for
     * @param writer    the writer
     */
    public void writePatternList(PartitionData partition, XMLWriter writer) {
        writer.writeText("");

        PartitionSubstitutionModel model = partition.getPartitionSubstitutionModel();

        String codonHeteroPattern = model.getCodonHeteroPattern();
        int partitionCount = model.getCodonPartitionCount();

        if (model.getDataType() == Nucleotides.INSTANCE && codonHeteroPattern != null && partitionCount > 1) {

            if (codonHeteroPattern.equals("112")) {
                writer.writeComment("The unique patterns for codon positions 1 & 2");
                writer.writeOpenTag(MergePatternsParser.MERGE_PATTERNS,
                        new Attribute[]{
                                new Attribute.Default<String>(XMLParser.ID, model.getPrefix(1) + partition.getPrefix() + SitePatternsParser.PATTERNS),
                        }
                );
                writePatternList(partition, 0, 3, writer);
                writePatternList(partition, 1, 3, writer);

                writer.writeCloseTag(MergePatternsParser.MERGE_PATTERNS);

                writer.writeComment("The unique patterns for codon positions 3");
                writer.writeOpenTag(MergePatternsParser.MERGE_PATTERNS,
                        new Attribute[]{
                                new Attribute.Default<String>(XMLParser.ID, model.getPrefix(2) + partition.getPrefix() + SitePatternsParser.PATTERNS),
                        }
                );

                writePatternList(partition, 2, 3, writer);

                writer.writeCloseTag(MergePatternsParser.MERGE_PATTERNS);

            } else {
                // pattern is 123
                // write pattern lists for all three codon positions
                for (int i = 1; i <= 3; i++) {
                    writer.writeComment("The unique patterns for codon positions " + i);
                    writer.writeOpenTag(MergePatternsParser.MERGE_PATTERNS,
                            new Attribute[]{
                                    new Attribute.Default<String>(XMLParser.ID, model.getPrefix(i) + partition.getPrefix() + SitePatternsParser.PATTERNS),
                            }
                    );

                    writePatternList(partition, i - 1, 3, writer);

                    writer.writeCloseTag(MergePatternsParser.MERGE_PATTERNS);
                }

            }
        } else {
            writePatternList(partition, 0, 1, writer);
        }
    }

    /**
     * Write a single pattern list
     *
     * @param partition the partition to write a pattern list for
     * @param offset    offset by
     * @param every     skip every
     * @param writer    the writer
     */
    private void writePatternList(PartitionData partition, int offset, int every, XMLWriter writer) {

        Alignment alignment = partition.getAlignment();
        int from = partition.getFromSite();
        int to = partition.getToSite();
        int partEvery = partition.getEvery();
        if (partEvery > 1 && every > 1) throw new IllegalArgumentException();

        if (from < 1) from = 1;
        every = Math.max(partEvery, every);

        from += offset;


        // this object is created solely to calculate the number of patterns in the alignment
        SitePatterns patterns = new SitePatterns(alignment, from - 1, to - 1, every);

        writer.writeComment("The unique patterns from " + from + " to " + (to > 0 ? to : "end") + ((every > 1) ? " every " + every : ""),
                "npatterns=" + patterns.getPatternCount());

        List<Attribute> attributes = new ArrayList<Attribute>();

        // no codon, unique patterns site patterns
        if (offset == 0 && every == 1)
            attributes.add(new Attribute.Default<String>(XMLParser.ID, partition.getPrefix() + SitePatternsParser.PATTERNS));

        attributes.add(new Attribute.Default<String>("from", "" + from));
        if (to >= 0) attributes.add(new Attribute.Default<String>("to", "" + to));

        if (every > 1) {
            attributes.add(new Attribute.Default<String>("every", "" + every));
        }

        // generate <patterns>
        writer.writeOpenTag(SitePatternsParser.PATTERNS, attributes);
        writer.writeIDref(AlignmentParser.ALIGNMENT, alignment.getId());
        writer.writeCloseTag(SitePatternsParser.PATTERNS);
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

        parameterPriorGenerator.writeParameterPriors(writer);

        for (PartitionTreeModel model : options.getPartitionTreeModels()) {
            PartitionTreePrior prior = model.getPartitionTreePrior();
            treePriorGenerator.writePriorLikelihoodReference(prior, model, writer);
            writer.writeText("");
        }

        for (PartitionTreePrior prior : options.getPartitionTreePriors()) {
            treePriorGenerator.writeEBSPVariableDemographicReference(prior, writer);
            writer.writeText("");
        }

        for (PartitionSubstitutionModel model : options.getPartitionSubstitutionModels()) {
            // e.g. <svsGeneralSubstitutionModel idref="locations.model" />
//            if (!(model.getLocationSubstType() == LocationSubstModelType.SYM_SUBST && (!model.isActivateBSSVS()))) {
            if (model.isActivateBSSVS()) {
                writer.writeIDref(GeneralSubstitutionModelParser.GENERAL_SUBSTITUTION_MODEL, model.getPrefix() + AbstractSubstitutionModel.MODEL);
                writer.writeText("");
            }
        }

        generateInsertionPoint(ComponentGenerator.InsertionPoint.IN_MCMC_PRIOR, writer);

        writer.writeCloseTag(CompoundLikelihoodParser.PRIOR);

        if (options.hasData()) {
            // write likelihood block
            writer.writeOpenTag(CompoundLikelihoodParser.LIKELIHOOD, new Attribute.Default<String>(XMLParser.ID, "likelihood"));

            treeLikelihoodGenerator.writeTreeLikelihoodReferences(writer);
            discreteTraitGenerator.writeAncestralTreeLikelihoodReferences(writer);
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
                substitutionModelGenerator, treeLikelihoodGenerator, discreteTraitGenerator);

        if (options.hasDiscreteTraitPartition()) {
            writer.writeComment("write discrete trait rate information to log");
        }

        for (PartitionData partition : options.dataPartitions) {
            if (partition.getTrait() != null && partition.getTrait().getTraitType() == TraitData.TraitType.DISCRETE) {
                logGenerator.writeDiscreteTraitLogToFile(writer, partition, branchRatesModelGenerator, substitutionModelGenerator);
            }
        }

        // write tree log to file
        logGenerator.writeTreeLogToFile(writer);

        writer.writeCloseTag("mcmc");
    }

}
