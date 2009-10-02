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
import dr.app.beauti.enumTypes.FixRateType;
import dr.app.beauti.enumTypes.TreePriorType;
import dr.app.beauti.options.*;
import dr.app.beauti.util.XMLWriter;
import dr.evolution.alignment.Alignment;
import dr.evolution.alignment.SitePatterns;
import dr.evolution.datatype.DataType;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.util.Taxa;
import dr.evolution.util.Taxon;
import dr.evolution.util.TaxonList;
import dr.evolution.util.Units;
import dr.evomodel.speciation.SpeciationLikelihood;
import dr.evomodel.speciation.TreePartitionCoalescent;
import dr.evomodel.tree.MonophylyStatistic;
import dr.evomodel.tree.TMRCAStatistic;
import dr.evomodel.tree.TreeModel;
import dr.evoxml.*;
import dr.inference.distribution.MixedDistributionLikelihood;
import dr.inference.model.*;
import dr.inference.operators.SimpleOperatorSchedule;
import dr.util.Attribute;
import dr.util.Version;
import dr.xml.XMLParser;

import java.io.Writer;
import java.util.*;

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

    private final TreePriorGenerator treePriorGenerator;
    private final TreeLikelihoodGenerator treeLikelihoodGenerator;
    private final SubstitutionModelGenerator substitutionModelGenerator;
    private final InitialTreeGenerator initialTreeGenerator;
    private final TreeModelGenerator treeModelGenerator;
    private final BranchRatesModelGenerator branchRatesModelGenerator;
    private final OperatorsGenerator operatorsGenerator;
    private final ParameterPriorGenerator parameterPriorGenerator;
    private final LogGenerator logGenerator;
    private final STARBEASTGenerator starEASTGeneratorGenerator;

    public BeastGenerator(BeautiOptions options, ComponentFactory[] components) {
        super(options, components);

        substitutionModelGenerator = new SubstitutionModelGenerator(options, components);
        treePriorGenerator = new TreePriorGenerator(options, components);
        treeLikelihoodGenerator = new TreeLikelihoodGenerator(options, components);

        initialTreeGenerator = new InitialTreeGenerator(options, components);
        treeModelGenerator = new TreeModelGenerator(options, components);
        branchRatesModelGenerator = new BranchRatesModelGenerator(options, components);

        operatorsGenerator = new OperatorsGenerator(options, components);
        parameterPriorGenerator = new ParameterPriorGenerator(options, components);
        logGenerator = new LogGenerator(options, components);

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
                    throw new IllegalArgumentException("For GMRF, tree model/tree prior combination not implemented by BEAST yet!" +
                            "\nIt is only available for single tree model partition for this release.");
                }
            }
        }

        //++++++++++++++++ clock model/tree model combination ++++++++++++++++++
        for (PartitionTreeModel model : options.getPartitionTreeModels()) {
            // clock model/tree model combination not implemented by BEAST yet
            validateClockTreeModelCombination(model);
        }


        //++++++++++++++++ Species tree ++++++++++++++++++
        if (options.starBEASTOptions.isSpeciesAnalysis()) {
//        	if (!(options.nodeHeightPrior == TreePriorType.SPECIES_BIRTH_DEATH || options.nodeHeightPrior == TreePriorType.SPECIES_YULE)) {
//        		//TODO: more species tree model
//        		throw new IllegalArgumentException("Species analysis requires to define species tree prior in Tree panel.");
//        	}
        }

        // add other tests and warnings here
        // Speciation model with dated tips
        // Sampling rates without dated tips or priors on rate or nodes

    }

    /**
     * Generate a beast xml file from these beast options
     *
     * @param w the writer
     */
    public void generateXML(Writer w) {

        XMLWriter writer = new XMLWriter(w);

        writer.writeText("<?xml version=\"1.0\" standalone=\"yes\"?>");
        writer.writeComment("Generated by BEAUTi " + version.getVersionString());
        writer.writeComment("      by Alexei J. Drummond and Andrew Rambaut");
        writer.writeComment("      Department of Computer Science, University of Auckland and");
        writer.writeComment("      Institute of Evolutionary Biology, University of Edinburgh");
        writer.writeComment("      http://beast.bio.ed.ac.uk/");
        writer.writeOpenTag("beast");
        writer.writeText("");

        // this gives any added implementations of the 'Component' interface a
        // chance to generate XML at this point in the BEAST file.
        generateInsertionPoint(ComponentGenerator.InsertionPoint.BEFORE_TAXA, writer);

        //++++++++++++++++ Taxon List ++++++++++++++++++
        writeTaxa(writer, options.taxonList);

        List<Taxa> taxonSets = options.taxonSets;
        if (taxonSets != null && taxonSets.size() > 0) {
            writeTaxonSets(writer, taxonSets); // TODO
        }

        if (options.allowDifferentTaxa) { // allow diff taxa for multi-gene
            writer.writeText("");
            writer.writeComment("List all taxons regarding each gene (file) for Multispecies Coalescent function");
            // write all taxa in each gene tree regarding each data partition,
            for (PartitionData partition : options.dataPartitions) {
                // do I need if (!alignments.contains(alignment)) {alignments.add(alignment);} ?
                writeDifferentTaxaForMultiGene(partition, writer);
            }
        }

        generateInsertionPoint(ComponentGenerator.InsertionPoint.AFTER_TAXA, writer);

        //++++++++++++++++ Alignments ++++++++++++++++++
        List<Alignment> alignments = new ArrayList<Alignment>();

        for (PartitionData partition : options.dataPartitions) {
            Alignment alignment = partition.getAlignment();
            if (!alignments.contains(alignment)) {
                alignments.add(alignment);
            }
        }

        if (!options.samplePriorOnly) {
            int index = 1;
            for (Alignment alignment : alignments) {
                if (alignments.size() > 1) {
                    //if (!options.allowDifferentTaxa) {
                    alignment.setId(AlignmentParser.ALIGNMENT + index);
                    //} else { // e.g. alignment_gene1
                    // alignment.setId("alignment_" + mulitTaxaTagName + index);
                    //}
                } else {
                    alignment.setId(AlignmentParser.ALIGNMENT);
                }
                writeAlignment(alignment, writer);
                index += 1;
                writer.writeText("");
            }

            generateInsertionPoint(ComponentGenerator.InsertionPoint.AFTER_SEQUENCES, writer);

            //++++++++++++++++ Pattern Lists ++++++++++++++++++
//            for (PartitionSubstitutionModel model : options.getPartitionSubstitutionModels()) {
//                writePatternList(model, writer);
            for (PartitionData partition : options.dataPartitions) { // Each PD has one TreeLikelihood
                writePatternList(partition, writer);
                writer.writeText("");
            }

            generateInsertionPoint(ComponentGenerator.InsertionPoint.AFTER_PATTERNS, writer);
        } else {
            Alignment alignment = alignments.get(0);
            alignment.setId(AlignmentParser.ALIGNMENT);
            writeAlignment(alignment, writer);
            writer.writeText("");

            generateInsertionPoint(ComponentGenerator.InsertionPoint.AFTER_SEQUENCES, writer);
            generateInsertionPoint(ComponentGenerator.InsertionPoint.AFTER_PATTERNS, writer);
        }

        //++++++++++++++++ Tree Prior Model ++++++++++++++++++
        for (PartitionTreePrior prior : options.getPartitionTreePriors()) {
            treePriorGenerator.writeTreePriorModel(prior, writer);
            writer.writeText("");
        }

        //++++++++++++++++ Starting Tree ++++++++++++++++++
        for (PartitionTreeModel model : options.getPartitionTreeModels()) {
            initialTreeGenerator.writeStartingTree(model, writer);
            writer.writeText("");
        }

        //++++++++++++++++ Tree Model +++++++++++++++++++

        for (PartitionTreeModel model : options.getPartitionTreeModels()) {
            treeModelGenerator.writeTreeModel(model, writer);
            writer.writeText("");
        }

        generateInsertionPoint(ComponentGenerator.InsertionPoint.AFTER_TREE_MODEL, writer);

        //++++++++++++++++ Tree Prior Likelihood ++++++++++++++++++
        for (PartitionTreeModel model : options.getPartitionTreeModels()) {
            PartitionTreePrior prior = model.getPartitionTreePrior();
            treePriorGenerator.writePriorLikelihood(prior, model, writer);
            writer.writeText("");
        }

        for (PartitionTreePrior prior : options.getPartitionTreePriors()) {
            treePriorGenerator.writeEBSPVariableDemographic(prior, writer);
        }

        generateInsertionPoint(ComponentGenerator.InsertionPoint.AFTER_TREE_PRIOR, writer);

        //++++++++++++++++ Branch Rates Model ++++++++++++++++++                
        for (PartitionClockModel model : options.getPartitionClockModels()) {
            branchRatesModelGenerator.writeBranchRatesModel(model, writer);
            writer.writeText("");
        }

        // write allClockRate for fix mean option in clock model panel
        if (options.clockModelOptions.getRateOptionClockModel() == FixRateType.FIX_MEAN) {
            writer.writeOpenTag(CompoundParameter.COMPOUND_PARAMETER, new Attribute[]{new Attribute.Default<String>(XMLParser.ID, "allClockRates")});
            for (PartitionClockModel model : options.getPartitionClockModels()) {
                branchRatesModelGenerator.writeAllClockRateRefs(model, writer);
            }
            writer.writeCloseTag(CompoundParameter.COMPOUND_PARAMETER);
            writer.writeText("");
        }

        //++++++++++++++++ Substitution Model ++++++++++++++++++
        for (PartitionSubstitutionModel model : options.getPartitionSubstitutionModels()) {
            substitutionModelGenerator.writeSubstitutionModel(model, writer);
            writer.writeText("");
        }

        generateInsertionPoint(ComponentGenerator.InsertionPoint.AFTER_SUBSTITUTION_MODEL, writer);

        //++++++++++++++++ Site Model ++++++++++++++++++
        for (PartitionSubstitutionModel model : options.getPartitionSubstitutionModels()) {            
            substitutionModelGenerator.writeSiteModel(model, writer); // site model
            substitutionModelGenerator.writeAllMus(model, writer); // allMus
            writer.writeText("");
        }

        generateInsertionPoint(ComponentGenerator.InsertionPoint.AFTER_SITE_MODEL, writer);

        //++++++++++++++++ Tree Likelihood ++++++++++++++++++
        for (PartitionData partition : options.dataPartitions) { // Each PD has one TreeLikelihood
            treeLikelihoodGenerator.writeTreeLikelihood(partition, writer);
            writer.writeText("");
        }

        generateInsertionPoint(ComponentGenerator.InsertionPoint.AFTER_TREE_LIKELIHOOD, writer);

        //++++++++++++++++ Traits ++++++++++++++++++
        // traits tag
        if (options.traitOptions.traits.size() > 0) {
            for (TraitGuesser trait : options.traitOptions.traits) {
                writeTraits(writer, trait, options.taxonList);
            }
            generateInsertionPoint(ComponentGenerator.InsertionPoint.AFTER_TRAITS, writer);
        }

        //++++++++++++++++  ++++++++++++++++++
        if (taxonSets != null && taxonSets.size() > 0) {
            //TODO: need to suit for multi-gene-tree
            writeTMRCAStatistics(writer);
        }

        //++++++++++++++++ Operators ++++++++++++++++++
        List<Operator> operators = options.selectOperators();
        operatorsGenerator.writeOperatorSchedule(operators, writer);
        writer.writeText("");

        generateInsertionPoint(ComponentGenerator.InsertionPoint.AFTER_OPERATORS, writer);

        //++++++++++++++++ MCMC ++++++++++++++++++
        // XMLWriter writer, List<PartitionSubstitutionModel> models,
        writeMCMC(writer);
        writer.writeText("");

        generateInsertionPoint(ComponentGenerator.InsertionPoint.AFTER_MCMC, writer);

        //++++++++++++++++  ++++++++++++++++++
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

        writer.writeCloseTag("beast");
        writer.flush();
    }

    /**
     * Generate a taxa block from these beast options
     *
     * @param writer    the writer
     * @param taxonList the taxon list to write
     */
    private void writeTaxa(XMLWriter writer, TaxonList taxonList) {
        // -1 (single taxa), 0 (1st gene of multi-taxa)

        writer.writeComment("The list of taxa analyse (can also include dates/ages).");
        writer.writeComment("ntax=" + taxonList.getTaxonCount());
        writer.writeOpenTag(TaxaParser.TAXA, new Attribute[]{new Attribute.Default<String>(XMLParser.ID, TaxaParser.TAXA)});


        boolean firstDate = true;
        for (int i = 0; i < taxonList.getTaxonCount(); i++) {
            Taxon taxon = taxonList.getTaxon(i);

            boolean hasDate = false;

            if (options.clockModelOptions.isTipCalibrated()) {
                hasDate = TaxonList.Utils.hasAttribute(taxonList, i, dr.evolution.util.Date.DATE);
            }


            writer.writeTag(TaxonParser.TAXON, new Attribute[]{new Attribute.Default<String>(XMLParser.ID, taxon.getId())}, !hasDate);


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
                writer.writeCloseTag(TaxonParser.TAXON);
            }

            generateInsertionPoint(ComponentGenerator.InsertionPoint.IN_TAXON, taxon, writer);
        }

        writer.writeCloseTag(TaxaParser.TAXA);
    }

    /**
     * Generate additional taxon sets
     *
     * @param writer    the writer
     * @param taxonSets a list of taxa to write
     */
    private void writeTaxonSets(XMLWriter writer, List<Taxa> taxonSets) {
        writer.writeText("");

        for (Taxa taxa : taxonSets) {
            writer.writeOpenTag(
                    TaxaParser.TAXA,
                    new Attribute[]{
                            new Attribute.Default<String>(XMLParser.ID, taxa.getId())
                    }
            );

            for (int j = 0; j < taxa.getTaxonCount(); j++) {
                writer.writeIDref(TaxonParser.TAXON, taxa.getTaxon(j).getId());
            }
            writer.writeCloseTag(TaxaParser.TAXA);
        }
    }


    /**
     * Determine and return the datatype description for these beast options
     * note that the datatype in XML may differ from the actual datatype
     *
     * @param alignment the alignment to get data type description of
     * @return description
     */

    private String getAlignmentDataTypeDescription(Alignment alignment) {
        String description;

        switch (alignment.getDataType().getType()) {
            case DataType.TWO_STATES:
            case DataType.COVARION:

                // TODO make this work
//                throw new RuntimeException("TO DO!");

                //switch (partition.getPartitionSubstitutionModel().binarySubstitutionModel) {
                //    case ModelOptions.BIN_COVARION:
                //        description = TwoStateCovarion.DESCRIPTION;
                //        break;
                //
                //    default:
                description = alignment.getDataType().getDescription();
                //}
                break;

            default:
                description = alignment.getDataType().getDescription();
        }

        return description;
    }


    public void writeDifferentTaxaForMultiGene(PartitionData dataPartition, XMLWriter writer) {
        String data = dataPartition.getName();
        Alignment alignment = dataPartition.getAlignment();

        writer.writeComment("gene name = " + data + ", ntax= " + alignment.getTaxonCount());
        writer.writeOpenTag(TaxaParser.TAXA, new Attribute[]{new Attribute.Default<String>(XMLParser.ID, data + "." + TaxaParser.TAXA)});

        for (int i = 0; i < alignment.getTaxonCount(); i++) {
            final Taxon taxon = alignment.getTaxon(i);
            writer.writeIDref(TaxonParser.TAXON, taxon.getId());
        }

        writer.writeCloseTag(TaxaParser.TAXA);
    }


    /**
     * Generate an alignment block from these beast options
     *
     * @param alignment the alignment to write
     * @param writer    the writer
     */
    public void writeAlignment(Alignment alignment, XMLWriter writer) {

        writer.writeText("");
        writer.writeComment("The sequence alignment (each sequence refers to a taxon above).");
        writer.writeComment("ntax=" + alignment.getTaxonCount() + " nchar=" + alignment.getSiteCount());
        if (options.samplePriorOnly) {
            writer.writeComment("Null sequences generated in order to sample from the prior only.");
        }

        writer.writeOpenTag(
                AlignmentParser.ALIGNMENT,
                new Attribute[]{
                        new Attribute.Default<String>(XMLParser.ID, alignment.getId()),
                        new Attribute.Default<String>("dataType", getAlignmentDataTypeDescription(alignment))
                }
        );

        for (int i = 0; i < alignment.getTaxonCount(); i++) {
            Taxon taxon = alignment.getTaxon(i);

            writer.writeOpenTag(SequenceParser.SEQUENCE);
            writer.writeIDref(TaxonParser.TAXON, taxon.getId());
            if (!options.samplePriorOnly) {
                writer.writeText(alignment.getAlignedSequenceString(i));
            } else {
                writer.writeText("N");
            }
            writer.writeCloseTag(SequenceParser.SEQUENCE);
        }
        writer.writeCloseTag(AlignmentParser.ALIGNMENT);
    }

    /**
     * Generate traits block regarding specific trait name (currently only <species>) from options
     *
     * @param writer     XMLWriter
     * @param trait      trait
     * @param traitType  traitType
     * @param taxonList  TaxonList
     */
    private void writeTraits(XMLWriter writer, TraitGuesser trait, TaxonList taxonList) {

        writer.writeText("");
        if (options.starBEASTOptions.isSpeciesAnalysis()) { // species
            writer.writeComment("Species definition: binds taxa, species and gene trees");
        }
        writer.writeComment("trait = " + trait.getTraitName() + " trait_type = " + trait.getTraitType());

        writer.writeOpenTag(trait.getTraitName(), new Attribute[]{new Attribute.Default<String>(XMLParser.ID, trait.getTraitName())});
        //new Attribute.Default<String>("traitType", traitType)});

        // write sub-tags for species
        if (options.starBEASTOptions.isSpeciesAnalysis()) { // species
            starEASTGeneratorGenerator.writeMultiSpecies(taxonList, writer);
        } // end write sub-tags for species

        writer.writeCloseTag(trait.getTraitName());

        if (options.starBEASTOptions.isSpeciesAnalysis()) { // species
            starEASTGeneratorGenerator.writeSTARBEAST(writer);
        }

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
                                new Attribute.Default<String>(XMLParser.ID, model.getPrefix(1) + partition.getName() + "." + SitePatternsParser.PATTERNS),
                        }
                );
//                for (PartitionData partition : options.dataPartitions) {
//                    if (partition.getPartitionSubstitutionModel() == model) {
                writePatternList(partition, 0, 3, writer);
                writePatternList(partition, 1, 3, writer);
//                    }
//                }
                writer.writeCloseTag(MergePatternsParser.MERGE_PATTERNS);

                writer.writeComment("The unique patterns for codon positions 3");
                writer.writeOpenTag(MergePatternsParser.MERGE_PATTERNS,
                        new Attribute[]{
                                new Attribute.Default<String>(XMLParser.ID, model.getPrefix(2) + partition.getName() + "." + SitePatternsParser.PATTERNS),
                        }
                );

//                for (PartitionData partition : options.dataPartitions) {
//                    if (partition.getPartitionSubstitutionModel() == model) {
                writePatternList(partition, 2, 3, writer);
//                    }
//                }

                writer.writeCloseTag(MergePatternsParser.MERGE_PATTERNS);

            } else {
                // pattern is 123
                // write pattern lists for all three codon positions
                for (int i = 1; i <= 3; i++) {
                    writer.writeComment("The unique patterns for codon positions " + i);
                    writer.writeOpenTag(MergePatternsParser.MERGE_PATTERNS,
                            new Attribute[]{
                                    new Attribute.Default<String>(XMLParser.ID, model.getPrefix(i) + partition.getName() + "." + SitePatternsParser.PATTERNS),
                            }
                    );

//                    for (PartitionData partition : options.dataPartitions) {
//                        if (partition.getPartitionSubstitutionModel() == model) {
                    writePatternList(partition, i - 1, 3, writer);
//                        }
//                    }

                    writer.writeCloseTag(MergePatternsParser.MERGE_PATTERNS);
                }

            }
        } else {
            //partitionCount = 1;
//            writer.writeComment("The unique patterns site patterns");
//            Alignment alignment = partition.getAlignment();

//            writer.writeOpenTag(SitePatternsParser.PATTERNS,
//                    new Attribute[]{
//                            new Attribute.Default<String>(XMLParser.ID, partition.getName() + "." + SitePatternsParser.PATTERNS),
//                    }
//            );
            writePatternList(partition, 0, 1, writer);
//            writer.writeIDref(AlignmentParser.ALIGNMENT, alignment.getId());
//            writer.writeCloseTag(SitePatternsParser.PATTERNS);

//            for (PartitionData partition : options.dataPartitions) {
//                if (partition.getPartitionSubstitutionModel() == model) {
//                    writePatternList(partition, 0, 1, writer);
//                }
//            }


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

        writer.writeComment("The unique patterns from " + from + " to " + (to > 0 ? to : "end") + ((every > 1) ? " every " + every : ""));

        // this object is created solely to calculate the number of patterns in the alignment
        SitePatterns patterns = new SitePatterns(alignment, from - 1, to - 1, every);
        writer.writeComment("npatterns=" + patterns.getPatternCount());

        List<Attribute> attributes = new ArrayList<Attribute>();

        // no codon, unique patterns site patterns
        if (offset == 0 && every == 1)
            attributes.add(new Attribute.Default<String>(XMLParser.ID, partition.getName() + "." + SitePatternsParser.PATTERNS));

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
     * Generate tmrca statistics
     *
     * @param writer the writer
     */
    public void writeTMRCAStatistics(XMLWriter writer) {

        writer.writeText("");
        for (Taxa taxa : options.taxonSets) {
            writer.writeOpenTag(
                    TMRCAStatistic.TMRCA_STATISTIC,
                    new Attribute[]{
                            new Attribute.Default<String>(XMLParser.ID, "tmrca(" + taxa.getId() + ")"),
                    }
            );
            writer.writeOpenTag(TMRCAStatistic.MRCA);
            writer.writeIDref(TaxaParser.TAXA, taxa.getId());
            writer.writeCloseTag(TMRCAStatistic.MRCA);
            writer.writeIDref(TreeModel.TREE_MODEL, TreeModel.TREE_MODEL);
            writer.writeCloseTag(TMRCAStatistic.TMRCA_STATISTIC);

            if (options.taxonSetsMono.get(taxa)) {
                writer.writeOpenTag(
                        MonophylyStatistic.MONOPHYLY_STATISTIC,
                        new Attribute[]{
                                new Attribute.Default<String>(XMLParser.ID, "monophyly(" + taxa.getId() + ")"),
                        });
                writer.writeOpenTag(MonophylyStatistic.MRCA);
                writer.writeIDref(TaxaParser.TAXA, taxa.getId());
                writer.writeCloseTag(MonophylyStatistic.MRCA);
                writer.writeIDref(TreeModel.TREE_MODEL, TreeModel.TREE_MODEL);
                writer.writeCloseTag(MonophylyStatistic.MONOPHYLY_STATISTIC);
            }
        }
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
     * @param writer  XMLWriter
     */
    public void writeMCMC(XMLWriter writer) {
        writer.writeComment("Define MCMC");
        writer.writeOpenTag(
                "mcmc",
                new Attribute[]{
                        new Attribute.Default<String>(XMLParser.ID, "mcmc"),
                        new Attribute.Default<Integer>("chainLength", options.chainLength),
                        new Attribute.Default<String>("autoOptimize", options.autoOptimize ? "true" : "false")
                });

        if (options.hasData()) {
            writer.writeOpenTag(CompoundLikelihood.POSTERIOR, new Attribute.Default<String>(XMLParser.ID, "posterior"));
        }

        // write prior block
        writer.writeOpenTag(CompoundLikelihood.PRIOR, new Attribute.Default<String>(XMLParser.ID, "prior"));

        if (options.starBEASTOptions.isSpeciesAnalysis()) { // species
            // coalescent prior
            writer.writeIDref(TreePartitionCoalescent.SPECIES_COALESCENT, TraitGuesser.Traits.TRAIT_SPECIES + "." + COALESCENT);
            // prior on population sizes
//            if (options.speciesTreePrior == TreePriorType.SPECIES_YULE) {
            writer.writeIDref(MixedDistributionLikelihood.DISTRIBUTION_LIKELIHOOD, SPOPS);
//            } else {
//                writer.writeIDref(SpeciesTreeBMPrior.STPRIOR, STP);
//            }
            // prior on species tree
            writer.writeIDref(SpeciationLikelihood.SPECIATION_LIKELIHOOD, SPECIATION_LIKE);
        }

        parameterPriorGenerator.writeParameterPriors(writer);

        for (PartitionTreeModel model : options.getPartitionTreeModels()) {
            PartitionTreePrior prior = model.getPartitionTreePrior();
            treePriorGenerator.writePriorLikelihoodReference(prior, model, writer);
            writer.writeText("");
        }

        for (PartitionTreePrior prior : options.getPartitionTreePriors()) {
            treePriorGenerator.writeEBSPVariableDemographicReference(prior, writer);
        }

        generateInsertionPoint(ComponentGenerator.InsertionPoint.IN_MCMC_PRIOR, writer);

        writer.writeCloseTag(CompoundLikelihood.PRIOR);

        if (options.hasData()) {
            // write likelihood block
            writer.writeOpenTag(CompoundLikelihood.LIKELIHOOD, new Attribute.Default<String>(XMLParser.ID, "likelihood"));

            treeLikelihoodGenerator.writeTreeLikelihoodReferences(writer);

            generateInsertionPoint(ComponentGenerator.InsertionPoint.IN_MCMC_LIKELIHOOD, writer);

            writer.writeCloseTag(CompoundLikelihood.LIKELIHOOD);

            writer.writeCloseTag(CompoundLikelihood.POSTERIOR);
        }

        writer.writeIDref(SimpleOperatorSchedule.OPERATOR_SCHEDULE, "operators");

        // write log to screen    	
        logGenerator.writeLogToScreen(writer, branchRatesModelGenerator);
        // write log to file
        logGenerator.writeLogToFile(writer, treePriorGenerator, branchRatesModelGenerator,
                        substitutionModelGenerator, treeLikelihoodGenerator);
        // write tree log to file
        logGenerator.writeTreeLogToFile(writer);

        writer.writeCloseTag("mcmc");
    }

}
