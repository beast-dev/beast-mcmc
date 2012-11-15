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

package dr.app.oldbeauti;

import dr.app.beast.BeastVersion;
import dr.app.beauti.generator.InitialTreeGenerator;
import dr.evolution.alignment.SitePatterns;
import dr.evolution.datatype.DataType;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.datatype.TwoStateCovarion;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxa;
import dr.evolution.util.Taxon;
import dr.evolution.util.TaxonList;
import dr.evolution.util.Units;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.coalescent.VariableDemographicModel;
import dr.evomodel.sitemodel.GammaSiteModel;
import dr.evomodel.sitemodel.SiteModel;
import dr.evomodel.speciation.BirthDeathGernhard08Model;
import dr.evomodel.substmodel.AminoAcidModelType;
import dr.evomodel.substmodel.NucModelType;
import dr.evomodel.tree.TreeModel;
import dr.evomodelxml.CSVExporterParser;
import dr.evomodelxml.branchratemodel.DiscretizedBranchRatesParser;
import dr.evomodelxml.branchratemodel.RandomLocalClockModelParser;
import dr.evomodelxml.branchratemodel.StrictClockBranchRatesParser;
import dr.evomodelxml.coalescent.*;
import dr.evomodelxml.coalescent.operators.SampleNonActiveGibbsOperatorParser;
import dr.evomodelxml.operators.ExchangeOperatorParser;
import dr.evomodelxml.operators.SubtreeSlideOperatorParser;
import dr.evomodelxml.operators.TreeBitMoveOperatorParser;
import dr.evomodelxml.operators.WilsonBaldingParser;
import dr.evomodelxml.sitemodel.GammaSiteModelParser;
import dr.evomodelxml.speciation.BirthDeathModelParser;
import dr.evomodelxml.speciation.SpeciationLikelihoodParser;
import dr.evomodelxml.speciation.YuleModelParser;
import dr.evomodelxml.substmodel.*;
import dr.evomodelxml.tree.*;
import dr.evomodelxml.treelikelihood.TreeLikelihoodParser;
import dr.evoxml.*;
import dr.inference.distribution.ExponentialDistributionModel;
import dr.inference.distribution.ExponentialMarkovModel;
import dr.inference.model.ParameterParser;
import dr.inferencexml.distribution.*;
import dr.inferencexml.loggers.ColumnsParser;
import dr.inferencexml.loggers.LoggerParser;
import dr.inferencexml.model.*;
import dr.inferencexml.operators.*;
import dr.util.Attribute;
import dr.util.Version;
import dr.xml.XMLParser;

import java.io.Writer;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

/**
 * This class holds all the data for the current BEAUti Document
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: BeastGenerator.java,v 1.4 2006/09/05 13:29:34 rambaut Exp $
 */
public class BeastGenerator extends BeautiOptions {

    private final static Version version = new BeastVersion();

    public BeastGenerator() {
        super();
    }

    /**
     * Checks various options to check they are valid. Throws IllegalArgumentExceptions with
     * descriptions of the problems.
     *
     * @throws IllegalArgumentException if there is a problem with the current settings
     */
    public void checkOptions() throws IllegalArgumentException {
        Set<String> ids = new HashSet<String>();

        ids.add(TaxaParser.TAXA);
        ids.add(AlignmentParser.ALIGNMENT);

        if (taxonList != null) {
            for (int i = 0; i < taxonList.getTaxonCount(); i++) {
                Taxon taxon = taxonList.getTaxon(i);
                if (ids.contains(taxon.getId())) {
                    throw new IllegalArgumentException("A taxon has the same id," + taxon.getId() +
                            "\nas another element (taxon, sequence, taxon set etc.):\nAll ids should be unique.");
                }
                ids.add(taxon.getId());
            }
        }

        for (Taxa taxa : taxonSets) {
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

        getPartionCount(codonHeteroPattern);
    }

    /**
     * Generate a beast xml file from these beast options
     *
     * @param w the writer
     */
    public void generateXML(Writer w) {

        XMLWriter writer = new XMLWriter(w);

        writer.writeText("<?xml version=\"1.0\" standalone=\"yes\"?>");
        writer.writeComment("Generated by BEAUti " + version.getVersionString());
        writer.writeComment("      by Alexei J. Drummond, Andrew Rambaut and Marc A. Suchard");
        writer.writeComment("      Department of Computer Science, University of Auckland,");
        writer.writeComment("      Institute of Evolutionary Biology, University of Edinburgh and");
        writer.writeComment("      David Geffen School of Medicine, University of California, Los Angeles");
        writer.writeComment("      http://beast.bio.ed.ac.uk/");
        writer.writeOpenTag("beast");
        writer.writeText("");
        writeTaxa(writer);

        if (taxonSets != null && taxonSets.size() > 0) {
            writeTaxonSets(writer);
        }

        if (alignment != null) {
            writeAlignment(writer);
            writePatternLists(writer);
        }

        writer.writeText("");
        writeNodeHeightPriorModel(writer);

        writer.writeText("");
        writeStartingTree(writer);
        writer.writeText("");
        writeTreeModel(writer);
        writer.writeText("");
        writeNodeHeightPrior(writer);
        if (nodeHeightPrior == LOGISTIC) {
            writer.writeText("");
            writeBooleanLikelihood(writer);
        } else if (nodeHeightPrior == SKYLINE) {
            writer.writeText("");
            writeExponentialMarkovLikelihood(writer);
        }

        writer.writeText("");
        writeBranchRatesModel(writer);

        if (alignment != null) {
            writer.writeText("");
            writeSubstitutionModel(writer);
            writer.writeText("");
            writeSiteModel(writer);
            writer.writeText("");
            writeTreeLikelihood(writer);
        }

        writer.writeText("");

        if (taxonSets != null && taxonSets.size() > 0) {
            writeTMRCAStatistics(writer);
        }

        ArrayList<Operator> operators = selectOperators();
        writeOperatorSchedule(operators, writer);
        writer.writeText("");
        writeMCMC(writer);
        writer.writeText("");
        writeTimerReport(writer);
        writer.writeText("");
        if (performTraceAnalysis) {
            writeTraceAnalysis(writer);
        }
        if (generateCSV) {
            writeAnalysisToCSVfile(writer);
        }

        writer.writeCloseTag("beast");
        writer.flush();
    }

    /**
     * Generate a taxa block from these beast options
     *
     * @param writer the writer
     */
    public void writeTaxa(XMLWriter writer) {

        writer.writeComment("The list of taxa to be analysed (can also include dates/ages).");
        writer.writeComment("ntax=" + taxonList.getTaxonCount());
        writer.writeOpenTag(TaxaParser.TAXA, new Attribute[]{new Attribute.Default<String>(XMLParser.ID, TaxaParser.TAXA)});

        boolean firstDate = true;
        for (int i = 0; i < taxonList.getTaxonCount(); i++) {
            Taxon taxon = taxonList.getTaxon(i);

            boolean hasDate = false;

            if (maximumTipHeight > 0.0) {
                hasDate = TaxonList.Utils.hasAttribute(taxonList, i, dr.evolution.util.Date.DATE);
            }

            writer.writeTag(TaxonParser.TAXON, new Attribute[]{new Attribute.Default<String>(XMLParser.ID, taxon.getId())}, !hasDate);

            if (hasDate) {
                dr.evolution.util.Date date = (dr.evolution.util.Date) taxon.getAttribute(dr.evolution.util.Date.DATE);

                if (firstDate) {
                    units = date.getUnits();
                    firstDate = false;
                } else {
                    if (units != date.getUnits()) {
                        System.err.println("Error: Units in dates do not match.");
                    }
                }

                Attribute[] attributes = {
                        new Attribute.Default<Double>(DateParser.VALUE, date.getTimeValue()),
                        new Attribute.Default<String>(DateParser.DIRECTION, date.isBackwards() ? DateParser.BACKWARDS : DateParser.FORWARDS),
                        new Attribute.Default<String>(DateParser.UNITS, Units.Utils.getDefaultUnitName(units))
                        /*,
                                                                                new Attribute.Default("origin", date.getOrigin()+"")*/
                };

                writer.writeTag(dr.evolution.util.Date.DATE, attributes, true);
                writer.writeCloseTag(TaxonParser.TAXON);
            }
        }

        writer.writeCloseTag(TaxaParser.TAXA);
    }

    /**
     * Generate additional taxon sets
     *
     * @param writer the writer
     */
    public void writeTaxonSets(XMLWriter writer) {

        writer.writeText("");
        for (Taxa taxa : taxonSets) {
            writer.writeOpenTag(
                    TaxaParser.TAXA,
                    new Attribute[]{
                            new Attribute.Default<String>(XMLParser.ID, taxa.getId())
                    }
            );

            for (int j = 0; j < taxa.getTaxonCount(); j++) {
                Taxon taxon = taxa.getTaxon(j);

                writer.writeTag(TaxonParser.TAXON, new Attribute[]{new Attribute.Default<String>(XMLParser.IDREF, taxon.getId())}, true);
            }
            writer.writeCloseTag(TaxaParser.TAXA);
        }
    }

    /**
     * Determine and return the datatype description for these beast options
     * note that the datatype in XML may differ from the actual datatype
     *
     * @return description
     */

    private String getAlignmentDataTypeDescription() {
        String description;

        switch (dataType) {
            case DataType.TWO_STATES:
            case DataType.COVARION:

                switch (binarySubstitutionModel) {
                    case BIN_COVARION:
                        description = TwoStateCovarion.DESCRIPTION;
                        break;

                    default:
                        description = alignment.getDataType().getDescription();
                }
                break;

            default:
                description = alignment.getDataType().getDescription();
        }

        return description;
    }


    /**
     * Generate an alignment block from these beast options
     *
     * @param writer the writer
     */
    public void writeAlignment(XMLWriter writer) {

        writer.writeText("");
        writer.writeComment("The sequence alignment (each sequence refers to a taxon above).");
        writer.writeComment("ntax=" + alignment.getTaxonCount() + " nchar=" + alignment.getSiteCount());
        if (samplePriorOnly) {
            writer.writeComment("Null sequences generated in order to sample from the prior only.");
        }


        writer.writeOpenTag(
                AlignmentParser.ALIGNMENT,
                new Attribute[]{
                        new Attribute.Default<String>(XMLParser.ID, AlignmentParser.ALIGNMENT),
                        new Attribute.Default<String>("dataType", getAlignmentDataTypeDescription())
                }
        );

        for (int i = 0; i < alignment.getTaxonCount(); i++) {
            Taxon taxon = alignment.getTaxon(i);

            writer.writeOpenTag("sequence");
            writer.writeTag(TaxonParser.TAXON, new Attribute[]{new Attribute.Default<String>(XMLParser.IDREF, taxon.getId())}, true);
            if (!samplePriorOnly) {
                writer.writeText(alignment.getAlignedSequenceString(i));
            } else {
                // 3 Ns written in case 3 codon positions selected...
                writer.writeText("NNN");
            }
            writer.writeCloseTag("sequence");
        }
        writer.writeCloseTag(AlignmentParser.ALIGNMENT);
    }

    /**
     * Write a demographic model
     *
     * @param writer the writer
     */
    public void writeNodeHeightPriorModel(XMLWriter writer) {

        String initialPopSize = null;

        if (nodeHeightPrior == CONSTANT) {

            writer.writeComment("A prior assumption that the population size has remained constant");
            writer.writeComment("throughout the time spanned by the genealogy.");
            writer.writeOpenTag(
                    ConstantPopulationModelParser.CONSTANT_POPULATION_MODEL,
                    new Attribute[]{
                            new Attribute.Default<String>(XMLParser.ID, "constant"),
                            new Attribute.Default<String>("units", Units.Utils.getDefaultUnitName(units))
                    }
            );

            writer.writeOpenTag(ConstantPopulationModelParser.POPULATION_SIZE);
            writeParameter("constant.popSize", writer);
            writer.writeCloseTag(ConstantPopulationModelParser.POPULATION_SIZE);
            writer.writeCloseTag(ConstantPopulationModelParser.CONSTANT_POPULATION_MODEL);

        } else if (nodeHeightPrior == EXPONENTIAL) {
            // generate an exponential prior tree

            writer.writeComment("A prior assumption that the population size has grown exponentially");
            writer.writeComment("throughout the time spanned by the genealogy.");
            writer.writeOpenTag(
                    ExponentialGrowthModelParser.EXPONENTIAL_GROWTH_MODEL,
                    new Attribute[]{
                            new Attribute.Default<String>(XMLParser.ID, "exponential"),
                            new Attribute.Default<String>("units", Units.Utils.getDefaultUnitName(units))
                    }
            );

            // write pop size socket
            writer.writeOpenTag(ExponentialGrowthModelParser.POPULATION_SIZE);
            writeParameter("exponential.popSize", writer);
            writer.writeCloseTag(ExponentialGrowthModelParser.POPULATION_SIZE);

            if (parameterization == GROWTH_RATE) {
                // write growth rate socket
                writer.writeOpenTag(ExponentialGrowthModelParser.GROWTH_RATE);
                writeParameter("exponential.growthRate", writer);
                writer.writeCloseTag(ExponentialGrowthModelParser.GROWTH_RATE);
            } else {
                // write doubling time socket
                writer.writeOpenTag(ExponentialGrowthModelParser.DOUBLING_TIME);
                writeParameter("exponential.doublingTime", writer);
                writer.writeCloseTag(ExponentialGrowthModelParser.DOUBLING_TIME);
            }

            writer.writeCloseTag(ExponentialGrowthModelParser.EXPONENTIAL_GROWTH_MODEL);
        } else if (nodeHeightPrior == LOGISTIC) {
            // generate an exponential prior tree

            writer.writeComment("A prior assumption that the population size has grown logistically");
            writer.writeComment("throughout the time spanned by the genealogy.");
            writer.writeOpenTag(
                    LogisticGrowthModelParser.LOGISTIC_GROWTH_MODEL,
                    new Attribute[]{
                            new Attribute.Default<String>(XMLParser.ID, "logistic"),
                            new Attribute.Default<String>("units", Units.Utils.getDefaultUnitName(units))
                    }
            );

            // write pop size socket
            writer.writeOpenTag(LogisticGrowthModelParser.POPULATION_SIZE);
            writeParameter("logistic.popSize", writer);
            writer.writeCloseTag(LogisticGrowthModelParser.POPULATION_SIZE);

            if (parameterization == GROWTH_RATE) {
                // write growth rate socket
                writer.writeOpenTag(LogisticGrowthModelParser.GROWTH_RATE);
                writeParameter("logistic.growthRate", writer);
                writer.writeCloseTag(LogisticGrowthModelParser.GROWTH_RATE);
            } else {
                // write doubling time socket
                writer.writeOpenTag(LogisticGrowthModelParser.DOUBLING_TIME);
                writeParameter("logistic.doublingTime", writer);
                writer.writeCloseTag(LogisticGrowthModelParser.DOUBLING_TIME);
            }

            // write logistic t50 socket
            writer.writeOpenTag(LogisticGrowthModelParser.TIME_50);
            writeParameter("logistic.t50", writer);
            writer.writeCloseTag(LogisticGrowthModelParser.TIME_50);

            writer.writeCloseTag(LogisticGrowthModelParser.LOGISTIC_GROWTH_MODEL);

            initialPopSize = "logistic.popSize";

        } else if (nodeHeightPrior == EXPANSION) {
            // generate an exponential prior tree

            writer.writeComment("A prior assumption that the population size has grown exponentially");
            writer.writeComment("from some ancestral population size in the past.");
            writer.writeOpenTag(
                    ExpansionModelParser.EXPANSION_MODEL,
                    new Attribute[]{
                            new Attribute.Default<String>(XMLParser.ID, "expansion"),
                            new Attribute.Default<String>("units", Units.Utils.getDefaultUnitName(units))
                    }
            );

            // write pop size socket
            writer.writeOpenTag(ExpansionModelParser.POPULATION_SIZE);
            writeParameter("expansion.popSize", writer);
            writer.writeCloseTag(ExpansionModelParser.POPULATION_SIZE);

            if (parameterization == GROWTH_RATE) {
                // write growth rate socket
                writer.writeOpenTag(ExpansionModelParser.GROWTH_RATE);
                writeParameter("expansion.growthRate", writer);
                writer.writeCloseTag(ExpansionModelParser.GROWTH_RATE);
            } else {
                // write doubling time socket
                writer.writeOpenTag(ExpansionModelParser.DOUBLING_TIME);
                writeParameter("expansion.doublingTime", writer);
                writer.writeCloseTag(ExpansionModelParser.DOUBLING_TIME);
            }

            // write ancestral proportion socket
            writer.writeOpenTag(ExpansionModelParser.ANCESTRAL_POPULATION_PROPORTION);
            writeParameter("expansion.ancestralProportion", writer);
            writer.writeCloseTag(ExpansionModelParser.ANCESTRAL_POPULATION_PROPORTION);

            writer.writeCloseTag(ExpansionModelParser.EXPANSION_MODEL);

            initialPopSize = "expansion.popSize";

        } else if (nodeHeightPrior == YULE) {
            writer.writeComment("A prior on the distribution node heights defined given");
            writer.writeComment("a Yule speciation process (a pure birth process).");
            writer.writeOpenTag(
                    YuleModelParser.YULE_MODEL,
                    new Attribute[]{
                            new Attribute.Default<String>(XMLParser.ID, "yule"),
                            new Attribute.Default<String>("units", Units.Utils.getDefaultUnitName(units))
                    }
            );

            writer.writeOpenTag(YuleModelParser.BIRTH_RATE);
            writeParameter("yule.birthRate", writer);
            writer.writeCloseTag(YuleModelParser.BIRTH_RATE);
            writer.writeCloseTag(YuleModelParser.YULE_MODEL);
        } else if (nodeHeightPrior == BIRTH_DEATH) {
            writer.writeComment("A prior on the distribution node heights defined given");
            writer.writeComment("a Birth-Death speciation process (Gernhard 2008).");
            writer.writeOpenTag(
                    BirthDeathGernhard08Model.BIRTH_DEATH_MODEL,
                    new Attribute[]{
                            new Attribute.Default<String>(XMLParser.ID, "birthDeath"),
                            new Attribute.Default<String>("units", Units.Utils.getDefaultUnitName(units))
                    }
            );

            writer.writeOpenTag(BirthDeathModelParser.BIRTHDIFF_RATE);
            writeParameter(BirthDeathModelParser.MEAN_GROWTH_RATE_PARAM_NAME, writer);
            writer.writeCloseTag(BirthDeathModelParser.BIRTHDIFF_RATE);
            writer.writeOpenTag(BirthDeathModelParser.RELATIVE_DEATH_RATE);
            writeParameter(BirthDeathModelParser.RELATIVE_DEATH_RATE_PARAM_NAME, writer);
            writer.writeCloseTag(BirthDeathModelParser.RELATIVE_DEATH_RATE);

            writer.writeCloseTag(BirthDeathGernhard08Model.BIRTH_DEATH_MODEL);
        }

        if (nodeHeightPrior != CONSTANT && nodeHeightPrior != EXPONENTIAL) {
            // If the node height prior is not one of these two then we need to simulate a
            // random starting tree under a constant size coalescent.

            writer.writeComment("This is a simple constant population size coalescent model");
            writer.writeComment("that is used to generate an initial tree for the chain.");
            writer.writeOpenTag(
                    ConstantPopulationModelParser.CONSTANT_POPULATION_MODEL,
                    new Attribute[]{
                            new Attribute.Default<String>(XMLParser.ID, "initialDemo"),
                            new Attribute.Default<String>("units", Units.Utils.getDefaultUnitName(units))
                    }
            );

            writer.writeOpenTag(ConstantPopulationModelParser.POPULATION_SIZE);
            if (initialPopSize != null) {
                writer.writeTag(ParameterParser.PARAMETER,
                        new Attribute[]{
                                new Attribute.Default<String>(XMLParser.IDREF, initialPopSize),
                        }, true);
            } else {
                writeParameter("initialDemo.popSize", 1, 100.0, Double.NaN, Double.NaN, writer);
            }
            writer.writeCloseTag(ConstantPopulationModelParser.POPULATION_SIZE);
            writer.writeCloseTag(ConstantPopulationModelParser.CONSTANT_POPULATION_MODEL);
        }

    }

    /**
     * Writes the pattern lists
     *
     * @param writer the writer
     */
    public void writePatternLists(XMLWriter writer) {

        partitionCount = getPartionCount(codonHeteroPattern);

        writer.writeText("");
        if (alignment.getDataType() == Nucleotides.INSTANCE && codonHeteroPattern != null && partitionCount > 1) {

            if (codonHeteroPattern.equals("112")) {
                writer.writeComment("The unique patterns for codon positions 1 & 2");
                writer.writeOpenTag(MergePatternsParser.MERGE_PATTERNS,
                        new Attribute[]{
                                new Attribute.Default<String>(XMLParser.ID, "patterns1+2"),
                        }
                );
                writePatternList(1, 3, writer);
                writePatternList(2, 3, writer);
                writer.writeCloseTag(MergePatternsParser.MERGE_PATTERNS);

                writePatternList(3, 3, writer);

            } else {
                // pattern is 123
                // write pattern lists for all three codon positions
                for (int i = 1; i <= 3; i++) {
                    writePatternList(i, 3, writer);
                }

            }
        } else {
            partitionCount = 1;
            writePatternList(-1, 0, writer);
        }
    }

    private int getPartionCount(String codonPattern) {

        if (codonPattern == null || codonPattern.equals("111")) {
            return 1;
        }
        if (codonPattern.equals("123")) {
            return 3;
        }
        if (codonPattern.equals("112")) {
            return 2;
        }
        throw new IllegalArgumentException("codonPattern must be one of '111', '112' or '123'");
    }

    /**
     * Write a single pattern list
     *
     * @param writer the writer
     * @param from   from site
     * @param every  skip every
     */
    private void writePatternList(int from, int every, XMLWriter writer) {

        String id = SitePatternsParser.PATTERNS;
        if (from < 1) {
            writer.writeComment("The unique patterns for all positions");
            from = 1;
        } else {
            writer.writeComment("The unique patterns for codon position " + from);
            id += Integer.toString(from);
        }

        SitePatterns patterns = new SitePatterns(alignment, from - 1, 0, every);
        writer.writeComment("npatterns=" + patterns.getPatternCount());
        if (every != 0) {
            writer.writeOpenTag(SitePatternsParser.PATTERNS,
                    new Attribute[]{
                            new Attribute.Default<String>(XMLParser.ID, id),
                            new Attribute.Default<String>("from", "" + from),
                            new Attribute.Default<String>("every", "" + every)
                    }
            );
        } else {
            writer.writeOpenTag(SitePatternsParser.PATTERNS,
                    new Attribute[]{
                            new Attribute.Default<String>(XMLParser.ID, id),
                            new Attribute.Default<String>("from", "" + from)
                    }
            );
        }

        writer.writeTag(AlignmentParser.ALIGNMENT, new Attribute.Default<String>(XMLParser.IDREF, AlignmentParser.ALIGNMENT), true);
        writer.writeCloseTag(SitePatternsParser.PATTERNS);
    }

    /**
     * Write tree model XML block.
     *
     * @param writer the writer
     */
    private void writeTreeModel(XMLWriter writer) {

        writer.writeTag(TreeModel.TREE_MODEL, new Attribute.Default<String>(XMLParser.ID, "treeModel"), false);

        if (userTree) {
            writer.writeTag("tree", new Attribute.Default<String>(XMLParser.IDREF, InitialTreeGenerator.STARTING_TREE), true);
        } else {
            writer.writeTag(CoalescentSimulatorParser.COALESCENT_TREE, new Attribute.Default<String>(XMLParser.IDREF, InitialTreeGenerator.STARTING_TREE), true);
        }

        writer.writeOpenTag(TreeModelParser.ROOT_HEIGHT);
        writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default<String>(XMLParser.ID, "treeModel.rootHeight"), true);
        writer.writeCloseTag(TreeModelParser.ROOT_HEIGHT);


        writer.writeOpenTag(TreeModelParser.NODE_HEIGHTS, new Attribute.Default<String>(TreeModelParser.INTERNAL_NODES, "true"));
        writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default<String>(XMLParser.ID, "treeModel.internalNodeHeights"), true);
        writer.writeCloseTag(TreeModelParser.NODE_HEIGHTS);

        writer.writeOpenTag(TreeModelParser.NODE_HEIGHTS,
                new Attribute[]{
                        new Attribute.Default<String>(TreeModelParser.INTERNAL_NODES, "true"),
                        new Attribute.Default<String>(TreeModelParser.ROOT_NODE, "true")
                });
        writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default<String>(XMLParser.ID, "treeModel.allInternalNodeHeights"), true);
        writer.writeCloseTag(TreeModelParser.NODE_HEIGHTS);

        if (clockModel == RANDOM_LOCAL_CLOCK) {
            writer.writeOpenTag(TreeModelParser.NODE_RATES,
                    new Attribute[]{
                            new Attribute.Default<String>(TreeModelParser.ROOT_NODE, "false"),
                            new Attribute.Default<String>(TreeModelParser.INTERNAL_NODES, "true"),
                            new Attribute.Default<String>(TreeModelParser.LEAF_NODES, "true")
                    });
            writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default<String>(XMLParser.ID, LOCAL_CLOCK + "." + "rates"), true);
            writer.writeCloseTag(TreeModelParser.NODE_RATES);

            writer.writeOpenTag(TreeModelParser.NODE_TRAITS,
                    new Attribute[]{
                            new Attribute.Default<String>(TreeModelParser.ROOT_NODE, "false"),
                            new Attribute.Default<String>(TreeModelParser.INTERNAL_NODES, "true"),
                            new Attribute.Default<String>(TreeModelParser.LEAF_NODES, "true")
                    });
            writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default<String>(XMLParser.ID, LOCAL_CLOCK + "." + "changes"), true);
            writer.writeCloseTag(TreeModelParser.NODE_TRAITS);
        }

        writer.writeCloseTag(TreeModel.TREE_MODEL);
    }

    /**
     * Writes the substitution model to XML.
     *
     * @param writer the writer
     */
    public void writeSubstitutionModel(XMLWriter writer) {


        switch (dataType) {
            case DataType.NUCLEOTIDES:
                // Jukes-Cantor model
                if (nucSubstitutionModel == JC) {
                    writer.writeComment("The JC substitution model (Jukes & Cantor, 1969)");
                    writer.writeOpenTag(
                            NucModelType.HKY.getXMLName(),
                            new Attribute[]{new Attribute.Default<String>(XMLParser.ID, "jc")}
                    );
                    writer.writeOpenTag(HKYParser.FREQUENCIES);
                    writer.writeOpenTag(
                            FrequencyModelParser.FREQUENCY_MODEL,
                            new Attribute[]{
                                    new Attribute.Default<String>("dataType", alignment.getDataType().getDescription())
                            }
                    );
                    writer.writeOpenTag(FrequencyModelParser.FREQUENCIES);
                    writer.writeTag(
                            ParameterParser.PARAMETER,
                            new Attribute[]{
                                    new Attribute.Default<String>(XMLParser.ID, "jc.frequencies"),
                                    new Attribute.Default<String>("value", "0.25 0.25 0.25 0.25")
                            },
                            true
                    );
                    writer.writeCloseTag(FrequencyModelParser.FREQUENCIES);

                    writer.writeCloseTag(FrequencyModelParser.FREQUENCY_MODEL);
                    writer.writeCloseTag(HKYParser.FREQUENCIES);

                    writer.writeOpenTag(HKYParser.KAPPA);
                    writeParameter("jc.kappa", 1, 1.0, Double.NaN, Double.NaN, writer);
                    writer.writeCloseTag(HKYParser.KAPPA);
                    writer.writeCloseTag(NucModelType.HKY.getXMLName());

                } else {
                    // Hasegawa Kishino and Yano 85 model
                    if (nucSubstitutionModel == HKY) {
                        if (unlinkedSubstitutionModel) {
                            for (int i = 1; i <= partitionCount; i++) {
                                writeHKYModel(i, writer);
                            }
                        } else {
                            writeHKYModel(-1, writer);
                        }
                    } else {
                        // General time reversible model
                        if (nucSubstitutionModel == GTR) {
                            if (unlinkedSubstitutionModel) {
                                for (int i = 1; i <= partitionCount; i++) {
                                    writeGTRModel(i, writer);
                                }
                            } else {
                                writeGTRModel(-1, writer);
                            }
                        }
                    }
                }
                break;

            case DataType.AMINO_ACIDS:
                // Amino Acid model
                String aaModel = "";

                switch (aaSubstitutionModel) {
                    case 0:
                        aaModel = AminoAcidModelType.BLOSUM_62.getXMLName();
                        break;
                    case 1:
                        aaModel = AminoAcidModelType.DAYHOFF.getXMLName();
                        break;
                    case 2:
                        aaModel = AminoAcidModelType.JTT.getXMLName();
                        break;
                    case 3:
                        aaModel = AminoAcidModelType.MT_REV_24.getXMLName();
                        break;
                    case 4:
                        aaModel = AminoAcidModelType.CP_REV_45.getXMLName();
                        break;
                    case 5:
                        aaModel = AminoAcidModelType.WAG.getXMLName();
                        break;
                }

                writer.writeComment("The " + aaModel + " substitution model");
                writer.writeTag(
                        EmpiricalAminoAcidModelParser.EMPIRICAL_AMINO_ACID_MODEL,
                        new Attribute[]{new Attribute.Default<String>(XMLParser.ID, "aa"),
                                new Attribute.Default<String>("type", aaModel)}, true
                );

                break;

            case DataType.TWO_STATES:
            case DataType.COVARION:

                switch (binarySubstitutionModel) {
                    case BIN_SIMPLE:
                        writeBinarySimpleModel(writer);
                        break;
                    case BIN_COVARION:
                        writeBinaryCovarionModel(writer);
                        break;
                }

                break;
        }
    }

    /**
     * Write the HKY model XML block.
     *
     * @param num    the model number
     * @param writer the writer
     */
    public void writeHKYModel(int num, XMLWriter writer) {
        String id = "hky";
        if (num > 0) {
            id += Integer.toString(num);
        }
        // Hasegawa Kishino and Yano 85 model
        writer.writeComment("The HKY substitution model (Hasegawa, Kishino & Yano, 1985)");
        writer.writeOpenTag(
                NucModelType.HKY.getXMLName(),
                new Attribute[]{new Attribute.Default<String>(XMLParser.ID, id)}
        );
        writer.writeOpenTag(HKYParser.FREQUENCIES);
        writer.writeOpenTag(
                FrequencyModelParser.FREQUENCY_MODEL,
                new Attribute[]{
                        new Attribute.Default<String>("dataType", alignment.getDataType().getDescription())
                }
        );
        writer.writeTag(AlignmentParser.ALIGNMENT, new Attribute[]{new Attribute.Default<String>(XMLParser.IDREF, AlignmentParser.ALIGNMENT)}, true);
        writer.writeOpenTag(FrequencyModelParser.FREQUENCIES);
        if (frequencyPolicy == ALLEQUAL)
            writeParameter(id + ".frequencies", 4, writer);
        else
            writeParameter(id + ".frequencies", 4, Double.NaN, Double.NaN, Double.NaN, writer);
        writer.writeCloseTag(FrequencyModelParser.FREQUENCIES);
        writer.writeCloseTag(FrequencyModelParser.FREQUENCY_MODEL);
        writer.writeCloseTag(HKYParser.FREQUENCIES);

        writer.writeOpenTag(HKYParser.KAPPA);
        writeParameter(id + ".kappa", writer);
        writer.writeCloseTag(HKYParser.KAPPA);
        writer.writeCloseTag(NucModelType.HKY.getXMLName());
    }

    /**
     * Write the GTR model XML block.
     *
     * @param num    the model number
     * @param writer the writer
     */
    public void writeGTRModel(int num, XMLWriter writer) {
        String id = "gtr";
        if (num > 0) {
            id += Integer.toString(num);
        }

        writer.writeComment("The general time reversible (GTR) substitution model");
        writer.writeOpenTag(
                GTRParser.GTR_MODEL,
                new Attribute[]{new Attribute.Default<String>(XMLParser.ID, id)}
        );
        writer.writeOpenTag(GTRParser.FREQUENCIES);
        writer.writeOpenTag(
                FrequencyModelParser.FREQUENCY_MODEL,
                new Attribute[]{
                        new Attribute.Default<String>("dataType", alignment.getDataType().getDescription())
                }
        );
        writer.writeTag(AlignmentParser.ALIGNMENT, new Attribute[]{new Attribute.Default<String>(XMLParser.IDREF, AlignmentParser.ALIGNMENT)}, true);
        writer.writeOpenTag(FrequencyModelParser.FREQUENCIES);
        if (frequencyPolicy == ALLEQUAL)
            writeParameter(id + ".frequencies", 4, writer);
        else
            writeParameter(id + ".frequencies", 4, Double.NaN, Double.NaN, Double.NaN, writer);
        writer.writeCloseTag(FrequencyModelParser.FREQUENCIES);
        writer.writeCloseTag(FrequencyModelParser.FREQUENCY_MODEL);
        writer.writeCloseTag(GTRParser.FREQUENCIES);

        writer.writeOpenTag(GTRParser.A_TO_C);
        writeParameter(id + ".ac", writer);
        writer.writeCloseTag(GTRParser.A_TO_C);

        writer.writeOpenTag(GTRParser.A_TO_G);
        writeParameter(id + ".ag", writer);
        writer.writeCloseTag(GTRParser.A_TO_G);

        writer.writeOpenTag(GTRParser.A_TO_T);
        writeParameter(id + ".at", writer);
        writer.writeCloseTag(GTRParser.A_TO_T);

        writer.writeOpenTag(GTRParser.C_TO_G);
        writeParameter(id + ".cg", writer);
        writer.writeCloseTag(GTRParser.C_TO_G);

        writer.writeOpenTag(GTRParser.G_TO_T);
        writeParameter(id + ".gt", writer);
        writer.writeCloseTag(GTRParser.G_TO_T);
        writer.writeCloseTag(GTRParser.GTR_MODEL);
    }


    /**
     * Write the Binary  simple model XML block.
     *
     * @param writer the writer
     */
    public void writeBinarySimpleModel(XMLWriter writer) {
        final String id = "bsimple";

        writer.writeComment("The Binary simple model (based on the general substitution model)");
        writer.writeOpenTag(
                BinarySubstitutionModelParser.BINARY_SUBSTITUTION_MODEL,
                new Attribute[]{new Attribute.Default<String>(XMLParser.ID, id)}
        );
        writer.writeOpenTag(GeneralSubstitutionModelParser.FREQUENCIES);
        writer.writeOpenTag(
                FrequencyModelParser.FREQUENCY_MODEL,
                new Attribute[]{
                        new Attribute.Default<String>("dataType", alignment.getDataType().getDescription())
                }
        );
        writer.writeTag(AlignmentParser.ALIGNMENT, new Attribute[]{new Attribute.Default<String>(XMLParser.IDREF, AlignmentParser.ALIGNMENT)}, true);
        writer.writeOpenTag(FrequencyModelParser.FREQUENCIES);
        writeParameter(id + ".frequencies", 2, Double.NaN, Double.NaN, Double.NaN, writer);
        writer.writeCloseTag(FrequencyModelParser.FREQUENCIES);
        writer.writeCloseTag(FrequencyModelParser.FREQUENCY_MODEL);
        writer.writeCloseTag(GeneralSubstitutionModelParser.FREQUENCIES);

        writer.writeCloseTag(BinarySubstitutionModelParser.BINARY_SUBSTITUTION_MODEL);
    }


    /**
     * Write the Binary covarion model XML block
     *
     * @param writer the writer
     */

    public void writeBinaryCovarionModel(XMLWriter writer) {
        String id = "bcov";

        writer.writeComment("The Binary covarion model");
        writer.writeOpenTag(
                BinaryCovarionModelParser.COVARION_MODEL,
                new Attribute[]{new Attribute.Default<String>(XMLParser.ID, id)}
        );

        writer.writeOpenTag(BinaryCovarionModelParser.FREQUENCIES);
        writeParameter(id + ".frequencies", 2, 0.5, 0.0, 1.0, writer);
        writer.writeCloseTag(BinaryCovarionModelParser.FREQUENCIES);

        writer.writeOpenTag(BinaryCovarionModelParser.HIDDEN_FREQUENCIES);
        writeParameter(id + ".hfrequencies", 2, 0.5, 0.0, 1.0, writer);
        writer.writeCloseTag(BinaryCovarionModelParser.HIDDEN_FREQUENCIES);

        writer.writeOpenTag(BinaryCovarionModelParser.ALPHA);
        writeParameter(id + ".alpha", writer);
        writer.writeCloseTag(BinaryCovarionModelParser.ALPHA);

        writer.writeOpenTag(BinaryCovarionModelParser.SWITCHING_RATE);
        writeParameter(id + ".s", writer);
        writer.writeCloseTag(BinaryCovarionModelParser.SWITCHING_RATE);

        writer.writeCloseTag(BinaryCovarionModelParser.COVARION_MODEL);
    }

    /**
     * Write the site model XML block.
     *
     * @param writer the writer
     */
    public void writeSiteModel(XMLWriter writer) {

        switch (dataType) {
            case DataType.NUCLEOTIDES:
                if (codonHeteroPattern != null) {
                    for (int i = 1; i <= partitionCount; i++) {
                        writeNucSiteModel(i, writer);
                    }
                    writer.println();
                    writer.writeOpenTag(CompoundParameterParser.COMPOUND_PARAMETER, new Attribute[]{new Attribute.Default<String>(XMLParser.ID, "allMus")});
                    for (int i = 1; i <= partitionCount; i++) {
                        writer.writeTag(ParameterParser.PARAMETER,
                                new Attribute[]{new Attribute.Default<String>(XMLParser.IDREF, SiteModel.SITE_MODEL + i + ".mu")}, true);
                    }
                    writer.writeCloseTag(CompoundParameterParser.COMPOUND_PARAMETER);
                } else {
                    writeNucSiteModel(-1, writer);
                }
                break;

            case DataType.AMINO_ACIDS:
                writeAASiteModel(writer);
                break;

            case DataType.TWO_STATES:
            case DataType.COVARION:
                writeTwoStateSiteModel(writer);
                break;

            default:
                throw new IllegalArgumentException("Unknown data type");
        }
    }

    /**
     * Write the nucleotide site model XML block.
     *
     * @param num    the model number
     * @param writer the writer
     */
    public void writeNucSiteModel(int num, XMLWriter writer) {

        String id = SiteModel.SITE_MODEL;
        if (num > 0) {
            id += Integer.toString(num);
        }

        writer.writeComment("site model");
        writer.writeOpenTag(GammaSiteModel.SITE_MODEL, new Attribute[]{new Attribute.Default<String>(XMLParser.ID, id)});


        writer.writeOpenTag(GammaSiteModelParser.SUBSTITUTION_MODEL);

        if (unlinkedSubstitutionModel) {
            switch (nucSubstitutionModel) {
                // JC cannot be unlinked because it has no parameters
                case JC:
                    writer.writeTag(NucModelType.HKY.getXMLName(), new Attribute.Default<String>(XMLParser.IDREF, "jc"), true);
                    break;
                case HKY:
                    writer.writeTag(NucModelType.HKY.getXMLName(), new Attribute.Default<String>(XMLParser.IDREF, "hky" + num), true);
                    break;
                case GTR:
                    writer.writeTag(GTRParser.GTR_MODEL, new Attribute.Default<String>(XMLParser.IDREF, "gtr" + num), true);
                    break;
                default:
                    throw new IllegalArgumentException("Unknown substitution model.");
            }
        } else {
            switch (nucSubstitutionModel) {
                case JC:
                    writer.writeTag(NucModelType.HKY.getXMLName(), new Attribute.Default<String>(XMLParser.IDREF, "jc"), true);
                    break;
                case HKY:
                    writer.writeTag(NucModelType.HKY.getXMLName(), new Attribute.Default<String>(XMLParser.IDREF, "hky"), true);
                    break;
                case GTR:
                    writer.writeTag(GTRParser.GTR_MODEL, new Attribute.Default<String>(XMLParser.IDREF, "gtr"), true);
                    break;
                default:
                    throw new IllegalArgumentException("Unknown substitution model.");
            }
        }
        writer.writeCloseTag(GammaSiteModelParser.SUBSTITUTION_MODEL);

        if (num != -1) {
            writer.writeOpenTag(GammaSiteModelParser.RELATIVE_RATE);
            writeParameter(id + ".mu", writer);
            writer.writeCloseTag(GammaSiteModelParser.RELATIVE_RATE);
        } else {
//            The actual mutation rate is now in the BranchRateModel so relativeRate can be missing
        }

        if (gammaHetero) {
            writer.writeOpenTag(GammaSiteModelParser.GAMMA_SHAPE, new Attribute.Default<String>(GammaSiteModelParser.GAMMA_CATEGORIES, "" + gammaCategories));
            if (num == -1 || unlinkedHeterogeneityModel) {
                writeParameter(id + ".alpha", writer);
            } else {
                // multiple partitions but linked heterogeneity
                if (num == 1) {
                    writeParameter(SiteModel.SITE_MODEL + "." + "alpha", writer);
                } else {
                    writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default<String>(XMLParser.IDREF, SiteModel.SITE_MODEL + "." + "alpha"), true);
                }
            }
            writer.writeCloseTag(GammaSiteModelParser.GAMMA_SHAPE);
        }

        if (invarHetero) {
            writer.writeOpenTag(GammaSiteModelParser.PROPORTION_INVARIANT);
            if (num == -1 || unlinkedHeterogeneityModel) {
                writeParameter(id + ".pInv", writer);
            } else {
                // multiple partitions but linked heterogeneity
                if (num == 1) {
                    writeParameter(SiteModel.SITE_MODEL + "." + "pInv", writer);
                } else {
                    writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default<String>(XMLParser.IDREF, SiteModel.SITE_MODEL + "." + "pInv"), true);
                }
            }
            writer.writeCloseTag(GammaSiteModelParser.PROPORTION_INVARIANT);
        }

        writer.writeCloseTag(GammaSiteModel.SITE_MODEL);
    }

    /**
     * Write the two states site model XML block.
     *
     * @param writer the writer
     */
    public void writeTwoStateSiteModel(XMLWriter writer) {

        String id = SiteModel.SITE_MODEL;

        writer.writeComment("site model");
        writer.writeOpenTag(GammaSiteModel.SITE_MODEL, new Attribute[]{new Attribute.Default<String>(XMLParser.ID, id)});


        writer.writeOpenTag(GammaSiteModelParser.SUBSTITUTION_MODEL);

        switch (binarySubstitutionModel) {
            case BIN_SIMPLE:
                //writer.writeTag(dr.evomodel.substmodel.GeneralSubstitutionModel.GENERAL_SUBSTITUTION_MODEL, new Attribute.Default<String>(XMLParser.IDREF, "bsimple"), true);
                writer.writeTag(BinarySubstitutionModelParser.BINARY_SUBSTITUTION_MODEL, new Attribute.Default<String>(XMLParser.IDREF, "bsimple"), true);
                break;
            case BIN_COVARION:
                writer.writeTag(BinaryCovarionModelParser.COVARION_MODEL, new Attribute.Default<String>(XMLParser.IDREF, "bcov"), true);
                break;
            default:
                throw new IllegalArgumentException("Unknown substitution model.");
        }

        writer.writeCloseTag(GammaSiteModelParser.SUBSTITUTION_MODEL);

        if (gammaHetero) {
            writer.writeOpenTag(GammaSiteModelParser.GAMMA_SHAPE, new Attribute.Default<String>(GammaSiteModelParser.GAMMA_CATEGORIES, "" + gammaCategories));
            writeParameter(id + ".alpha", writer);
            writer.writeCloseTag(GammaSiteModelParser.GAMMA_SHAPE);
        }

        if (invarHetero) {
            writer.writeOpenTag(GammaSiteModelParser.PROPORTION_INVARIANT);
            writeParameter(id + ".pInv", writer);
            writer.writeCloseTag(GammaSiteModelParser.PROPORTION_INVARIANT);
        }

        writer.writeCloseTag(GammaSiteModel.SITE_MODEL);
    }


    /**
     * Write the AA site model XML block.
     *
     * @param writer the writer
     */
    public void writeAASiteModel(XMLWriter writer) {

        writer.writeComment("site model");
        writer.writeOpenTag(GammaSiteModel.SITE_MODEL, new Attribute[]{
                new Attribute.Default<String>(XMLParser.ID, SiteModel.SITE_MODEL)});


        writer.writeOpenTag(GammaSiteModelParser.SUBSTITUTION_MODEL);
        writer.writeTag(EmpiricalAminoAcidModelParser.EMPIRICAL_AMINO_ACID_MODEL,
                new Attribute.Default<String>(XMLParser.IDREF, "aa"), true);
        writer.writeCloseTag(GammaSiteModelParser.SUBSTITUTION_MODEL);

//            The actual mutation rate is now in the BranchRateModel so relativeRate can be missing

        if (gammaHetero) {
            writer.writeOpenTag(GammaSiteModelParser.GAMMA_SHAPE, new Attribute.Default<String>(GammaSiteModelParser.GAMMA_CATEGORIES, "" + gammaCategories));
            writeParameter(SiteModel.SITE_MODEL + "." + "alpha", writer);
            writer.writeCloseTag(GammaSiteModelParser.GAMMA_SHAPE);
        }

        if (invarHetero) {
            writer.writeOpenTag(GammaSiteModelParser.PROPORTION_INVARIANT);
            writeParameter(SiteModel.SITE_MODEL + "." + "pInv", writer);
            writer.writeCloseTag(GammaSiteModelParser.PROPORTION_INVARIANT);
        }

        writer.writeCloseTag(GammaSiteModel.SITE_MODEL);
    }


    /**
     * Write the relaxed clock branch rates block.
     *
     * @param writer the writer
     */
    public void writeBranchRatesModel(XMLWriter writer) {
        if (clockModel == STRICT_CLOCK) {
            if (fixedSubstitutionRate) {
                fixParameter("clock.rate", meanSubstitutionRate);
            }

            writer.writeComment("The strict clock (Uniform rates across branches)");
            writer.writeOpenTag(
                    StrictClockBranchRatesParser.STRICT_CLOCK_BRANCH_RATES,
                    new Attribute[]{new Attribute.Default<String>(XMLParser.ID, BranchRateModel.BRANCH_RATES)}
            );
            writer.writeOpenTag("rate");

            writeParameter("clock.rate", writer);
            writer.writeCloseTag("rate");
            writer.writeCloseTag(StrictClockBranchRatesParser.STRICT_CLOCK_BRANCH_RATES);
        } else if (clockModel == RANDOM_LOCAL_CLOCK) {
            if (fixedSubstitutionRate) {
                fixParameter("clock.rate", meanSubstitutionRate);
            }

            writer.writeComment("The random local clock model (Drummond & Suchard, 2007)");
            writer.writeOpenTag(
                    RandomLocalClockModelParser.LOCAL_BRANCH_RATES,
                    new Attribute[]{
                            new Attribute.Default<String>(XMLParser.ID, BranchRateModel.BRANCH_RATES),
                            new Attribute.Default<String>("ratesAreMultipliers", "false")
                    }
            );
            writer.writeTag(TreeModel.TREE_MODEL, new Attribute.Default<String>(XMLParser.IDREF, "treeModel"), true);

            writer.writeOpenTag("rates");
            writer.writeTag("parameter", new Attribute.Default<String>(XMLParser.IDREF, LOCAL_CLOCK + "." + "rates"), true);
            writer.writeCloseTag("rates");

            writer.writeOpenTag("rateIndicator");
            writer.writeTag("parameter", new Attribute.Default<String>(XMLParser.IDREF, LOCAL_CLOCK + "." + "changes"), true);
            writer.writeCloseTag("rateIndicator");

            writer.writeOpenTag("clockRate");
            writeParameter("clock.rate", writer);
            writer.writeCloseTag("clockRate");

            writer.writeCloseTag(RandomLocalClockModelParser.LOCAL_BRANCH_RATES);

            writer.writeText("");
            writer.writeOpenTag(
                    SumStatisticParser.SUM_STATISTIC,
                    new Attribute[]{
                            new Attribute.Default<String>(XMLParser.ID, "rateChanges"),
                            new Attribute.Default<String>("name", "rateChangeCount"),
                            new Attribute.Default<String>("elementwise", "true"),
                    }
            );
            writer.writeTag("parameter", new Attribute.Default<String>(XMLParser.IDREF, LOCAL_CLOCK + "." + "changes"), true);
            writer.writeCloseTag(SumStatisticParser.SUM_STATISTIC);

            writer.writeText("");

            writer.writeOpenTag(
                    RateStatisticParser.RATE_STATISTIC,
                    new Attribute[]{
                            new Attribute.Default<String>(XMLParser.ID, "meanRate"),
                            new Attribute.Default<String>("name", "meanRate"),
                            new Attribute.Default<String>("mode", "mean"),
                            new Attribute.Default<String>("internal", "true"),
                            new Attribute.Default<String>("external", "true")
                    }
            );
            writer.writeTag(TreeModel.TREE_MODEL, new Attribute.Default<String>(XMLParser.IDREF, "treeModel"), true);
            writer.writeTag(RandomLocalClockModelParser.LOCAL_BRANCH_RATES, new Attribute.Default<String>(XMLParser.IDREF, BranchRateModel.BRANCH_RATES), true);
            writer.writeCloseTag(RateStatisticParser.RATE_STATISTIC);

            writer.writeText("");
            writer.writeOpenTag(
                    RateStatisticParser.RATE_STATISTIC,
                    new Attribute[]{
                            new Attribute.Default<String>(XMLParser.ID, RateStatisticParser.COEFFICIENT_OF_VARIATION),
                            new Attribute.Default<String>("name", RateStatisticParser.COEFFICIENT_OF_VARIATION),
                            new Attribute.Default<String>("mode", RateStatisticParser.COEFFICIENT_OF_VARIATION),
                            new Attribute.Default<String>("internal", "true"),
                            new Attribute.Default<String>("external", "true")
                    }
            );
            writer.writeTag(TreeModel.TREE_MODEL, new Attribute.Default<String>(XMLParser.IDREF, "treeModel"), true);
            writer.writeTag(RandomLocalClockModelParser.LOCAL_BRANCH_RATES, new Attribute.Default<String>(XMLParser.IDREF, BranchRateModel.BRANCH_RATES), true);
            writer.writeCloseTag(RateStatisticParser.RATE_STATISTIC);

            writer.writeText("");
            writer.writeOpenTag(
                    RateCovarianceStatisticParser.RATE_COVARIANCE_STATISTIC,
                    new Attribute[]{
                            new Attribute.Default<String>(XMLParser.ID, "covariance"),
                            new Attribute.Default<String>("name", "covariance")
                    }
            );
            writer.writeTag(TreeModel.TREE_MODEL, new Attribute.Default<String>(XMLParser.IDREF, "treeModel"), true);
            writer.writeTag(RandomLocalClockModelParser.LOCAL_BRANCH_RATES, new Attribute.Default<String>(XMLParser.IDREF, BranchRateModel.BRANCH_RATES), true);
            writer.writeCloseTag(RateCovarianceStatisticParser.RATE_COVARIANCE_STATISTIC);

        } else {
            writer.writeComment("The uncorrelated relaxed clock (Drummond, Ho, Phillips & Rambaut, 2006)");
            writer.writeOpenTag(
                    DiscretizedBranchRatesParser.DISCRETIZED_BRANCH_RATES,
                    new Attribute[]{new Attribute.Default<String>(XMLParser.ID, BranchRateModel.BRANCH_RATES)}
            );
            writer.writeTag(TreeModel.TREE_MODEL, new Attribute.Default<String>(XMLParser.IDREF, "treeModel"), true);
            writer.writeOpenTag("distribution");
            if (clockModel == UNCORRELATED_EXPONENTIAL) {
                if (fixedSubstitutionRate) {
                    fixParameter(UCED_MEAN, meanSubstitutionRate);
                }

                final String eModelName = ExponentialDistributionModel.EXPONENTIAL_DISTRIBUTION_MODEL;
                writer.writeOpenTag(eModelName);
                writer.writeOpenTag("mean");
                writeParameter(UCED_MEAN, writer);
                writer.writeCloseTag("mean");
                writer.writeCloseTag(eModelName);
            } else if (clockModel == UNCORRELATED_LOGNORMAL) {
                if (fixedSubstitutionRate) {
                    fixParameter(UCLD_MEAN, meanSubstitutionRate);
                }

                writer.writeOpenTag(LogNormalDistributionModelParser.LOGNORMAL_DISTRIBUTION_MODEL,
                        new Attribute.Default<String>(LogNormalDistributionModelParser.MEAN_IN_REAL_SPACE, "true"));
                writer.writeOpenTag("mean");
                writeParameter(UCLD_MEAN, writer);
                writer.writeCloseTag("mean");
                writer.writeOpenTag("stdev");
                writeParameter(UCLD_STDEV, writer);
                writer.writeCloseTag("stdev");
                writer.writeCloseTag(LogNormalDistributionModelParser.LOGNORMAL_DISTRIBUTION_MODEL);
            } else {
                throw new RuntimeException("Unrecognised relaxed clock model");
            }
            writer.writeCloseTag("distribution");
            writer.writeOpenTag("rateCategories");
            int categoryCount = (alignment.getSequenceCount() - 1) * 2;
            writeParameter("branchRates.categories", categoryCount, writer);
            writer.writeCloseTag("rateCategories");
            writer.writeCloseTag(DiscretizedBranchRatesParser.DISCRETIZED_BRANCH_RATES);

            writer.writeText("");
            writer.writeOpenTag(
                    RateStatisticParser.RATE_STATISTIC,
                    new Attribute[]{
                            new Attribute.Default<String>(XMLParser.ID, "meanRate"),
                            new Attribute.Default<String>("name", "meanRate"),
                            new Attribute.Default<String>("mode", "mean"),
                            new Attribute.Default<String>("internal", "true"),
                            new Attribute.Default<String>("external", "true")
                    }
            );
            writer.writeTag(TreeModel.TREE_MODEL, new Attribute.Default<String>(XMLParser.IDREF, "treeModel"), true);
            writer.writeTag(DiscretizedBranchRatesParser.DISCRETIZED_BRANCH_RATES, new Attribute.Default<String>(XMLParser.IDREF, BranchRateModel.BRANCH_RATES), true);
            writer.writeCloseTag(RateStatisticParser.RATE_STATISTIC);

            writer.writeText("");
            writer.writeOpenTag(
                    RateStatisticParser.RATE_STATISTIC,
                    new Attribute[]{
                            new Attribute.Default<String>(XMLParser.ID, RateStatisticParser.COEFFICIENT_OF_VARIATION),
                            new Attribute.Default<String>("name", RateStatisticParser.COEFFICIENT_OF_VARIATION),
                            new Attribute.Default<String>("mode", RateStatisticParser.COEFFICIENT_OF_VARIATION),
                            new Attribute.Default<String>("internal", "true"),
                            new Attribute.Default<String>("external", "true")
                    }
            );
            writer.writeTag(TreeModel.TREE_MODEL, new Attribute.Default<String>(XMLParser.IDREF, "treeModel"), true);
            writer.writeTag(DiscretizedBranchRatesParser.DISCRETIZED_BRANCH_RATES, new Attribute.Default<String>(XMLParser.IDREF, BranchRateModel.BRANCH_RATES), true);
            writer.writeCloseTag(RateStatisticParser.RATE_STATISTIC);

            writer.writeText("");
            writer.writeOpenTag(
                    RateCovarianceStatisticParser.RATE_COVARIANCE_STATISTIC,
                    new Attribute[]{
                            new Attribute.Default<String>(XMLParser.ID, "covariance"),
                            new Attribute.Default<String>("name", "covariance")
                    }
            );
            writer.writeTag(TreeModel.TREE_MODEL, new Attribute.Default<String>(XMLParser.IDREF, "treeModel"), true);
            writer.writeTag(DiscretizedBranchRatesParser.DISCRETIZED_BRANCH_RATES, new Attribute.Default<String>(XMLParser.IDREF, BranchRateModel.BRANCH_RATES), true);
            writer.writeCloseTag(RateCovarianceStatisticParser.RATE_COVARIANCE_STATISTIC);
        }
    }

    /**
     * Write the prior on node heights (coalescent or speciational models)
     *
     * @param writer the writer
     */
    public void writeNodeHeightPrior(XMLWriter writer) {
        if (nodeHeightPrior == YULE || nodeHeightPrior == BIRTH_DEATH) {
            // generate a speciational process

            writer.writeOpenTag(
                    SpeciationLikelihoodParser.SPECIATION_LIKELIHOOD,
                    new Attribute[]{
                            new Attribute.Default<String>(XMLParser.ID, "speciation")
                    }
            );

            // write pop size socket
            writer.writeOpenTag(SpeciationLikelihoodParser.MODEL);
            writeNodeHeightPriorModelRef(writer);
            writer.writeCloseTag(SpeciationLikelihoodParser.MODEL);
            writer.writeOpenTag(SpeciationLikelihoodParser.TREE);
            writer.writeTag(TreeModel.TREE_MODEL, new Attribute.Default<String>(XMLParser.IDREF, "treeModel"), true);
            writer.writeCloseTag(SpeciationLikelihoodParser.TREE);

            writer.writeCloseTag(SpeciationLikelihoodParser.SPECIATION_LIKELIHOOD);

        } else if (nodeHeightPrior == SKYLINE) {
            // generate a Bayesian skyline plot

            writer.writeOpenTag(
                    BayesianSkylineLikelihoodParser.SKYLINE_LIKELIHOOD,
                    new Attribute[]{
                            new Attribute.Default<String>(XMLParser.ID, "skyline"),
                            new Attribute.Default<String>("linear", skylineModel == LINEAR_SKYLINE ? "true" : "false")
                    }
            );

            // write pop size socket
            writer.writeOpenTag(BayesianSkylineLikelihoodParser.POPULATION_SIZES);
            if (skylineModel == LINEAR_SKYLINE) {
                writeParameter("skyline.popSize", skylineGroupCount + 1, writer);
            } else {
                writeParameter("skyline.popSize", skylineGroupCount, writer);
            }
            writer.writeCloseTag(BayesianSkylineLikelihoodParser.POPULATION_SIZES);

            // write group size socket
            writer.writeOpenTag(BayesianSkylineLikelihoodParser.GROUP_SIZES);
            writeParameter("skyline.groupSize", skylineGroupCount, writer);
            writer.writeCloseTag(BayesianSkylineLikelihoodParser.GROUP_SIZES);

            writer.writeOpenTag(CoalescentLikelihoodParser.POPULATION_TREE);
            writer.writeTag(TreeModel.TREE_MODEL, new Attribute.Default<String>(XMLParser.IDREF, "treeModel"), true);
            writer.writeCloseTag(CoalescentLikelihoodParser.POPULATION_TREE);

            writer.writeCloseTag(BayesianSkylineLikelihoodParser.SKYLINE_LIKELIHOOD);
        } else if (nodeHeightPrior == EXTENDED_SKYLINE) {
            final String tagName = VariableDemographicModelParser.MODEL_NAME;

            writer.writeOpenTag(
                    tagName,
                    new Attribute[]{
                            new Attribute.Default<String>(XMLParser.ID, VariableDemographicModelParser.demoElementName),
                            new Attribute.Default<String>(VariableDemographicModelParser.TYPE, extendedSkylineModel)
                    }
            );

            writer.writeOpenTag(VariableDemographicModelParser.POPULATION_SIZES);
            final int nTax = taxonList.getTaxonCount();
            final int nPops = nTax - (extendedSkylineModel.equals(VariableDemographicModel.Type.STEPWISE.toString()) ? 1 : 0);
            writeParameter(VariableDemographicModelParser.demoElementName + ".popSize", nPops, writer);
            writer.writeCloseTag(VariableDemographicModelParser.POPULATION_SIZES);

            writer.writeOpenTag(VariableDemographicModelParser.INDICATOR_PARAMETER);
            writeParameter(VariableDemographicModelParser.demoElementName + ".indicators", nPops - 1, writer);
            writer.writeCloseTag(VariableDemographicModelParser.INDICATOR_PARAMETER);

            writer.writeOpenTag(VariableDemographicModelParser.POPULATION_TREES);

            writer.writeOpenTag(VariableDemographicModelParser.POP_TREE);
            writer.writeTag(TreeModel.TREE_MODEL, new Attribute.Default<String>(XMLParser.IDREF, "treeModel"), true);
            writer.writeCloseTag(VariableDemographicModelParser.POP_TREE);

            writer.writeCloseTag(VariableDemographicModelParser.POPULATION_TREES);

            writer.writeCloseTag(tagName);

            writer.writeOpenTag(CoalescentLikelihoodParser.COALESCENT_LIKELIHOOD, new Attribute.Default<String>(XMLParser.ID, "coalescent"));
            writer.writeOpenTag(CoalescentLikelihoodParser.MODEL);
            writer.writeTag(tagName, new Attribute.Default<String>(XMLParser.IDREF, VariableDemographicModelParser.demoElementName), true);
            writer.writeCloseTag(CoalescentLikelihoodParser.MODEL);
            writer.writeComment("Take population Tree from demographic");
            writer.writeCloseTag(CoalescentLikelihoodParser.COALESCENT_LIKELIHOOD);

            writer.writeOpenTag(SumStatisticParser.SUM_STATISTIC,
                    new Attribute[]{
                            new Attribute.Default<String>(XMLParser.ID, VariableDemographicModelParser.demoElementName + ".populationSizeChanges"),
                            new Attribute.Default<String>("elementwise", "true")
                    });
            writer.writeTag(ParameterParser.PARAMETER,
                    new Attribute.Default<String>(XMLParser.IDREF, VariableDemographicModelParser.demoElementName + ".indicators"), true);
            writer.writeCloseTag(SumStatisticParser.SUM_STATISTIC);
            writer.writeOpenTag(ExponentialDistributionModel.EXPONENTIAL_DISTRIBUTION_MODEL,
                    new Attribute[]{
                            new Attribute.Default<String>(XMLParser.ID, VariableDemographicModelParser.demoElementName + ".populationMeanDist")
                            //,new Attribute.Default<String>("elementwise", "true")
                    });
            writer.writeOpenTag(DistributionModelParser.MEAN);
            writer.writeTag(ParameterParser.PARAMETER,
                    new Attribute[]{
                            new Attribute.Default<String>(XMLParser.ID, VariableDemographicModelParser.demoElementName + ".populationMean"),
                            new Attribute.Default<String>("value", "1")}, true);
            writer.writeCloseTag(DistributionModelParser.MEAN);
            writer.writeCloseTag(ExponentialDistributionModel.EXPONENTIAL_DISTRIBUTION_MODEL);

        } else {
            // generate a coalescent process

            writer.writeOpenTag(
                    CoalescentLikelihoodParser.COALESCENT_LIKELIHOOD,
                    new Attribute[]{new Attribute.Default<String>(XMLParser.ID, "coalescent")}
            );
            writer.writeOpenTag(CoalescentLikelihoodParser.MODEL);
            writeNodeHeightPriorModelRef(writer);
            writer.writeCloseTag(CoalescentLikelihoodParser.MODEL);
            writer.writeOpenTag(CoalescentLikelihoodParser.POPULATION_TREE);
            writer.writeTag(TreeModel.TREE_MODEL, new Attribute.Default<String>(XMLParser.IDREF, "treeModel"), true);
            writer.writeCloseTag(CoalescentLikelihoodParser.POPULATION_TREE);
            writer.writeCloseTag(CoalescentLikelihoodParser.COALESCENT_LIKELIHOOD);
        }
    }

    /**
     * Write the boolean likelihood
     *
     * @param writer the writer
     */
    public void writeBooleanLikelihood(XMLWriter writer) {
        writer.writeOpenTag(
                BooleanLikelihoodParser.BOOLEAN_LIKELIHOOD,
                new Attribute[]{new Attribute.Default<String>(XMLParser.ID, "booleanLikelihood1")}
        );
        writer.writeOpenTag(
                TestStatisticParser.TEST_STATISTIC,
                new Attribute[]{
                        new Attribute.Default<String>(XMLParser.ID, "test1"),
                        new Attribute.Default<String>("name", "test1")
                }
        );
        writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default<String>(XMLParser.IDREF, "logistic.t50"), true);
        writer.writeOpenTag("lessThan");
        writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default<String>(XMLParser.IDREF, "treeModel.rootHeight"), true);
        writer.writeCloseTag("lessThan");
        writer.writeCloseTag(TestStatisticParser.TEST_STATISTIC);
        writer.writeCloseTag(BooleanLikelihoodParser.BOOLEAN_LIKELIHOOD);
    }

    public void writeExponentialMarkovLikelihood(XMLWriter writer) {
        writer.writeOpenTag(
                ExponentialMarkovModel.EXPONENTIAL_MARKOV_MODEL,
                new Attribute[]{new Attribute.Default<String>(XMLParser.ID, "eml1"),
                        new Attribute.Default<String>("jeffreys", "true")}
        );
        writer.writeOpenTag(ExponentialMarkovModelParser.CHAIN_PARAMETER);
        writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default<String>(XMLParser.IDREF, "skyline.popSize"), true);
        writer.writeCloseTag(ExponentialMarkovModelParser.CHAIN_PARAMETER);
        writer.writeCloseTag(ExponentialMarkovModel.EXPONENTIAL_MARKOV_MODEL);
    }


    /**
     * Write the tree likelihood XML block.
     *
     * @param writer the writer
     */
    public void writeTreeLikelihood(XMLWriter writer) {

        boolean nucs = alignment.getDataType() == Nucleotides.INSTANCE;
        if (nucs && codonHeteroPattern != null) {
            for (int i = 1; i <= partitionCount; i++) {
                writeTreeLikelihood(i, writer);
            }
        } else {
            writeTreeLikelihood(-1, writer);
        }
    }


    /**
     * Determine and return the datatype description for these beast options
     * note that the datatype in XML may differ from the actual datatype
     *
     * @return description
     */

    private Boolean useAmbiguities() {
        Boolean useAmbiguities = false;

        switch (dataType) {
            case DataType.TWO_STATES:
            case DataType.COVARION:

                switch (binarySubstitutionModel) {
                    case BIN_COVARION:
                        useAmbiguities = true;
                        break;

                    default:
                }
                break;

            default:
                useAmbiguities = false;
        }

        return useAmbiguities;
    }

    /**
     * Write the tree likelihood XML block.
     *
     * @param num    the likelihood number
     * @param writer the writer
     */
    public void writeTreeLikelihood(int num, XMLWriter writer) {

        if (num > 0) {
            writer.writeOpenTag(
                    TreeLikelihoodParser.TREE_LIKELIHOOD,
                    new Attribute[]{
                            new Attribute.Default<String>(XMLParser.ID, TreeLikelihoodParser.TREE_LIKELIHOOD + num),
                            new Attribute.Default<Boolean>(TreeLikelihoodParser.USE_AMBIGUITIES, useAmbiguities())}
            );
            if (codonHeteroPattern.equals("112")) {
                if (num == 1) {
                    writer.writeTag(SitePatternsParser.PATTERNS, new Attribute[]{new Attribute.Default<String>(XMLParser.IDREF, "patterns1+2")}, true);
                } else {
                    writer.writeTag(SitePatternsParser.PATTERNS, new Attribute[]{new Attribute.Default<String>(XMLParser.IDREF, "patterns3")}, true);
                }
            } else {
                writer.writeTag(SitePatternsParser.PATTERNS, new Attribute[]{new Attribute.Default<String>(XMLParser.IDREF, SitePatternsParser.PATTERNS + num)}, true);
            }
            writer.writeTag(TreeModel.TREE_MODEL, new Attribute[]{new Attribute.Default<String>(XMLParser.IDREF, "treeModel")}, true);
            writer.writeTag(GammaSiteModel.SITE_MODEL, new Attribute[]{new Attribute.Default<String>(XMLParser.IDREF, SiteModel.SITE_MODEL + num)}, true);
        } else {
            writer.writeOpenTag(
                    TreeLikelihoodParser.TREE_LIKELIHOOD,
                    new Attribute[]{
                            new Attribute.Default<String>(XMLParser.ID, TreeLikelihoodParser.TREE_LIKELIHOOD),
                            new Attribute.Default<Boolean>(TreeLikelihoodParser.USE_AMBIGUITIES, useAmbiguities())
                    }
            );
            writer.writeTag(SitePatternsParser.PATTERNS, new Attribute[]{new Attribute.Default<String>(XMLParser.IDREF, SitePatternsParser.PATTERNS)}, true);
            writer.writeTag(TreeModel.TREE_MODEL, new Attribute[]{new Attribute.Default<String>(XMLParser.IDREF, "treeModel")}, true);
            writer.writeTag(GammaSiteModel.SITE_MODEL, new Attribute[]{new Attribute.Default<String>(XMLParser.IDREF, SiteModel.SITE_MODEL)}, true);
        }
        if (clockModel == STRICT_CLOCK) {
            writer.writeTag(StrictClockBranchRatesParser.STRICT_CLOCK_BRANCH_RATES, new Attribute[]{new Attribute.Default<String>(XMLParser.IDREF, BranchRateModel.BRANCH_RATES)}, true);
        } else {
            writer.writeTag(DiscretizedBranchRatesParser.DISCRETIZED_BRANCH_RATES, new Attribute[]{new Attribute.Default<String>(XMLParser.IDREF, BranchRateModel.BRANCH_RATES)}, true);
        }

        writer.writeCloseTag(TreeLikelihoodParser.TREE_LIKELIHOOD);
    }

    /**
     * Generate tmrca statistics
     *
     * @param writer the writer
     */
    public void writeTMRCAStatistics(XMLWriter writer) {

        writer.writeText("");
        for (Taxa taxa : taxonSets) {
            writer.writeOpenTag(
                    TMRCAStatisticParser.TMRCA_STATISTIC,
                    new Attribute[]{
                            new Attribute.Default<String>(XMLParser.ID, "tmrca(" + taxa.getId() + ")"),
                    }
            );
            writer.writeOpenTag(TMRCAStatisticParser.MRCA);
            writer.writeTag(TaxaParser.TAXA, new Attribute[]{new Attribute.Default<String>(XMLParser.IDREF, taxa.getId())}, true);
            writer.writeCloseTag(TMRCAStatisticParser.MRCA);
            writer.writeTag(TreeModel.TREE_MODEL, new Attribute[]{new Attribute.Default<String>(XMLParser.IDREF, "treeModel")}, true);
            writer.writeCloseTag(TMRCAStatisticParser.TMRCA_STATISTIC);

            if (taxonSetsMono.get(taxa)) {
                writer.writeOpenTag(
                        MonophylyStatisticParser.MONOPHYLY_STATISTIC,
                        new Attribute[]{
                                new Attribute.Default<String>(XMLParser.ID, "monophyly(" + taxa.getId() + ")"),
                        });
                writer.writeOpenTag(MonophylyStatisticParser.MRCA);
                writer.writeTag(TaxaParser.TAXA, new Attribute[]{new Attribute.Default<String>(XMLParser.IDREF, taxa.getId())}, true);
                writer.writeCloseTag(MonophylyStatisticParser.MRCA);
                writer.writeTag(TreeModel.TREE_MODEL, new Attribute[]{new Attribute.Default<String>(XMLParser.IDREF, "treeModel")}, true);
                writer.writeCloseTag(MonophylyStatisticParser.MONOPHYLY_STATISTIC);
            }
        }
    }

    /**
     * Write the operator schedule XML block.
     *
     * @param operators the list of operators
     * @param writer    the writer
     */
    public void writeOperatorSchedule(ArrayList<Operator> operators, XMLWriter writer) {
        writer.writeOpenTag(
                SimpleOperatorScheduleParser.OPERATOR_SCHEDULE,
                new Attribute[]{new Attribute.Default<String>(XMLParser.ID, "operators")}
        );

        for (Operator operator : operators) {
            if (operator.weight > 0. && operator.inUse)
                writeOperator(operator, writer);
        }

        writer.writeCloseTag(SimpleOperatorScheduleParser.OPERATOR_SCHEDULE);
    }

    private void writeOperator(Operator operator, XMLWriter writer) {
        if (operator.type.equals(SCALE)) {
            writeScaleOperator(operator, writer);
        } else if (operator.type.equals(RANDOM_WALK)) {
            writeRandomWalkOperator(operator, writer);
        } else if (operator.type.equals(INTEGER_RANDOM_WALK)) {
            writeIntegerRandomWalkOperator(operator, writer);
        } else if (operator.type.equals(UP_DOWN)) {
            writeUpDownOperator(operator, writer);
        } else if (operator.type.equals(SCALE_ALL)) {
            writeScaleAllOperator(operator, writer);
        } else if (operator.type.equals(CENTERED_SCALE)) {
            writeCenteredOperator(operator, writer);
        } else if (operator.type.equals(DELTA_EXCHANGE)) {
            writeDeltaOperator(operator, writer);
        } else if (operator.type.equals(INTEGER_DELTA_EXCHANGE)) {
            writeIntegerDeltaOperator(operator, writer);
        } else if (operator.type.equals(SWAP)) {
            writeSwapOperator(operator, writer);
        } else if (operator.type.equals(BITFLIP)) {
            writeBitFlipOperator(operator, writer);
        } else if (operator.type.equals(TREE_BIT_MOVE)) {
            writeTreeBitMoveOperator(operator, writer);
        } else if (operator.type.equals(UNIFORM)) {
            writeUniformOperator(operator, writer);
        } else if (operator.type.equals(INTEGER_UNIFORM)) {
            writeIntegerUniformOperator(operator, writer);
        } else if (operator.type.equals(SUBTREE_SLIDE)) {
            writeSubtreeSlideOperator(operator, writer);
        } else if (operator.type.equals(NARROW_EXCHANGE)) {
            writeNarrowExchangeOperator(operator, writer);
        } else if (operator.type.equals(WIDE_EXCHANGE)) {
            writeWideExchangeOperator(operator, writer);
        } else if (operator.type.equals(WILSON_BALDING)) {
            writeWilsonBaldingOperator(operator, writer);
        } else if (operator.type.equals(SAMPLE_NONACTIVE)) {
            writeSampleNonActiveOperator(operator, writer);
        } else if (operator.type.equals(SCALE_WITH_INDICATORS)) {
            writeScaleWithIndicatorsOperator(operator, writer);
        }
    }

    private Attribute getRef(String name) {
        return new Attribute.Default<String>(XMLParser.IDREF, name);
    }

    private void writeParameterRefByName(XMLWriter writer, String name) {
        writer.writeTag(ParameterParser.PARAMETER, getRef(name), true);
    }

    private void writeParameter1Ref(XMLWriter writer, Operator operator) {
        writeParameterRefByName(writer, operator.parameter1.getName());
    }

    private void writeScaleOperator(Operator operator, XMLWriter writer) {
        writer.writeOpenTag(
                ScaleOperatorParser.SCALE_OPERATOR,
                new Attribute[]{
                        new Attribute.Default<Double>(ScaleOperatorParser.SCALE_FACTOR, operator.tuning),
                        new Attribute.Default<Double>("weight", operator.weight),
                });
        writeParameter1Ref(writer, operator);
        writer.writeCloseTag(ScaleOperatorParser.SCALE_OPERATOR);
    }

    private void writeRandomWalkOperator(Operator operator, XMLWriter writer) {
        writer.writeOpenTag(
                "randomWalkOperator",
                new Attribute[]{
                        new Attribute.Default<Double>("windowSize", operator.tuning),
                        new Attribute.Default<Double>("weight", operator.weight)
                });
        writeParameter1Ref(writer, operator);
        writer.writeCloseTag("randomWalkOperator");
    }

    private void writeIntegerRandomWalkOperator(Operator operator, XMLWriter writer) {
        writer.writeOpenTag(
                RandomWalkIntegerOperatorParser.RANDOM_WALK_INTEGER_OPERATOR,
                new Attribute[]{
                        new Attribute.Default<Double>("windowSize", operator.tuning),
                        new Attribute.Default<Double>("weight", operator.weight)
                });
        writeParameter1Ref(writer, operator);
        writer.writeCloseTag(RandomWalkIntegerOperatorParser.RANDOM_WALK_INTEGER_OPERATOR);
    }

    private void writeScaleAllOperator(Operator operator, XMLWriter writer) {
        writer.writeOpenTag(
                ScaleOperatorParser.SCALE_OPERATOR,
                new Attribute[]{
                        new Attribute.Default<Double>(ScaleOperatorParser.SCALE_FACTOR, operator.tuning),
                        new Attribute.Default<String>(ScaleOperatorParser.SCALE_ALL, "true"),
                        new Attribute.Default<Double>("weight", operator.weight),
                });
        writer.writeOpenTag(CompoundParameterParser.COMPOUND_PARAMETER);
        writeParameter1Ref(writer, operator);
        writer.writeTag(ParameterParser.PARAMETER, new Attribute[]{new Attribute.Default<String>(XMLParser.IDREF, operator.parameter2.getName())}, true);
        writer.writeCloseTag(CompoundParameterParser.COMPOUND_PARAMETER);
        writer.writeCloseTag(ScaleOperatorParser.SCALE_OPERATOR);
    }

    private void writeUpDownOperator(Operator operator, XMLWriter writer) {
        writer.writeOpenTag(UpDownOperatorParser.UP_DOWN_OPERATOR,
                new Attribute[]{
                        new Attribute.Default<Double>(UpDownOperatorParser.SCALE_FACTOR, operator.tuning),
                        new Attribute.Default<Double>("weight", operator.weight),
                }
        );

        writer.writeOpenTag(UpDownOperatorParser.UP);
        writeParameter1Ref(writer, operator);
        writer.writeCloseTag(UpDownOperatorParser.UP);

        writer.writeOpenTag(UpDownOperatorParser.DOWN);
        writer.writeTag(ParameterParser.PARAMETER, new Attribute[]{new Attribute.Default<String>(XMLParser.IDREF, operator.parameter2.getName())}, true);
        writer.writeCloseTag(UpDownOperatorParser.DOWN);

        writer.writeCloseTag(UpDownOperatorParser.UP_DOWN_OPERATOR);
    }

    private void writeCenteredOperator(Operator operator, XMLWriter writer) {
        writer.writeOpenTag(CenteredScaleOperatorParser.CENTERED_SCALE,
                new Attribute[]{
                        new Attribute.Default<Double>(CenteredScaleOperatorParser.SCALE_FACTOR, operator.tuning),
                        new Attribute.Default<Double>("weight", operator.weight),
                }
        );
        writeParameter1Ref(writer, operator);
        writer.writeCloseTag(CenteredScaleOperatorParser.CENTERED_SCALE);
    }

    private void writeDeltaOperator(Operator operator, XMLWriter writer) {
        partitionCount = getPartionCount(codonHeteroPattern);

        if (operator.name.equals("Relative rates") && codonHeteroPattern.equals("112")) {
            writer.writeOpenTag(DeltaExchangeOperatorParser.DELTA_EXCHANGE,
                    new Attribute[]{
                            new Attribute.Default<Double>(DeltaExchangeOperatorParser.DELTA, operator.tuning),
                            new Attribute.Default<String>(DeltaExchangeOperatorParser.PARAMETER_WEIGHTS, "2 1"),
                            new Attribute.Default<Double>("weight", operator.weight),
                    }
            );
        } else {
            writer.writeOpenTag(DeltaExchangeOperatorParser.DELTA_EXCHANGE,
                    new Attribute[]{
                            new Attribute.Default<Double>(DeltaExchangeOperatorParser.DELTA, operator.tuning),
                            new Attribute.Default<Double>("weight", operator.weight),
                    }
            );
        }

        writeParameter1Ref(writer, operator);
        writer.writeCloseTag(DeltaExchangeOperatorParser.DELTA_EXCHANGE);
    }

    private void writeIntegerDeltaOperator(Operator operator, XMLWriter writer) {
        writer.writeOpenTag(DeltaExchangeOperatorParser.DELTA_EXCHANGE,
                new Attribute[]{
                        new Attribute.Default<String>(DeltaExchangeOperatorParser.DELTA, Integer.toString((int) operator.tuning)),
                        new Attribute.Default<String>("integer", "true"),
                        new Attribute.Default<Double>("weight", operator.weight),
                        new Attribute.Default<String>("autoOptimize", "false")
                }
        );
        writeParameter1Ref(writer, operator);
        writer.writeCloseTag(DeltaExchangeOperatorParser.DELTA_EXCHANGE);
    }

    private void writeSwapOperator(Operator operator, XMLWriter writer) {
        writer.writeOpenTag(SwapOperatorParser.SWAP_OPERATOR,
                new Attribute[]{
                        new Attribute.Default<String>("size", Integer.toString((int) operator.tuning)),
                        new Attribute.Default<Double>("weight", operator.weight),
                        new Attribute.Default<String>("autoOptimize", "false")
                }
        );
        writeParameter1Ref(writer, operator);
        writer.writeCloseTag(SwapOperatorParser.SWAP_OPERATOR);
    }

    private void writeBitFlipOperator(Operator operator, XMLWriter writer) {
        writer.writeOpenTag(BitFlipOperatorParser.BIT_FLIP_OPERATOR,
                new Attribute[]{
                        new Attribute.Default<Double>("weight", operator.weight),
                }
        );
        writeParameter1Ref(writer, operator);
        writer.writeCloseTag(BitFlipOperatorParser.BIT_FLIP_OPERATOR);
    }

    private void writeTreeBitMoveOperator(Operator operator, XMLWriter writer) {
        writer.writeOpenTag(TreeBitMoveOperatorParser.BIT_MOVE_OPERATOR,
                new Attribute[]{
                        new Attribute.Default<Double>("weight", operator.weight),
                }
        );
        writer.writeTag(TreeModel.TREE_MODEL, new Attribute[]{new Attribute.Default<String>(XMLParser.IDREF, "treeModel")}, true);
        writer.writeCloseTag(TreeBitMoveOperatorParser.BIT_MOVE_OPERATOR);
    }

    private void writeUniformOperator(Operator operator, XMLWriter writer) {
        writer.writeOpenTag("uniformOperator",
                new Attribute.Default<Double>("weight", operator.weight));
        writeParameter1Ref(writer, operator);
        writer.writeCloseTag("uniformOperator");
    }

    private void writeIntegerUniformOperator(Operator operator, XMLWriter writer) {
        writer.writeOpenTag("uniformIntegerOperator",
                new Attribute.Default<Double>("weight", operator.weight));
        writeParameter1Ref(writer, operator);
        writer.writeCloseTag("uniformIntegerOperator");
    }

    private void writeNarrowExchangeOperator(Operator operator, XMLWriter writer) {
        writer.writeOpenTag(ExchangeOperatorParser.NARROW_EXCHANGE,
                new Attribute.Default<Double>("weight", operator.weight));
        writer.writeTag(TreeModel.TREE_MODEL, new Attribute[]{new Attribute.Default<String>(XMLParser.IDREF, "treeModel")}, true);
        writer.writeCloseTag(ExchangeOperatorParser.NARROW_EXCHANGE);
    }

    private void writeWideExchangeOperator(Operator operator, XMLWriter writer) {
        writer.writeOpenTag(ExchangeOperatorParser.WIDE_EXCHANGE,
                new Attribute.Default<Double>("weight", operator.weight));
        writer.writeTag(TreeModel.TREE_MODEL, new Attribute[]{new Attribute.Default<String>(XMLParser.IDREF, "treeModel")}, true);
        writer.writeCloseTag(ExchangeOperatorParser.WIDE_EXCHANGE);
    }

    private void writeWilsonBaldingOperator(Operator operator, XMLWriter writer) {
        writer.writeOpenTag(WilsonBaldingParser.WILSON_BALDING,
                new Attribute.Default<Double>("weight", operator.weight));
        writer.writeTag(TreeModel.TREE_MODEL, new Attribute[]{new Attribute.Default<String>(XMLParser.IDREF, "treeModel")}, true);
        if (nodeHeightPrior == CONSTANT) {
            writeNodeHeightPriorModelRef(writer);
        }
        writer.writeCloseTag(WilsonBaldingParser.WILSON_BALDING);
    }

    private void writeSampleNonActiveOperator(Operator operator, XMLWriter writer) {
        writer.writeOpenTag(SampleNonActiveGibbsOperatorParser.SAMPLE_NONACTIVE_GIBBS_OPERATOR,
                new Attribute.Default<Double>("weight", operator.weight));

        writer.writeOpenTag(SampleNonActiveGibbsOperatorParser.DISTRIBUTION);
        writeParameterRefByName(writer, operator.name);
        writer.writeCloseTag(SampleNonActiveGibbsOperatorParser.DISTRIBUTION);

        writer.writeOpenTag(SampleNonActiveGibbsOperatorParser.DATA_PARAMETER);
        writeParameter1Ref(writer, operator);
        writer.writeCloseTag(SampleNonActiveGibbsOperatorParser.DATA_PARAMETER);

        writer.writeOpenTag(SampleNonActiveGibbsOperatorParser.INDICATOR_PARAMETER);
        writeParameterRefByName(writer, operator.parameter2.getName());
        writer.writeCloseTag(SampleNonActiveGibbsOperatorParser.INDICATOR_PARAMETER);

        writer.writeCloseTag(SampleNonActiveGibbsOperatorParser.SAMPLE_NONACTIVE_GIBBS_OPERATOR);
    }

    private void writeScaleWithIndicatorsOperator(Operator operator, XMLWriter writer) {
        writer.writeOpenTag(
                ScaleOperatorParser.SCALE_OPERATOR,
                new Attribute[]{
                        new Attribute.Default<Double>(ScaleOperatorParser.SCALE_FACTOR, operator.tuning),
                        new Attribute.Default<Double>("weight", operator.weight),
                });
        writeParameter1Ref(writer, operator);
        writer.writeOpenTag(ScaleOperatorParser.INDICATORS, new Attribute.Default<String>(ScaleOperatorParser.PICKONEPROB, "1.0"));
        writeParameterRefByName(writer, operator.parameter2.getName());
        writer.writeCloseTag(ScaleOperatorParser.INDICATORS);
        writer.writeCloseTag(ScaleOperatorParser.SCALE_OPERATOR);
    }

    private void writeSubtreeSlideOperator(Operator operator, XMLWriter writer) {
        writer.writeOpenTag(SubtreeSlideOperatorParser.SUBTREE_SLIDE,
                new Attribute[]{
                        new Attribute.Default<Double>("size", operator.tuning),
                        new Attribute.Default<String>("gaussian", "true"),
                        new Attribute.Default<Double>("weight", operator.weight)
                }
        );
        writer.writeTag(TreeModel.TREE_MODEL, new Attribute[]{new Attribute.Default<String>(XMLParser.IDREF, "treeModel")}, true);
        writer.writeCloseTag(SubtreeSlideOperatorParser.SUBTREE_SLIDE);
    }

    /**
     * Write the timer report block.
     *
     * @param writer the writer
     */
    public void writeTimerReport(XMLWriter writer) {
        writer.writeOpenTag("report");
        writer.writeOpenTag("property", new Attribute.Default<String>("name", "timer"));
        writer.writeTag("object", new Attribute.Default<String>(XMLParser.IDREF, "mcmc"), true);
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
                        new Attribute.Default<String>("fileName", logFileName)
                },
                true
        );
    }

    public void writeAnalysisToCSVfile(XMLWriter writer) {
        if (nodeHeightPrior == EXTENDED_SKYLINE) {
            writer.writeOpenTag(EBSPAnalysisParser.VD_ANALYSIS, new Attribute[]{
                    new Attribute.Default<String>(XMLParser.ID, "demographic.analysis"),
                    new Attribute.Default<Double>(EBSPAnalysisParser.BURN_IN, 0.1)}
            );

            writer.writeOpenTag(EBSPAnalysisParser.LOG_FILE_NAME);
            writer.writeText(logFileName);
            writer.writeCloseTag(EBSPAnalysisParser.LOG_FILE_NAME);

            writer.writeOpenTag(EBSPAnalysisParser.TREE_FILE_NAMES);
            writer.writeOpenTag(EBSPAnalysisParser.TREE_LOG);
            writer.writeText(treeFileName);
            writer.writeCloseTag(EBSPAnalysisParser.TREE_LOG);
            writer.writeCloseTag(EBSPAnalysisParser.TREE_FILE_NAMES);

            writer.writeOpenTag(EBSPAnalysisParser.MODEL_TYPE);
            writer.writeText(extendedSkylineModel);
            writer.writeCloseTag(EBSPAnalysisParser.MODEL_TYPE);

            writer.writeOpenTag(EBSPAnalysisParser.POPULATION_FIRST_COLUMN);
            writer.writeText(VariableDemographicModelParser.demoElementName + ".popSize" + 1);
            writer.writeCloseTag(EBSPAnalysisParser.POPULATION_FIRST_COLUMN);

            writer.writeOpenTag(EBSPAnalysisParser.INDICATORS_FIRST_COLUMN);
            writer.writeText(VariableDemographicModelParser.demoElementName + ".indicators" + 1);
            writer.writeCloseTag(EBSPAnalysisParser.INDICATORS_FIRST_COLUMN);

            writer.writeCloseTag(EBSPAnalysisParser.VD_ANALYSIS);

            writer.writeOpenTag(CSVExporterParser.CSV_EXPORT,
                    new Attribute[]{
                            new Attribute.Default<String>(CSVExporterParser.FILE_NAME,
                                    logFileName.subSequence(0, logFileName.length() - 4) + ".csv"),
                            new Attribute.Default<String>(CSVExporterParser.SEPARATOR, ",")
                    });
            writer.writeOpenTag(CSVExporterParser.COLUMNS);
            writer.writeTag(EBSPAnalysisParser.VD_ANALYSIS,
                    new Attribute[]{new Attribute.Default<String>(XMLParser.IDREF, "demographic.analysis")}, true);
            writer.writeCloseTag(CSVExporterParser.COLUMNS);
            writer.writeCloseTag(CSVExporterParser.CSV_EXPORT);
        }

    }

    /**
     * Write the MCMC block.
     *
     * @param writer the writer
     */
    public void writeMCMC(XMLWriter writer) {
        writer.writeOpenTag(
                "mcmc",
                new Attribute[]{
                        new Attribute.Default<String>(XMLParser.ID, "mcmc"),
                        new Attribute.Default<Integer>("chainLength", chainLength),
                        new Attribute.Default<String>("autoOptimize", autoOptimize ? "true" : "false")
                });

        if (alignment != null) {
            // we have data...
            writer.writeOpenTag(CompoundLikelihoodParser.POSTERIOR, new Attribute.Default<String>(XMLParser.ID, "posterior"));
        }

        // write prior block
        writer.writeOpenTag(CompoundLikelihoodParser.PRIOR, new Attribute.Default<String>(XMLParser.ID, "prior"));

        writeParameterPriors(writer);


        if (nodeHeightPrior == YULE || nodeHeightPrior == BIRTH_DEATH) {
            writer.writeTag(SpeciationLikelihoodParser.SPECIATION_LIKELIHOOD, new Attribute.Default<String>(XMLParser.IDREF, "speciation"), true);
        } else if (nodeHeightPrior == SKYLINE) {
            writer.writeTag(BayesianSkylineLikelihoodParser.SKYLINE_LIKELIHOOD, new Attribute.Default<String>(XMLParser.IDREF, "skyline"), true);
        } else {
            writer.writeTag(CoalescentLikelihoodParser.COALESCENT_LIKELIHOOD, new Attribute.Default<String>(XMLParser.IDREF, "coalescent"), true);
        }

        if (nodeHeightPrior == LOGISTIC) {
            writer.writeTag(BooleanLikelihoodParser.BOOLEAN_LIKELIHOOD, new Attribute.Default<String>(XMLParser.IDREF, "booleanLikelihood1"), true);
        }

        if (nodeHeightPrior == SKYLINE) {
            writer.writeTag(ExponentialMarkovModel.EXPONENTIAL_MARKOV_MODEL, new Attribute.Default<String>(XMLParser.IDREF, "eml1"), true);
        }
        if (nodeHeightPrior == EXTENDED_SKYLINE) {
            writer.writeOpenTag(MixedDistributionLikelihoodParser.DISTRIBUTION_LIKELIHOOD);

            writer.writeOpenTag(MixedDistributionLikelihoodParser.DISTRIBUTION0);
            writer.writeTag(ExponentialDistributionModel.EXPONENTIAL_DISTRIBUTION_MODEL,
                    new Attribute.Default<String>(XMLParser.IDREF, "demographic.populationMeanDist"), true);
            writer.writeCloseTag(MixedDistributionLikelihoodParser.DISTRIBUTION0);

            writer.writeOpenTag(MixedDistributionLikelihoodParser.DISTRIBUTION1);
            writer.writeTag(ExponentialDistributionModel.EXPONENTIAL_DISTRIBUTION_MODEL,
                    new Attribute.Default<String>(XMLParser.IDREF, "demographic.populationMeanDist"), true);
            writer.writeCloseTag(MixedDistributionLikelihoodParser.DISTRIBUTION1);

            writer.writeOpenTag(MixedDistributionLikelihoodParser.DATA);
            writer.writeTag(ParameterParser.PARAMETER,
                    new Attribute.Default<String>(XMLParser.IDREF, "demographic.popSize"), true);
            writer.writeCloseTag(MixedDistributionLikelihoodParser.DATA);

            writer.writeOpenTag(MixedDistributionLikelihoodParser.INDICATORS);
            writer.writeTag(ParameterParser.PARAMETER,
                    new Attribute.Default<String>(XMLParser.IDREF, "demographic.indicators"), true);
            writer.writeCloseTag(MixedDistributionLikelihoodParser.INDICATORS);

            writer.writeCloseTag(MixedDistributionLikelihoodParser.DISTRIBUTION_LIKELIHOOD);
        }
        writer.writeCloseTag(CompoundLikelihoodParser.PRIOR);

        if (alignment != null) {
            // write likelihood block
            writer.writeOpenTag(CompoundLikelihoodParser.LIKELIHOOD, new Attribute.Default<String>(XMLParser.ID, "likelihood"));

            boolean nucs = alignment.getDataType() == Nucleotides.INSTANCE;
            if (nucs && codonHeteroPattern != null) {
                for (int i = 1; i <= partitionCount; i++) {
                    writer.writeTag(TreeLikelihoodParser.TREE_LIKELIHOOD, new Attribute.Default<String>(XMLParser.IDREF, TreeLikelihoodParser.TREE_LIKELIHOOD + i), true);
                }
            } else {
                writer.writeTag(TreeLikelihoodParser.TREE_LIKELIHOOD, new Attribute.Default<String>(XMLParser.IDREF, TreeLikelihoodParser.TREE_LIKELIHOOD), true);
            }

            writer.writeCloseTag(CompoundLikelihoodParser.LIKELIHOOD);


            writer.writeCloseTag(CompoundLikelihoodParser.POSTERIOR);
        }

        writer.writeTag(SimpleOperatorScheduleParser.OPERATOR_SCHEDULE, new Attribute.Default<String>(XMLParser.IDREF, "operators"), true);

        // write log to screen
        writer.writeOpenTag(LoggerParser.LOG,
                new Attribute[]{
                        new Attribute.Default<String>(XMLParser.ID, "screenLog"),
                        new Attribute.Default<String>(LoggerParser.LOG_EVERY, echoEvery + "")
                });
        writeScreenLog(writer);
        writer.writeCloseTag(LoggerParser.LOG);

        // write log to file
        if (logFileName == null) {
            logFileName = fileNameStem + ".log";
        }
        writer.writeOpenTag(LoggerParser.LOG,
                new Attribute[]{
                        new Attribute.Default<String>(XMLParser.ID, "fileLog"),
                        new Attribute.Default<String>(LoggerParser.LOG_EVERY, logEvery + ""),
                        new Attribute.Default<String>(LoggerParser.FILE_NAME, logFileName)
                });
        writeLog(writer);
        writer.writeCloseTag(LoggerParser.LOG);

        // write tree log to file
        if (treeFileName == null) {
            if (substTreeLog) {
                treeFileName = fileNameStem + "(time).trees";
            } else {
                treeFileName = fileNameStem + ".trees";
            }
        }
        writer.writeOpenTag(TreeLoggerParser.LOG_TREE,
                new Attribute[]{
                        new Attribute.Default<String>(XMLParser.ID, "treeFileLog"),
                        new Attribute.Default<String>(TreeLoggerParser.LOG_EVERY, logEvery + ""),
                        new Attribute.Default<String>(TreeLoggerParser.NEXUS_FORMAT, "true"),
                        new Attribute.Default<String>(TreeLoggerParser.FILE_NAME, treeFileName),
                        new Attribute.Default<String>(TreeLoggerParser.SORT_TRANSLATION_TABLE, "true")
                });
        writer.writeTag(TreeModel.TREE_MODEL, new Attribute.Default<String>(XMLParser.IDREF, "treeModel"), true);
        if (clockModel != STRICT_CLOCK) {
            writer.writeTag(DiscretizedBranchRatesParser.DISCRETIZED_BRANCH_RATES, new Attribute[]{new Attribute.Default<String>(XMLParser.IDREF, BranchRateModel.BRANCH_RATES)}, true);
        }
        if (alignment != null) {
            // we have data...
            writer.writeTag("posterior", new Attribute.Default<String>(XMLParser.IDREF, "posterior"), true);
        }
        writer.writeCloseTag(TreeLoggerParser.LOG_TREE);

//        if (mapTreeLog) {
//            // write tree log to file
//            if (mapTreeFileName == null) {
//                mapTreeFileName = fileNameStem + ".MAP.tree";
//            }
//            writer.writeOpenTag("logML",
//                    new Attribute[] {
//                        new Attribute.Default<String>(TreeLogger.FILE_NAME, mapTreeFileName)
//                    });
//            writer.writeOpenTag("ml");
//            writer.writeTag(CompoundLikelihood.POSTERIOR, new Attribute.Default<String>(XMLParser.IDREF, "posterior"), true);
//            writer.writeCloseTag("ml");
//            writer.writeOpenTag("column", new Attribute[] {
//                        new Attribute.Default<String>("label", "MAP tree")
//                    });
//            writer.writeTag(TreeModel.TREE_MODEL, new Attribute.Default<String>(XMLParser.IDREF, "treeModel"), true);
//            writer.writeCloseTag("column");
//            writer.writeCloseTag("logML");
//        }

        if (substTreeLog) {
            // write tree log to file
            if (substTreeFileName == null) {
                substTreeFileName = fileNameStem + "(subst).trees";
            }
            writer.writeOpenTag(TreeLoggerParser.LOG_TREE,
                    new Attribute[]{
                            new Attribute.Default<String>(XMLParser.ID, "substTreeFileLog"),
                            new Attribute.Default<String>(TreeLoggerParser.LOG_EVERY, logEvery + ""),
                            new Attribute.Default<String>(TreeLoggerParser.NEXUS_FORMAT, "true"),
                            new Attribute.Default<String>(TreeLoggerParser.FILE_NAME, substTreeFileName),
                            new Attribute.Default<String>(TreeLoggerParser.BRANCH_LENGTHS, TreeLoggerParser.SUBSTITUTIONS)
                    });
            writer.writeTag(TreeModel.TREE_MODEL, new Attribute.Default<String>(XMLParser.IDREF, "treeModel"), true);
            if (clockModel == STRICT_CLOCK) {
                writer.writeTag(StrictClockBranchRatesParser.STRICT_CLOCK_BRANCH_RATES, new Attribute[]{new Attribute.Default<String>(XMLParser.IDREF, BranchRateModel.BRANCH_RATES)}, true);
            } else {
                writer.writeTag(DiscretizedBranchRatesParser.DISCRETIZED_BRANCH_RATES, new Attribute[]{new Attribute.Default<String>(XMLParser.IDREF, BranchRateModel.BRANCH_RATES)}, true);
            }
            writer.writeCloseTag(TreeLoggerParser.LOG_TREE);
        }

        writer.writeCloseTag("mcmc");
    }

    /**
     * Write the priors for each parameter
     *
     * @param writer the writer
     */
    private void writeParameterPriors(XMLWriter writer) {
        boolean first = true;

        for (Map.Entry<Taxa, Boolean> taxa : taxonSetsMono.entrySet()) {
            if ( taxa.getValue() ) {
                if (first) {
                    writer.writeOpenTag(BooleanLikelihoodParser.BOOLEAN_LIKELIHOOD);
                    first = false;
                }
                final String taxaRef = "monophyly(" + taxa.getKey().getId() + ")";
                final Attribute.Default attr = new Attribute.Default<String>(XMLParser.IDREF, taxaRef);
                writer.writeTag(MonophylyStatisticParser.MONOPHYLY_STATISTIC, new Attribute[]{attr}, true);
            }
        }
        if (!first) {
            writer.writeCloseTag(BooleanLikelihoodParser.BOOLEAN_LIKELIHOOD);
        }

        ArrayList<Parameter> parameters = selectParameters();
        for (Parameter parameter : parameters) {
            if (parameter.priorType != PriorType.NONE) {
                if (parameter.priorType != PriorType.UNIFORM_PRIOR || parameter.isNodeHeight) {
                    writeParameterPrior(parameter, writer);
                }
            }
        }

    }

    /**
     * Write the priors for each parameter
     *
     * @param parameter the parameter
     * @param writer    the writer
     */
    private void writeParameterPrior(Parameter parameter, XMLWriter writer) {
        switch (parameter.priorType) {
            case UNIFORM_PRIOR:
                writer.writeOpenTag(PriorParsers.UNIFORM_PRIOR,
                        new Attribute[]{
                                new Attribute.Default<String>(PriorParsers.LOWER, "" + parameter.uniformLower),
                                new Attribute.Default<String>(PriorParsers.UPPER, "" + parameter.uniformUpper)
                        });
                writeParameterIdref(writer, parameter);
                writer.writeCloseTag(PriorParsers.UNIFORM_PRIOR);
                break;
            case EXPONENTIAL_PRIOR:
                writer.writeOpenTag(PriorParsers.EXPONENTIAL_PRIOR,
                        new Attribute[]{
                                new Attribute.Default<String>(PriorParsers.MEAN, "" + parameter.exponentialMean),
                                new Attribute.Default<String>(PriorParsers.OFFSET, "" + parameter.exponentialOffset)
                        });
                writeParameterIdref(writer, parameter);
                writer.writeCloseTag(PriorParsers.EXPONENTIAL_PRIOR);
                break;
            case NORMAL_PRIOR:
                writer.writeOpenTag(PriorParsers.NORMAL_PRIOR,
                        new Attribute[]{
                                new Attribute.Default<String>(PriorParsers.MEAN, "" + parameter.normalMean),
                                new Attribute.Default<String>(PriorParsers.STDEV, "" + parameter.normalStdev)
                        });
                writeParameterIdref(writer, parameter);
                writer.writeCloseTag(PriorParsers.NORMAL_PRIOR);
                break;
            case LOGNORMAL_PRIOR:
                writer.writeOpenTag(PriorParsers.LOG_NORMAL_PRIOR,
                        new Attribute[]{
                                new Attribute.Default<String>(PriorParsers.MEAN, "" + parameter.logNormalMean),
                                new Attribute.Default<String>(PriorParsers.STDEV, "" + parameter.logNormalStdev),
                                new Attribute.Default<String>(PriorParsers.OFFSET, "" + parameter.logNormalOffset),

                                // this is to be implemented...
                                new Attribute.Default<String>(PriorParsers.MEAN_IN_REAL_SPACE, "false")
                        });
                writeParameterIdref(writer, parameter);
                writer.writeCloseTag(PriorParsers.LOG_NORMAL_PRIOR);
                break;
            case GAMMA_PRIOR:
                writer.writeOpenTag(PriorParsers.GAMMA_PRIOR,
                        new Attribute[]{
                                new Attribute.Default<String>(PriorParsers.SHAPE, "" + parameter.gammaAlpha),
                                new Attribute.Default<String>(PriorParsers.SCALE, "" + parameter.gammaBeta),
                                new Attribute.Default<String>(PriorParsers.OFFSET, "" + parameter.gammaOffset)
                        });
                writeParameterIdref(writer, parameter);
                writer.writeCloseTag(PriorParsers.GAMMA_PRIOR);
                break;
            case JEFFREYS_PRIOR:
                writer.writeOpenTag(OneOnXPriorParser.ONE_ONE_X_PRIOR);
                writeParameterIdref(writer, parameter);
                writer.writeCloseTag(OneOnXPriorParser.ONE_ONE_X_PRIOR);
                break;
            case POISSON_PRIOR:
                writer.writeOpenTag(PriorParsers.POISSON_PRIOR,
                        new Attribute[]{
                                new Attribute.Default<String>(PriorParsers.MEAN, "" + parameter.poissonMean),
                                new Attribute.Default<String>(PriorParsers.OFFSET, "" + parameter.poissonOffset)
                        });
                writeParameterIdref(writer, parameter);
                writer.writeCloseTag(PriorParsers.POISSON_PRIOR);
                break;
            case TRUNC_NORMAL_PRIOR:
                writer.writeOpenTag(PriorParsers.UNIFORM_PRIOR,
                        new Attribute[]{
                                new Attribute.Default<String>(PriorParsers.LOWER, "" + parameter.uniformLower),
                                new Attribute.Default<String>(PriorParsers.UPPER, "" + parameter.uniformUpper)
                        });
                writeParameterIdref(writer, parameter);
                writer.writeCloseTag(PriorParsers.UNIFORM_PRIOR);
                writer.writeOpenTag(PriorParsers.NORMAL_PRIOR,
                        new Attribute[]{
                                new Attribute.Default<String>(PriorParsers.MEAN, "" + parameter.normalMean),
                                new Attribute.Default<String>(PriorParsers.STDEV, "" + parameter.normalStdev)
                        });
                writeParameterIdref(writer, parameter);
                writer.writeCloseTag(PriorParsers.NORMAL_PRIOR);
                break;
            default:
                throw new IllegalArgumentException("Unknown priorType");
        }
    }

    private void writeParameterIdref(XMLWriter writer, Parameter parameter) {
        if (parameter.isStatistic) {
            writer.writeTag("statistic", new Attribute[]{new Attribute.Default<String>(XMLParser.IDREF, parameter.getName())}, true);
        } else {
            writer.writeTag(ParameterParser.PARAMETER, new Attribute[]{new Attribute.Default<String>(XMLParser.IDREF, parameter.getName())}, true);
        }
    }

    private void writeSumStatisticColumn(XMLWriter writer, String name, String label) {
        writer.writeOpenTag(ColumnsParser.COLUMN,
                new Attribute[]{
                        new Attribute.Default<String>(ColumnsParser.LABEL, label),
                        new Attribute.Default<String>(ColumnsParser.DECIMAL_PLACES, "0"),
                        new Attribute.Default<String>(ColumnsParser.WIDTH, "12")
                }
        );
        writer.writeTag(SumStatisticParser.SUM_STATISTIC, new Attribute.Default<String>(XMLParser.IDREF, name), true);
        writer.writeCloseTag(ColumnsParser.COLUMN);
    }

    /**
     * Write the log
     *
     * @param writer the writer
     */
    private void writeScreenLog(XMLWriter writer) {
        if (alignment != null) {
            writer.writeOpenTag(ColumnsParser.COLUMN,
                    new Attribute[]{
                            new Attribute.Default<String>(ColumnsParser.LABEL, "Posterior"),
                            new Attribute.Default<String>(ColumnsParser.DECIMAL_PLACES, "4"),
                            new Attribute.Default<String>(ColumnsParser.WIDTH, "12")
                    }
            );
            writer.writeTag(CompoundLikelihoodParser.POSTERIOR, new Attribute.Default<String>(XMLParser.IDREF, "posterior"), true);
            writer.writeCloseTag(ColumnsParser.COLUMN);
        }

        writer.writeOpenTag(ColumnsParser.COLUMN,
                new Attribute[]{
                        new Attribute.Default<String>(ColumnsParser.LABEL, "Prior"),
                        new Attribute.Default<String>(ColumnsParser.DECIMAL_PLACES, "4"),
                        new Attribute.Default<String>(ColumnsParser.WIDTH, "12")
                }
        );
        writer.writeTag(CompoundLikelihoodParser.PRIOR, new Attribute.Default<String>(XMLParser.IDREF, "prior"), true);
        writer.writeCloseTag(ColumnsParser.COLUMN);

        if (alignment != null) {
            writer.writeOpenTag(ColumnsParser.COLUMN,
                    new Attribute[]{
                            new Attribute.Default<String>(ColumnsParser.LABEL, "Likelihood"),
                            new Attribute.Default<String>(ColumnsParser.DECIMAL_PLACES, "4"),
                            new Attribute.Default<String>(ColumnsParser.WIDTH, "12")
                    }
            );
            writer.writeTag(CompoundLikelihoodParser.LIKELIHOOD, new Attribute.Default<String>(XMLParser.IDREF, "likelihood"), true);
            writer.writeCloseTag(ColumnsParser.COLUMN);
        }

        writer.writeOpenTag(ColumnsParser.COLUMN,
                new Attribute[]{
                        new Attribute.Default<String>(ColumnsParser.LABEL, "Root Height"),
                        new Attribute.Default<String>(ColumnsParser.SIGNIFICANT_FIGURES, "6"),
                        new Attribute.Default<String>(ColumnsParser.WIDTH, "12")
                }
        );
        writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default<String>(XMLParser.IDREF, "treeModel.rootHeight"), true);
        writer.writeCloseTag(ColumnsParser.COLUMN);

        writer.writeOpenTag(ColumnsParser.COLUMN,
                new Attribute[]{
                        new Attribute.Default<String>(ColumnsParser.LABEL, "Rate"),
                        new Attribute.Default<String>(ColumnsParser.SIGNIFICANT_FIGURES, "6"),
                        new Attribute.Default<String>(ColumnsParser.WIDTH, "12")
                }
        );
        if (clockModel == STRICT_CLOCK) {
            writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default<String>(XMLParser.IDREF, "clock.rate"), true);
        } else {
            writer.writeTag(RateStatisticParser.RATE_STATISTIC, new Attribute.Default<String>(XMLParser.IDREF, "meanRate"), true);
        }
        writer.writeCloseTag(ColumnsParser.COLUMN);

        if (clockModel == RANDOM_LOCAL_CLOCK) {
            writeSumStatisticColumn(writer, "rateChanges", "Rate Changes");
        }

        // I think this is too much info for the screen - it is all in the log file.
//		if (alignment != null) {
//			boolean nucs = alignment.getDataType() == Nucleotides.INSTANCE;
//			if (nucs && codonHeteroPattern != null) {
//				if (codonHeteroPattern.equals("112")) {
//					writer.writeOpenTag(ColumnsParser.COLUMN,
//							new Attribute[] {
//									new Attribute.Default<String>(ColumnsParser.LABEL, "L(codon pos 1+2)"),
//									new Attribute.Default<String>(ColumnsParser.DECIMAL_PLACES, "4"),
//									new Attribute.Default<String>(ColumnsParser.WIDTH, "12")
//							}
//					);
//					writer.writeTag(OldTreeLikelihoodParser.TREE_LIKELIHOOD, new Attribute.Default<String>(XMLParser.IDREF,"treeLikelihood1"), true);
//					writer.writeCloseTag(ColumnsParser.COLUMN);
//					writer.writeOpenTag(ColumnsParser.COLUMN,
//							new Attribute[] {
//									new Attribute.Default<String>(ColumnsParser.LABEL, "L(codon pos 3)"),
//									new Attribute.Default<String>(ColumnsParser.DECIMAL_PLACES, "4"),
//									new Attribute.Default<String>(ColumnsParser.WIDTH, "12")
//							}
//					);
//					writer.writeTag(OldTreeLikelihoodParser.TREE_LIKELIHOOD, new Attribute.Default<String>(XMLParser.IDREF,"treeLikelihood2"), true);
//					writer.writeCloseTag(ColumnsParser.COLUMN);
//				} else if (codonHeteroPattern.equals("123")) {
//					for (int i =1; i <= 3; i++) {
//						writer.writeOpenTag(ColumnsParser.COLUMN,
//								new Attribute[] {
//										new Attribute.Default<String>(ColumnsParser.LABEL, "L(codon pos "+i+")"),
//										new Attribute.Default<String>(ColumnsParser.DECIMAL_PLACES, "4"),
//										new Attribute.Default<String>(ColumnsParser.WIDTH, "12")
//								}
//						);
//						writer.writeTag(OldTreeLikelihoodParser.TREE_LIKELIHOOD, new Attribute.Default<String>(XMLParser.IDREF,OldTreeLikelihoodParser.TREE_LIKELIHOOD + i), true);
//						writer.writeCloseTag(ColumnsParser.COLUMN);
//					}
//				}
//			} else {
//				writer.writeOpenTag(ColumnsParser.COLUMN,
//						new Attribute[] {
//								new Attribute.Default<String>(ColumnsParser.LABEL, "L(tree)"),
//								new Attribute.Default<String>(ColumnsParser.DECIMAL_PLACES, "4"),
//								new Attribute.Default<String>(ColumnsParser.WIDTH, "12")
//						}
//				);
//				writer.writeTag(OldTreeLikelihoodParser.TREE_LIKELIHOOD, new Attribute.Default<String>(XMLParser.IDREF,OldTreeLikelihoodParser.TREE_LIKELIHOOD), true);
//				writer.writeCloseTag(ColumnsParser.COLUMN);
//			}
//		}
//		if (nodeHeightPrior == YULE || nodeHeightPrior == BIRTH_DEATH) {
//			writer.writeOpenTag(ColumnsParser.COLUMN,
//					new Attribute[] {
//							new Attribute.Default<String>(ColumnsParser.LABEL, "L(speciation)"),
//							new Attribute.Default<String>(ColumnsParser.DECIMAL_PLACES, "4"),
//							new Attribute.Default<String>(ColumnsParser.WIDTH, "12")
//					}
//			);
//			writer.writeTag(SpeciationLikelihoodParser.SPECIATION_LIKELIHOOD, new Attribute.Default<String>(XMLParser.IDREF, "speciation"), true);
//		} else {
//			writer.writeOpenTag(ColumnsParser.COLUMN,
//					new Attribute[] {
//							new Attribute.Default<String>(ColumnsParser.LABEL, "L(coalecent)"),
//							new Attribute.Default<String>(ColumnsParser.DECIMAL_PLACES, "4"),
//							new Attribute.Default<String>(ColumnsParser.WIDTH, "12")
//					}
//			);
//			if (nodeHeightPrior == SKYLINE) {
//				writer.writeTag(BayesianSkylineLikelihoodParser.SKYLINE_LIKELIHOOD, new Attribute.Default<String>(XMLParser.IDREF, "skyline"), true);
//			} else {
//				writer.writeTag(CoalescentLikelihoodParser.COALESCENT_LIKELIHOOD, new Attribute.Default<String>(XMLParser.IDREF,"coalescent"), true);
//			}
//		}
//		writer.writeCloseTag(ColumnsParser.COLUMN);

    }

    /**
     * Write the log
     *
     * @param writer the writer
     */
    private void writeLog(XMLWriter writer) {
        if (alignment != null) {
            writer.writeTag(CompoundLikelihoodParser.POSTERIOR, new Attribute.Default<String>(XMLParser.IDREF, "posterior"), true);
        }
        writer.writeTag(CompoundLikelihoodParser.PRIOR, new Attribute.Default<String>(XMLParser.IDREF, "prior"), true);
        if (alignment != null) {
            writer.writeTag(CompoundLikelihoodParser.LIKELIHOOD, new Attribute.Default<String>(XMLParser.IDREF, "likelihood"), true);
        }

        // As of v1.4.2, always write the rate parameter even if fixed...
        //if (!fixedSubstitutionRate) {
        if (clockModel == STRICT_CLOCK) {
            writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default<String>(XMLParser.IDREF, "clock.rate"), true);
        } else {
            writer.writeTag(RateStatisticParser.RATE_STATISTIC, new Attribute.Default<String>(XMLParser.IDREF, "meanRate"), true);
        }
        //}

        writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default<String>(XMLParser.IDREF, "treeModel.rootHeight"), true);

        for (Taxa taxa : taxonSets) {
            writer.writeTag("tmrcaStatistic", new Attribute[]{new Attribute.Default<String>(XMLParser.IDREF, "tmrca(" + taxa.getId() + ")")}, true);
        }

        if (nodeHeightPrior == CONSTANT) {
            writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default<String>(XMLParser.IDREF, "constant.popSize"), true);
        } else if (nodeHeightPrior == EXPONENTIAL) {
            writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default<String>(XMLParser.IDREF, "exponential.popSize"), true);
            if (parameterization == GROWTH_RATE) {
                writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default<String>(XMLParser.IDREF, "exponential.growthRate"), true);
            } else {
                writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default<String>(XMLParser.IDREF, "exponential.doublingTime"), true);
            }
        } else if (nodeHeightPrior == LOGISTIC) {
            writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default<String>(XMLParser.IDREF, "logistic.popSize"), true);
            if (parameterization == GROWTH_RATE) {
                writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default<String>(XMLParser.IDREF, "logistic.growthRate"), true);
            } else {
                writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default<String>(XMLParser.IDREF, "logistic.doublingTime"), true);
            }
            writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default<String>(XMLParser.IDREF, "logistic.t50"), true);
        } else if (nodeHeightPrior == EXPANSION) {
            writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default<String>(XMLParser.IDREF, "expansion.popSize"), true);
            if (parameterization == GROWTH_RATE) {
                writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default<String>(XMLParser.IDREF, "expansion.growthRate"), true);
            } else {
                writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default<String>(XMLParser.IDREF, "expansion.doublingTime"), true);
            }
            writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default<String>(XMLParser.IDREF, "expansion.ancestralProportion"), true);
        } else if (nodeHeightPrior == SKYLINE) {
            writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default<String>(XMLParser.IDREF, "skyline.popSize"), true);
            writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default<String>(XMLParser.IDREF, "skyline.groupSize"), true);
        } else if (nodeHeightPrior == EXTENDED_SKYLINE) {
            writeSumStatisticColumn(writer, "demographic.populationSizeChanges", "popSize_changes");
            writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default<String>(XMLParser.IDREF, "demographic.populationMean"), true);
            writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default<String>(XMLParser.IDREF, "demographic.popSize"), true);
            writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default<String>(XMLParser.IDREF, "demographic.indicators"), true);
        } else if (nodeHeightPrior == YULE) {
            writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default<String>(XMLParser.IDREF, "yule.birthRate"), true);
        } else if (nodeHeightPrior == BIRTH_DEATH) {
            writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default<String>(XMLParser.IDREF, BirthDeathModelParser.MEAN_GROWTH_RATE_PARAM_NAME), true);
            writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default<String>(XMLParser.IDREF, BirthDeathModelParser.RELATIVE_DEATH_RATE_PARAM_NAME), true);
        }

        if (alignment != null) {
            switch (dataType) {
                case DataType.NUCLEOTIDES:
                    if (partitionCount > 1) {
                        for (int i = 1; i <= partitionCount; i++) {
                            writer.writeTag(ParameterParser.PARAMETER,
                                    new Attribute.Default<String>(XMLParser.IDREF, SiteModel.SITE_MODEL + i + ".mu"), true);
                        }
                    }
                    switch (nucSubstitutionModel) {
                        case HKY:
                            if (partitionCount > 1 && unlinkedSubstitutionModel) {
                                for (int i = 1; i <= partitionCount; i++) {
                                    writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default<String>(XMLParser.IDREF, "hky" + i + ".kappa"), true);
                                }
                            } else {
                                writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default<String>(XMLParser.IDREF, "hky.kappa"), true);
                            }
                            break;

                        case GTR:
                            if (partitionCount > 1 && unlinkedSubstitutionModel) {
                                for (int i = 1; i <= partitionCount; i++) {
                                    writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default<String>(XMLParser.IDREF, "gtr" + i + ".ac"), true);
                                    writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default<String>(XMLParser.IDREF, "gtr" + i + ".ag"), true);
                                    writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default<String>(XMLParser.IDREF, "gtr" + i + ".at"), true);
                                    writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default<String>(XMLParser.IDREF, "gtr" + i + ".cg"), true);
                                    writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default<String>(XMLParser.IDREF, "gtr" + i + ".gt"), true);
                                }
                            } else {
                                writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default<String>(XMLParser.IDREF, "gtr.ac"), true);
                                writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default<String>(XMLParser.IDREF, "gtr.ag"), true);
                                writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default<String>(XMLParser.IDREF, "gtr.at"), true);
                                writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default<String>(XMLParser.IDREF, "gtr.cg"), true);
                                writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default<String>(XMLParser.IDREF, "gtr.gt"), true);
                            }
                            break;
                    }
                    break;//NUCLEOTIDES

                case DataType.AMINO_ACIDS:
                    break;//AMINO_ACIDS

                case DataType.TWO_STATES:
                case DataType.COVARION:

                    switch (binarySubstitutionModel) {
                        case BIN_SIMPLE:
                            break;
                        case BIN_COVARION:
                            writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default<String>(XMLParser.IDREF, "bcov.alpha"), true);
                            writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default<String>(XMLParser.IDREF, "bcov.s"), true);
                            writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default<String>(XMLParser.IDREF, "bcov.frequencies"), true);
                            writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default<String>(XMLParser.IDREF, "bcov.hfrequencies"), true);
                            break;

                    }
                    break;//BINARY
            }

            if (gammaHetero) {
                if (partitionCount > 1 && unlinkedHeterogeneityModel) {
                    for (int i = 1; i <= partitionCount; i++) {
                        writer.writeTag(ParameterParser.PARAMETER,
                                new Attribute.Default<String>(XMLParser.IDREF, SiteModel.SITE_MODEL + i + ".alpha"), true);
                    }
                } else {
                    writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default<String>(XMLParser.IDREF, SiteModel.SITE_MODEL + "." + "alpha"), true);
                }
            }

            if (invarHetero) {
                if (partitionCount > 1 && unlinkedHeterogeneityModel) {
                    for (int i = 1; i <= partitionCount; i++) {
                        writer.writeTag(ParameterParser.PARAMETER,
                                new Attribute.Default<String>(XMLParser.IDREF, SiteModel.SITE_MODEL + i + ".pInv"), true);
                    }
                } else {
                    writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default<String>(XMLParser.IDREF, SiteModel.SITE_MODEL + "." + "pInv"), true);
                }
            }
        }

        if (clockModel != STRICT_CLOCK) {
//			if (!fixedSubstitutionRate) {
            if (clockModel == UNCORRELATED_EXPONENTIAL) {
                writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default<String>(XMLParser.IDREF, UCED_MEAN), true);
            } else if (clockModel == UNCORRELATED_LOGNORMAL) {
                writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default<String>(XMLParser.IDREF, UCLD_MEAN), true);
                writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default<String>(XMLParser.IDREF, UCLD_STDEV), true);
            }
//			}
            writer.writeTag(RateStatisticParser.RATE_STATISTIC, new Attribute.Default<String>(XMLParser.IDREF, RateStatisticParser.COEFFICIENT_OF_VARIATION), true);
            writer.writeTag(RateCovarianceStatisticParser.RATE_COVARIANCE_STATISTIC, new Attribute.Default<String>(XMLParser.IDREF, "covariance"), true);

            if (clockModel == RANDOM_LOCAL_CLOCK) {
                writer.writeTag(SumStatisticParser.SUM_STATISTIC, new Attribute.Default<String>(XMLParser.IDREF, "rateChanges"), true);
            }
        }

        if (alignment != null) {
            boolean nucs = alignment.getDataType() == Nucleotides.INSTANCE;
            if (nucs && partitionCount > 1) {
                for (int i = 1; i <= partitionCount; i++) {
                    writer.writeTag(TreeLikelihoodParser.TREE_LIKELIHOOD, new Attribute.Default<String>(XMLParser.IDREF, TreeLikelihoodParser.TREE_LIKELIHOOD + i), true);
                }
            } else
                writer.writeTag(TreeLikelihoodParser.TREE_LIKELIHOOD, new Attribute.Default<String>(XMLParser.IDREF, TreeLikelihoodParser.TREE_LIKELIHOOD), true);
        }
        if (nodeHeightPrior == YULE || nodeHeightPrior == BIRTH_DEATH) {
            writer.writeTag(SpeciationLikelihoodParser.SPECIATION_LIKELIHOOD, new Attribute.Default<String>(XMLParser.IDREF, "speciation"), true);
        } else if (nodeHeightPrior == SKYLINE) {
            writer.writeTag(BayesianSkylineLikelihoodParser.SKYLINE_LIKELIHOOD, new Attribute.Default<String>(XMLParser.IDREF, "skyline"), true);
        } else {
            writer.writeTag(CoalescentLikelihoodParser.COALESCENT_LIKELIHOOD, new Attribute.Default<String>(XMLParser.IDREF, "coalescent"), true);
        }


    }

    /**
     * fix a parameter
     *
     * @param id    the id
     * @param value the value
     */
    public void fixParameter(String id, double value) {
        Parameter parameter = parameters.get(id);
        if (parameter == null) {
            throw new IllegalArgumentException("parameter with name, " + id + ", is unknown");
        }
        parameter.isFixed = true;
        parameter.initial = value;
    }


    private String multiDimensionValue(int dimension, double value) {
        String multi = "";

        multi += value + "";
        for (int i = 2; i <= dimension; i++)
            multi += " " + value;

        return multi;
    }

    /**
     * write a parameter
     *
     * @param id     the id
     * @param writer the writer
     */
    public void writeParameter(String id, XMLWriter writer) {
        Parameter parameter = parameters.get(id);
        if (parameter == null) {
            throw new IllegalArgumentException("parameter with name, " + id + ", is unknown");
        }
        if (parameter.isFixed) {
            writeParameter(id, 1, parameter.initial, Double.NaN, Double.NaN, writer);
        } else {
            if (parameter.priorType == PriorType.UNIFORM_PRIOR || parameter.priorType == PriorType.TRUNC_NORMAL_PRIOR) {
                writeParameter(id, 1, parameter.initial, parameter.uniformLower, parameter.uniformUpper, writer);
            } else {
                writeParameter(id, 1, parameter.initial, parameter.lower, parameter.upper, writer);
            }
        }
    }

    /**
     * write a parameter
     *
     * @param id        the id
     * @param dimension the dimension
     * @param writer    the writer
     */
    public void writeParameter(String id, int dimension, XMLWriter writer) {
        Parameter parameter = parameters.get(id);
        if (parameter == null) {
            throw new IllegalArgumentException("parameter with name, " + id + ", is unknown");
        }
        if (parameter.isFixed) {
            writeParameter(id, dimension, parameter.initial, Double.NaN, Double.NaN, writer);
        } else if (parameter.priorType == PriorType.UNIFORM_PRIOR || parameter.priorType == PriorType.TRUNC_NORMAL_PRIOR) {
            writeParameter(id, dimension, parameter.initial, parameter.uniformLower, parameter.uniformUpper, writer);
        } else {
            writeParameter(id, dimension, parameter.initial, parameter.lower, parameter.upper, writer);
        }
    }

    /**
     * write a parameter
     *
     * @param id        the id
     * @param dimension the dimension
     * @param value     the value
     * @param lower     the lower bound
     * @param upper     the upper bound
     * @param writer    the writer
     */
    public void writeParameter(String id, int dimension, double value, double lower, double upper, XMLWriter writer) {
        ArrayList<Attribute.Default> attributes = new ArrayList<Attribute.Default>();
        attributes.add(new Attribute.Default<String>(XMLParser.ID, id));
        if (dimension > 1) {
            attributes.add(new Attribute.Default<String>("dimension", dimension + ""));
        }
        if (!Double.isNaN(value)) {
            attributes.add(new Attribute.Default<String>("value", multiDimensionValue(dimension, value)));
        }
        if (!Double.isNaN(lower)) {
            attributes.add(new Attribute.Default<String>("lower", multiDimensionValue(dimension, lower)));
        }
        if (!Double.isNaN(upper)) {
            attributes.add(new Attribute.Default<String>("upper", multiDimensionValue(dimension, upper)));
        }

        Attribute[] attrArray = new Attribute[attributes.size()];
        for (int i = 0; i < attrArray.length; i++) {
            attrArray[i] = attributes.get(i);
        }

        writer.writeTag(ParameterParser.PARAMETER, attrArray, true);
    }

    /**
     * Generate XML for the starting tree
     *
     * @param writer the writer
     */
    public void writeStartingTree(XMLWriter writer) {
        if (userTree) {
            writeUserTree(tree, writer);
        } else if (upgmaStartingTree) {
            // generate a upgma starting tree
            writer.writeComment("Construct a rough-and-ready UPGMA tree as an starting tree");
            Parameter rootHeight = getParameter("treeModel.rootHeight");
            if (rootHeight.priorType != PriorType.NONE) {
                writer.writeOpenTag(
                        UPGMATreeParser.UPGMA_TREE,
                        new Attribute[]{
                                new Attribute.Default<String>(XMLParser.ID, InitialTreeGenerator.STARTING_TREE),
                                new Attribute.Default<String>(UPGMATreeParser.ROOT_HEIGHT, "" + rootHeight.initial)
                        }
                );
            } else {
                writer.writeOpenTag(
                        UPGMATreeParser.UPGMA_TREE,
                        new Attribute[]{
                                new Attribute.Default<String>(XMLParser.ID, InitialTreeGenerator.STARTING_TREE)
                        }
                );
            }
            writer.writeOpenTag(
                    DistanceMatrixParser.DISTANCE_MATRIX,
                    new Attribute[]{
                            new Attribute.Default<String>(DistanceMatrixParser.CORRECTION, "JC")
                    }
            );
            writer.writeOpenTag(SitePatternsParser.PATTERNS);
            writer.writeTag(AlignmentParser.ALIGNMENT, new Attribute[]{new Attribute.Default<String>(XMLParser.IDREF, AlignmentParser.ALIGNMENT)}, true);
            writer.writeCloseTag(SitePatternsParser.PATTERNS);
            writer.writeCloseTag(DistanceMatrixParser.DISTANCE_MATRIX);
            writer.writeCloseTag(UPGMATreeParser.UPGMA_TREE);
        } else {
            // generate a coalescent tree
            writer.writeComment("Generate a random starting tree under the coalescent process");
            Parameter rootHeight = getParameter("treeModel.rootHeight");
            if (rootHeight.priorType != PriorType.NONE) {
                writer.writeOpenTag(
                        CoalescentSimulatorParser.COALESCENT_TREE,
                        new Attribute[]{
                                new Attribute.Default<String>(XMLParser.ID, InitialTreeGenerator.STARTING_TREE),
                                new Attribute.Default<String>(TreeModelParser.ROOT_HEIGHT, "" + rootHeight.initial)
                        }
                );
            } else {
                writer.writeOpenTag(
                        CoalescentSimulatorParser.COALESCENT_TREE,
                        new Attribute[]{
                                new Attribute.Default<String>(XMLParser.ID, InitialTreeGenerator.STARTING_TREE)
                        }
                );
            }

            Attribute[] taxaAttribute = {new Attribute.Default<String>(XMLParser.IDREF, TaxaParser.TAXA)};
            if (taxonSets.size() > 0) {
                writer.writeOpenTag(CoalescentSimulatorParser.CONSTRAINED_TAXA);
                writer.writeTag(TaxaParser.TAXA, taxaAttribute, true);
                for (Taxa taxonSet : taxonSets) {
                    Parameter statistic = statistics.get(taxonSet);

                    Attribute mono = new Attribute.Default<Boolean>(CoalescentSimulatorParser.IS_MONOPHYLETIC, taxonSetsMono.get(taxonSet));

                    writer.writeOpenTag(CoalescentSimulatorParser.TMRCA_CONSTRAINT, mono);

                    writer.writeTag(TaxaParser.TAXA,
                            new Attribute[]{new Attribute.Default<String>(XMLParser.IDREF, taxonSet.getId())}, true);
                    if (statistic.isNodeHeight) {
                        if (statistic.priorType == PriorType.UNIFORM_PRIOR || statistic.priorType == PriorType.TRUNC_NORMAL_PRIOR) {
                            writer.writeOpenTag(UniformDistributionModelParser.UNIFORM_DISTRIBUTION_MODEL);
                            writer.writeTag(UniformDistributionModelParser.LOWER, new Attribute[]{}, "" + statistic.uniformLower, true);
                            writer.writeTag(UniformDistributionModelParser.UPPER, new Attribute[]{}, "" + statistic.uniformUpper, true);
                            writer.writeCloseTag(UniformDistributionModelParser.UNIFORM_DISTRIBUTION_MODEL);
                        }
                    }

                    writer.writeCloseTag(CoalescentSimulatorParser.TMRCA_CONSTRAINT);
                }
                writer.writeCloseTag(CoalescentSimulatorParser.CONSTRAINED_TAXA);
            } else {
                writer.writeTag(TaxaParser.TAXA, taxaAttribute, true);
            }

            writeInitialDemoModelRef(writer);
            writer.writeCloseTag(CoalescentSimulatorParser.COALESCENT_TREE);
        }
    }

    public void writeInitialDemoModelRef(XMLWriter writer) {
        if (nodeHeightPrior == CONSTANT) {
            writer.writeTag(ConstantPopulationModelParser.CONSTANT_POPULATION_MODEL, new Attribute[]{new Attribute.Default<String>(XMLParser.IDREF, "constant")}, true);
        } else if (nodeHeightPrior == EXPONENTIAL) {
            writer.writeTag(ExponentialGrowthModelParser.EXPONENTIAL_GROWTH_MODEL, new Attribute[]{new Attribute.Default<String>(XMLParser.IDREF, "exponential")}, true);
        } else {
            writer.writeTag(ConstantPopulationModelParser.CONSTANT_POPULATION_MODEL, new Attribute[]{new Attribute.Default<String>(XMLParser.IDREF, "initialDemo")}, true);
        }
    }

    public void writeNodeHeightPriorModelRef(XMLWriter writer) {
        if (nodeHeightPrior == CONSTANT) {
            writer.writeTag(ConstantPopulationModelParser.CONSTANT_POPULATION_MODEL, new Attribute[]{new Attribute.Default<String>(XMLParser.IDREF, "constant")}, true);
        } else if (nodeHeightPrior == EXPONENTIAL) {
            writer.writeTag(ExponentialGrowthModelParser.EXPONENTIAL_GROWTH_MODEL, new Attribute[]{new Attribute.Default<String>(XMLParser.IDREF, "exponential")}, true);
        } else if (nodeHeightPrior == LOGISTIC) {
            writer.writeTag(LogisticGrowthModelParser.LOGISTIC_GROWTH_MODEL, new Attribute[]{new Attribute.Default<String>(XMLParser.IDREF, "logistic")}, true);
        } else if (nodeHeightPrior == EXPANSION) {
            writer.writeTag(ExpansionModelParser.EXPANSION_MODEL, new Attribute[]{new Attribute.Default<String>(XMLParser.IDREF, "expansion")}, true);
        } else if (nodeHeightPrior == SKYLINE) {
            writer.writeTag(BayesianSkylineLikelihoodParser.SKYLINE_LIKELIHOOD, new Attribute[]{new Attribute.Default<String>(XMLParser.IDREF, "skyline")}, true);
        } else if (nodeHeightPrior == YULE) {
            writer.writeTag(YuleModelParser.YULE_MODEL, new Attribute[]{new Attribute.Default<String>(XMLParser.IDREF, "yule")}, true);
        } else if (nodeHeightPrior == BIRTH_DEATH) {
            writer.writeTag(BirthDeathGernhard08Model.BIRTH_DEATH_MODEL, new Attribute[]{new Attribute.Default<String>(XMLParser.IDREF, "birthDeath")}, true);
        } else {
            throw new RuntimeException("No coalescent model has been specified so cannot refer to it");
        }
    }

    /**
     * Generate XML for the user tree
     *
     * @param tree   the user tree
     * @param writer the writer
     */
    private void writeUserTree(Tree tree, XMLWriter writer) {

        writer.writeComment("The starting tree.");
        writer.writeOpenTag(
                "tree",
                new Attribute[]{
                        new Attribute.Default<String>("height", InitialTreeGenerator.STARTING_TREE),
                        new Attribute.Default<String>("usingDates", (maximumTipHeight > 0 ? "true" : "false"))
                }
        );
        writeNode(tree, tree.getRoot(), writer);
        writer.writeCloseTag("tree");
    }

    /**
     * Generate XML for the node of a user tree.
     *
     * @param tree   the user tree
     * @param node   the current node
     * @param writer the writer
     */
    private void writeNode(Tree tree, NodeRef node, XMLWriter writer) {

        writer.writeOpenTag(
                "node",
                new Attribute[]{new Attribute.Default<String>("height", "" + tree.getNodeHeight(node))}
        );

        if (tree.getChildCount(node) == 0) {
            writer.writeTag(TaxonParser.TAXON, new Attribute[]{new Attribute.Default<String>(XMLParser.IDREF, tree.getNodeTaxon(node).getId())}, true);
        }
        for (int i = 0; i < tree.getChildCount(node); i++) {
            writeNode(tree, tree.getChild(node, i), writer);
        }
        writer.writeCloseTag("node");
    }
}

