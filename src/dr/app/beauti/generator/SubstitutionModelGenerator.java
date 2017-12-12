/*
 * SubstitutionModelGenerator.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

import dr.app.beauti.options.*;
import dr.evomodel.substmodel.nucleotide.GTR;
import dr.evomodel.substmodel.nucleotide.NucModelType;
import dr.app.beauti.components.ComponentFactory;
import dr.app.beauti.types.FrequencyPolicyType;
import dr.app.beauti.types.MicroSatModelType;
import dr.app.beauti.util.XMLWriter;
import dr.evolution.alignment.Alignment;
import dr.evolution.datatype.DataType;
import dr.evomodelxml.substmodel.BinaryCovarionModelParser;
import dr.evomodelxml.substmodel.BinarySubstitutionModelParser;
import dr.evomodelxml.substmodel.EmpiricalAminoAcidModelParser;
import dr.evomodelxml.substmodel.FrequencyModelParser;
import dr.evomodelxml.substmodel.GTRParser;
import dr.evomodelxml.substmodel.GeneralSubstitutionModelParser;
import dr.evomodelxml.substmodel.HKYParser;
import dr.evomodelxml.substmodel.TN93Parser;
import dr.evomodelxml.siteratemodel.GammaSiteModelParser;
import dr.inference.model.StatisticParser;
import dr.oldevomodel.substmodel.AsymmetricQuadraticModel;
import dr.oldevomodel.substmodel.LinearBiasModel;
import dr.oldevomodel.substmodel.TwoPhaseModel;
import dr.evoxml.AlignmentParser;
import dr.evoxml.MicrosatelliteParser;
import dr.inference.model.ParameterParser;
import dr.oldevomodelxml.substmodel.*;
import dr.util.Attribute;
import dr.xml.XMLParser;

import java.util.List;

/**
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @author Walter Xie
 */
public class SubstitutionModelGenerator extends Generator {

    public SubstitutionModelGenerator(BeautiOptions options, ComponentFactory[] components) {
        super(options, components);
    }

    /**
     * Writes the substitution model to XML.
     *
     * @param writer the writer
     * @param model  the partition model to write in BEAST XML
     */
    public void writeSubstitutionSiteModel(PartitionSubstitutionModel model, XMLWriter writer) {

        DataType dataType = model.getDataType();
        String dataTypeDescription = dataType.getDescription();

        switch (dataType.getType()) {
            case DataType.NUCLEOTIDES:
                if (model.isUnlinkedSubstitutionModel()) {
                    for (int i = 1; i <= model.getCodonPartitionCount(); i++) {
                        switch (model.getNucSubstitutionModel()) {
                            case JC:
                                writeJCModel(i, writer, model);
                                break;
                            case HKY:
                                writeHKYModel(i, writer, model);
                                break;
                            case TN93:
                                writeTN93Model(i, writer, model);
                                break;
                            case GTR:
                                writeGTRModel(i, writer, model);
                                break;
                            default:
                                throw new IllegalArgumentException("unknown substition model type");
                        }
                    }
                } else {
                    switch (model.getNucSubstitutionModel()) {
                        case JC:
                            writeJCModel(-1, writer, model);
                            break;
                        case HKY:
                            writeHKYModel(-1, writer, model);
                            break;
                        case TN93:
                            writeTN93Model(-1, writer, model);
                            break;
                        case GTR:
                            writeGTRModel(-1, writer, model);
                            break;
                        default:
                            throw new IllegalArgumentException("unknown substitution model type");
                    }
                }

                //****************** Site Model *****************
                if (model.getCodonPartitionCount() > 1) { //model.getCodonHeteroPattern() != null) {
                    for (int i = 1; i <= model.getCodonPartitionCount(); i++) {
                        writeNucSiteModel(i, writer, model);
                    }
                    writer.println();
                } else {
                    writeNucSiteModel(-1, writer, model);
                }
                break;

            case DataType.AMINO_ACIDS:
                // Amino Acid model
                String aaModel = model.getAaSubstitutionModel().getXMLName();

                writer.writeComment("The " + aaModel + " substitution model");
                writer.writeTag(
                        EmpiricalAminoAcidModelParser.EMPIRICAL_AMINO_ACID_MODEL,
                        new Attribute[]{new Attribute.Default<String>(XMLParser.ID, model.getPrefix() + "aa"),
                                new Attribute.Default<String>("type", aaModel)}, true
                );

                //****************** Site Model *****************
                writeAASiteModel(writer, model);

                break;

            case DataType.TWO_STATES:
            case DataType.COVARION:

                switch (model.getBinarySubstitutionModel()) {
                    case BIN_DOLLO:
                        return;
                    case BIN_SIMPLE:
                        writeBinarySimpleModel(writer, model);
                        break;
                    case BIN_COVARION:
                        writeBinaryCovarionModel(writer, model);
                        break;
                }

                //****************** Site Model *****************
                writeTwoStateSiteModel(writer, model);
                break;

            case DataType.GENERAL:
            case DataType.CONTINUOUS:
                //handled by component
                break;

            case DataType.MICRO_SAT:
                writeMicrosatSubstModel(model, writer);
                break;

            default:
                throw new IllegalArgumentException("Unknown data type");
        }
    }

    /**
     * Write the JC model XML block.
     *
     * @param num    the model number
     * @param writer the writer
     * @param model  the partition model to write in BEAST XML
     */
    public void writeJCModel(int num, XMLWriter writer, PartitionSubstitutionModel model) {

        String prefix = model.getPrefix(num);

        writer.writeComment("The JC substitution model (Jukes & Cantor, 1969)");
        writer.writeOpenTag(NucModelType.HKY.getXMLName(),
                new Attribute[]{new Attribute.Default<String>(XMLParser.ID, prefix + "jc")}
        );
        writer.writeOpenTag(HKYParser.FREQUENCIES);
        writeFrequencyModelDNA(writer, model, num);
        writer.writeCloseTag(HKYParser.FREQUENCIES);

        writer.writeOpenTag(HKYParser.KAPPA);
        writeParameter("", 1, 1.0, Double.NaN, Double.NaN, writer);
        writer.writeCloseTag(HKYParser.KAPPA);

        writer.writeCloseTag(NucModelType.HKY.getXMLName());
    }

    /**
     * Write the HKY model XML block.
     *
     * @param num    the model number
     * @param writer the writer
     * @param model  the partition model to write in BEAST XML
     */
    public void writeHKYModel(int num, XMLWriter writer, PartitionSubstitutionModel model) {

        String prefix = model.getPrefix(num);

        // Hasegawa Kishino and Yano 85 model
        writer.writeComment("The HKY substitution model (Hasegawa, Kishino & Yano, 1985)");
        writer.writeOpenTag(NucModelType.HKY.getXMLName(),
                new Attribute[]{new Attribute.Default<String>(XMLParser.ID, prefix + "hky")}
        );
        writer.writeOpenTag(HKYParser.FREQUENCIES);
        writeFrequencyModelDNA(writer, model, num);
        writer.writeCloseTag(HKYParser.FREQUENCIES);

        writeParameter(num, HKYParser.KAPPA, "kappa", model, writer);
        writer.writeCloseTag(NucModelType.HKY.getXMLName());
    }

    /**
     * Write the TN93 model XML block.
     *
     * @param num    the model number
     * @param writer the writer
     * @param model  the partition model to write in BEAST XML
     */
    public void writeTN93Model(int num, XMLWriter writer, PartitionSubstitutionModel model) {

        String prefix = model.getPrefix(num);

        // TN93
        writer.writeComment("The TN93 substitution model");
        writer.writeOpenTag(NucModelType.TN93.getXMLName(),
                new Attribute[]{new Attribute.Default<String>(XMLParser.ID, prefix + "tn93")}
        );
        writer.writeOpenTag(HKYParser.FREQUENCIES);
        writeFrequencyModelDNA(writer, model, num);
        writer.writeCloseTag(HKYParser.FREQUENCIES);

        writeParameter(num, TN93Parser.KAPPA1, "kappa1", model, writer);
        writeParameter(num, TN93Parser.KAPPA2, "kappa2", model, writer);
        writer.writeCloseTag(NucModelType.TN93.getXMLName());
    }

    /**
     * Write the GTR model XML block.
     *
     * @param num    the model number
     * @param writer the writer
     * @param model  the partition model to write in BEAST XML
     */
    public void writeGTRModel(int num, XMLWriter writer, PartitionSubstitutionModel model) {

        String prefix = model.getPrefix(num);

        writer.writeComment("The general time reversible (GTR) substitution model");
        writer.writeOpenTag(GTRParser.GTR_MODEL,
                new Attribute[]{new Attribute.Default<String>(XMLParser.ID, prefix + "gtr")}
        );
        writer.writeOpenTag(GTRParser.FREQUENCIES);
        writeFrequencyModelDNA(writer, model, num);
        writer.writeCloseTag(GTRParser.FREQUENCIES);

        if (!options.classicOperatorsAndPriors && options.NEW_GTR_PARAMETERIZATION) {
            Parameter parameter = model.getParameter(model.getPrefixCodon(num) + PartitionSubstitutionModel.GTR_RATES);
            String prefix1 = model.getPrefix(num);
            writer.writeOpenTag(GTRParser.RATES);
            // fix the initial value to give maintained sum
            double initialValue = parameter.maintainedSum / parameter.dimension;
            writeParameter(prefix1 + PartitionSubstitutionModel.GTR_RATES, 6, initialValue, 0.0, Double.NaN, writer);
            writer.writeCloseTag(GTRParser.RATES);
        } else {
            writeParameter(num, GTRParser.A_TO_C, PartitionSubstitutionModel.GTR_RATE_NAMES[0], model, writer);
            writeParameter(num, GTRParser.A_TO_G, PartitionSubstitutionModel.GTR_RATE_NAMES[1], model, writer);
            writeParameter(num, GTRParser.A_TO_T, PartitionSubstitutionModel.GTR_RATE_NAMES[2], model, writer);
            writeParameter(num, GTRParser.C_TO_G, PartitionSubstitutionModel.GTR_RATE_NAMES[3], model, writer);
            writeParameter(num, GTRParser.G_TO_T, PartitionSubstitutionModel.GTR_RATE_NAMES[4], model, writer);
        }

        writer.writeCloseTag(GTRParser.GTR_MODEL);
    }

    // write frequencies for DNA data

    private void writeFrequencyModelDNA(XMLWriter writer, PartitionSubstitutionModel model, int num) {
        String dataTypeDescription = model.getDataType().getDescription();
        String prefix = model.getPrefix(num);
        writer.writeOpenTag(
                FrequencyModelParser.FREQUENCY_MODEL,
                new Attribute[]{
                        new Attribute.Default<String>("dataType", dataTypeDescription)
                }
        );

        writeAlignmentRefInFrequencies(writer, model, prefix, num);

        writer.writeOpenTag(FrequencyModelParser.FREQUENCIES);
        switch (model.getFrequencyPolicy()) {
            case ALLEQUAL:
            case ESTIMATED:
                if (num == -1 || model.isUnlinkedFrequencyModel()) { // single partition, or multiple partitions unlinked frequency
                    writer.writeTag(ParameterParser.PARAMETER, new Attribute[]{
                            new Attribute.Default<String>(XMLParser.ID, prefix + "frequencies"),
                            new Attribute.Default<String>(ParameterParser.VALUE, "0.25 0.25 0.25 0.25")}, true);
                } else { // multiple partitions but linked frequency
                    if (num == 1) {
                        writer.writeTag(ParameterParser.PARAMETER, new Attribute[]{
                                new Attribute.Default<String>(XMLParser.ID, model.getPrefix() + "frequencies"),
                                new Attribute.Default<String>(ParameterParser.VALUE, "0.25 0.25 0.25 0.25")}, true);
                    } else {
                        writer.writeIDref(ParameterParser.PARAMETER, model.getPrefix() + "frequencies");
                    }
                }
                break;

            case EMPIRICAL:
                if (num == -1 || model.isUnlinkedFrequencyModel()) { // single partition, or multiple partitions unlinked frequency
                    writeParameter(prefix + "frequencies", 4, Double.NaN, Double.NaN, Double.NaN, writer);
                } else { // multiple partitions but linked frequency
                    if (num == 1) {
                        writeParameter(model.getPrefix() + "frequencies", 4, Double.NaN, Double.NaN, Double.NaN, writer);
                    } else {
                        writer.writeIDref(ParameterParser.PARAMETER, model.getPrefix() + "frequencies");
                    }
                }
                break;
        }
        writer.writeCloseTag(FrequencyModelParser.FREQUENCIES);
        writer.writeCloseTag(FrequencyModelParser.FREQUENCY_MODEL);
    }

    // adding mergePatterns or alignment ref for EMPIRICAL

    private void writeAlignmentRefInFrequencies(XMLWriter writer, PartitionSubstitutionModel model, String prefix, int num) {
        if (model.getFrequencyPolicy() == FrequencyPolicyType.EMPIRICAL) {
            if (model.getDataType().getType() == DataType.NUCLEOTIDES && model.getCodonPartitionCount() > 1 && model.isUnlinkedSubstitutionModel()) {
                writeCodonPatternsRef(prefix, num, model.getCodonPartitionCount(), writer);

                // get the data partition for this substitution model.
               AbstractPartitionData partition = options.getDataPartitions(model).get(0);

               // for empirical frequencies use the entire alignment
               if (partition instanceof PartitionData) {
                   Alignment alignment = ((PartitionData)partition).getAlignment();
                   writer.writeIDref(AlignmentParser.ALIGNMENT, alignment.getId());
               } else {
                   throw new IllegalArgumentException("Partition is missing a data partition");
               }
            } else {
                for (AbstractPartitionData partition : options.getDataPartitions(model)) { //?
                    writer.writeIDref(AlignmentParser.ALIGNMENT, partition.getTaxonList().getId());
                }
            }
        }
    }

    /**
     * Write the Binary  simple model XML block.
     *
     * @param writer the writer
     * @param model  the partition model to write in BEAST XML
     */
    public void writeBinarySimpleModel(XMLWriter writer, PartitionSubstitutionModel model) {
        String dataTypeDescription = model.getDataType().getDescription();

        String prefix = model.getPrefix();

        writer.writeComment("The Binary simple model (based on the general substitution model)");
        writer.writeOpenTag(
                BinarySubstitutionModelParser.BINARY_SUBSTITUTION_MODEL,
                new Attribute[]{new Attribute.Default<String>(XMLParser.ID, prefix + "bsimple")}
        );
        writer.writeOpenTag(GeneralSubstitutionModelParser.FREQUENCIES);
        writer.writeOpenTag(
                FrequencyModelParser.FREQUENCY_MODEL,
                new Attribute[]{
                        new Attribute.Default<String>("dataType", dataTypeDescription)
                }
        );

        writeAlignmentRefInFrequencies(writer, model, prefix, -1);

        writeFrequencyModelBinary(writer, model, prefix);

        writer.writeCloseTag(FrequencyModelParser.FREQUENCY_MODEL);
        writer.writeCloseTag(GeneralSubstitutionModelParser.FREQUENCIES);

        writer.writeCloseTag(BinarySubstitutionModelParser.BINARY_SUBSTITUTION_MODEL);
    }

    /**
     * Write the Binary covarion model XML block
     *
     * @param writer the writer
     * @param model  the partition model to write
     */
    public void writeBinaryCovarionModel(XMLWriter writer, PartitionSubstitutionModel model) {
//        String dataTypeDescription = TwoStateCovarion.INSTANCE.getDescription(); // dataType="twoStateCovarion" for COVARION_MODEL
        String prefix = model.getPrefix();

        writer.writeComment("The Binary covarion model");
        writer.writeOpenTag(
                BinaryCovarionModelParser.COVARION_MODEL,
                new Attribute[]{new Attribute.Default<String>(XMLParser.ID, prefix + "bcov")}
        );

        // merge patterns then get frequencies.

        if (model.getFrequencyPolicy() == FrequencyPolicyType.EMPIRICAL) {
            List<AbstractPartitionData> partitions = options.getDataPartitions(model);
            Alignment alignment = ((PartitionData) partitions.get(0)).getAlignment();
//            Patterns patterns = new Patterns(partitions.get(0).getAlignment());
//            for (int i = 1; i < partitions.size(); i++) {
//                patterns.addPatterns(partitions.get(i).getAlignment());
//            }
            double[] frequencies = alignment.getStateFrequencies();
            writer.writeOpenTag(FrequencyModelParser.FREQUENCIES);
            writer.writeTag(ParameterParser.PARAMETER, new Attribute[]{
                    new Attribute.Default<String>(XMLParser.ID, prefix + "frequencies"),
                    new Attribute.Default<String>(ParameterParser.VALUE, frequencies[0] + " " + frequencies[1])}, true);
            writer.writeCloseTag(FrequencyModelParser.FREQUENCIES);

        } else {
            writeFrequencyModelBinary(writer, model, prefix);
        }

        writeParameter(BinaryCovarionModelParser.HIDDEN_FREQUENCIES,
                prefix + "hfrequencies", 2, 0.5, 0.0, 1.0, writer); // hfrequencies also 0.5 0.5

        writeParameter(BinaryCovarionModelParser.ALPHA, "bcov.alpha", model, writer);
        writeParameter(BinaryCovarionModelParser.SWITCHING_RATE, "bcov.s", model, writer);

        writer.writeCloseTag(BinaryCovarionModelParser.COVARION_MODEL);
    }

    // write frequencies for binary data

    private void writeFrequencyModelBinary(XMLWriter writer, PartitionSubstitutionModel model, String prefix) {
        writer.writeOpenTag(FrequencyModelParser.FREQUENCIES);
        switch (model.getFrequencyPolicy()) {
            case ALLEQUAL:
            case ESTIMATED:
                writer.writeTag(ParameterParser.PARAMETER, new Attribute[]{
                        new Attribute.Default<String>(XMLParser.ID, prefix + "frequencies"),
                        new Attribute.Default<String>(ParameterParser.VALUE, "0.5 0.5")}, true);
                break;

            case EMPIRICAL:
                writeParameter(prefix + "frequencies", 2, Double.NaN, Double.NaN, Double.NaN, writer);
                break;
        }
//        writeParameter(prefix + "frequencies", 2, Double.NaN, Double.NaN, Double.NaN, writer);
        writer.writeCloseTag(FrequencyModelParser.FREQUENCIES);
    }

    public void writeLog(PartitionSubstitutionModel model, XMLWriter writer) {

        int codonPartitionCount = model.getCodonPartitionCount();

        switch (model.getDataType().getType()) {
            case DataType.NUCLEOTIDES:

// THIS IS DONE BY ALLMUS logging in BeastGenerator
// single partition use clock.rate, no "mu"; multi-partition use "allmus"
//                if (codonPartitionCount > 1) {
//
//                    for (int i = 1; i <= codonPartitionCount; i++) {
//                        writer.writeIDref(ParameterParser.PARAMETER, model.getPrefix(i) + "mu");
//                    }
//                }
                switch (model.getNucSubstitutionModel()) {
                    case HKY:
                        if (codonPartitionCount > 1 && model.isUnlinkedSubstitutionModel()) {
                            for (int i = 1; i <= codonPartitionCount; i++) {
                                writer.writeIDref(ParameterParser.PARAMETER, model.getPrefix(i) + "kappa");
                            }
                        } else {
                            writer.writeIDref(ParameterParser.PARAMETER, model.getPrefix() + "kappa");
                        }
                        break;

                    case TN93:
                        if (codonPartitionCount > 1 && model.isUnlinkedSubstitutionModel()) {
                            for (int i = 1; i <= codonPartitionCount; i++) {
                                writer.writeIDref(ParameterParser.PARAMETER, model.getPrefix(i) + "kappa1");
                                writer.writeIDref(ParameterParser.PARAMETER, model.getPrefix(i) + "kappa2");
                            }
                        } else {
                            writer.writeIDref(ParameterParser.PARAMETER, model.getPrefix() + "kappa1");
                            writer.writeIDref(ParameterParser.PARAMETER, model.getPrefix() + "kappa2");
                        }
                        break;

                    case GTR:
                        if (codonPartitionCount > 1 && model.isUnlinkedSubstitutionModel()) {
                            for (int i = 1; i <= codonPartitionCount; i++) {
                                if (!options.classicOperatorsAndPriors && options.NEW_GTR_PARAMETERIZATION) {
                                    writer.writeIDref(ParameterParser.PARAMETER, model.getPrefix(i) + PartitionSubstitutionModel.GTR_RATES);
                                } else {
                                    for (String rateName : PartitionSubstitutionModel.GTR_RATE_NAMES) {
                                        writer.writeIDref(ParameterParser.PARAMETER, model.getPrefix(i) + rateName);
                                    }
                                }
                            }
                        } else {
                            for (String rateName : PartitionSubstitutionModel.GTR_RATE_NAMES) {
                                writer.writeIDref(ParameterParser.PARAMETER, model.getPrefix() + rateName);
                            }
                        }
                        break;
                }

                if (model.getFrequencyPolicy() == FrequencyPolicyType.ESTIMATED) {
                    if (codonPartitionCount > 1 && model.isUnlinkedSubstitutionModel() && model.isUnlinkedFrequencyModel()) {
                        for (int i = 1; i <= codonPartitionCount; i++) {
                            writer.writeIDref(ParameterParser.PARAMETER, model.getPrefix(i) + "frequencies");
                        }
                    } else {
                        writer.writeIDref(ParameterParser.PARAMETER, model.getPrefix() + "frequencies");
                    }
                }
                break;//NUCLEOTIDES

            case DataType.AMINO_ACIDS:
                break;//AMINO_ACIDS

            case DataType.TWO_STATES:
            case DataType.COVARION:
                String prefix = model.getPrefix();
                switch (model.getBinarySubstitutionModel()) {
                    case BIN_SIMPLE:
                    case BIN_DOLLO:
                        break;
                    case BIN_COVARION:
                        writeParameterRef(prefix + "bcov.alpha", writer);
                        writeParameterRef(prefix + "bcov.s", writer);
                        writeParameterRef(prefix + "frequencies", writer);
                        writeParameterRef(prefix + "hfrequencies", writer);
                        break;

                }
                break;//BINARY

            case DataType.GENERAL:
            case DataType.CONTINUOUS:
                // these datatypes are handled by components
                break;

            case DataType.MICRO_SAT:
                writeMicrosatSubstModelParameterRef(model, writer);
                break;

            default:
                throw new IllegalArgumentException("Unknown data type");
        }

        if (model.isGammaHetero()) {
            if (codonPartitionCount > 1 && model.isUnlinkedHeterogeneityModel()) {
                for (int i = 1; i <= codonPartitionCount; i++) {
                    writer.writeIDref(ParameterParser.PARAMETER, model.getPrefix(i) + "alpha");
                }
            } else {
                writer.writeIDref(ParameterParser.PARAMETER, model.getPrefix() + "alpha");
            }
        }

        if (model.isInvarHetero()) {
            if (codonPartitionCount > 1 && model.isUnlinkedHeterogeneityModel()) {
                for (int i = 1; i <= codonPartitionCount; i++) {
                    writer.writeIDref(ParameterParser.PARAMETER, model.getPrefix(i) + "pInv");
                }
            } else {
                writer.writeIDref(ParameterParser.PARAMETER, model.getPrefix() + "pInv");
            }
        }
    }

    public void writeMicrosatSubstModelParameterRef(PartitionSubstitutionModel model, XMLWriter writer) {
        if (model.getRatePorportion() == MicroSatModelType.RateProportionality.EQUAL_RATE) {

        } else if (model.getRatePorportion() == MicroSatModelType.RateProportionality.PROPORTIONAL_RATE) {
            writeParameterRef(model.getPrefix() + "propLinear", writer);
        } else if (model.getRatePorportion() == MicroSatModelType.RateProportionality.ASYM_QUAD) {

        }
        if (model.getMutationBias() == MicroSatModelType.MutationalBias.UNBIASED) {

        } else if (model.getMutationBias() == MicroSatModelType.MutationalBias.CONSTANT_BIAS) {
            writeParameterRef(model.getPrefix() + "biasConst", writer);
        } else if (model.getMutationBias() == MicroSatModelType.MutationalBias.LINEAR_BIAS) {
            writeParameterRef(model.getPrefix() + "biasConst", writer);
            writeParameterRef(model.getPrefix() + "biasLinear", writer);
        }
        if (model.getPhase() == MicroSatModelType.Phase.ONE_PHASE) {

        } else if (model.getPhase() == MicroSatModelType.Phase.TWO_PHASE) {
            writeParameterRef(model.getPrefix() + "geomDist", writer);
        } else if (model.getPhase() == MicroSatModelType.Phase.TWO_PHASE_STAR) {
            writeParameterRef(model.getPrefix() + "geomDist", writer);
            writeParameterRef(model.getPrefix() + "onePhaseProb", writer);
        }
    }

    /**
     * Write the nucleotide site model XML block.
     *
     * @param num    the model number
     * @param writer the writer
     * @param model  the partition model to write in BEAST XML
     */
    private void writeNucSiteModel(int num, XMLWriter writer, PartitionSubstitutionModel model) {

        String prefix = model.getPrefix(num);
        String prefix2 = model.getPrefix();

        writer.writeComment("site model");
        writer.writeOpenTag(GammaSiteModelParser.SITE_MODEL,
                new Attribute[]{new Attribute.Default<String>(XMLParser.ID, prefix + GammaSiteModelParser.SITE_MODEL)});


        writer.writeOpenTag(GammaSiteModelParser.SUBSTITUTION_MODEL);

        if (model.isUnlinkedSubstitutionModel()) {
            switch (model.getNucSubstitutionModel()) {

                // JC cannot be unlinked because it has no parameters
                case JC:
                    writer.writeIDref(NucModelType.HKY.getXMLName(), prefix + "jc");
                    break;
                case HKY:
                    writer.writeIDref(NucModelType.HKY.getXMLName(), prefix + "hky");
                    break;
                case GTR:
                    writer.writeIDref(GTRParser.GTR_MODEL, prefix + "gtr");
                    break;
                case TN93:
                    writer.writeIDref(NucModelType.TN93.getXMLName(), prefix + "tn93");
                    break;
                default:
                    throw new IllegalArgumentException("Unknown substitution model.");
            }

        } else {

            switch (model.getNucSubstitutionModel()) {
                case JC:
                    writer.writeIDref(NucModelType.HKY.getXMLName(), prefix2 + "jc");
                    break;
                case HKY:
                    writer.writeIDref(NucModelType.HKY.getXMLName(), prefix2 + "hky");
                    break;
                case GTR:
                    writer.writeIDref(GTRParser.GTR_MODEL, prefix2 + "gtr");
                    break;
                case TN93:
                    writer.writeIDref(NucModelType.TN93.getXMLName(), prefix2 + "tn93");
                    break;
                default:
                    throw new IllegalArgumentException("Unknown substitution model.");
            }
        }

        writer.writeCloseTag(GammaSiteModelParser.SUBSTITUTION_MODEL);

        if (!options.classicOperatorsAndPriors && options.NEW_RELATIVE_RATE_PARAMETERIZATION) {
            Parameter parameter;
            if (model.hasCodonPartitions()) {
                parameter = model.getParameter(model.getPrefixCodon(num) + "nu");
            } else {
                parameter = model.getParameter("nu");
            }
            if (parameter.getParent() != null && parameter.getParent().getSubParameters().size() > 0) {
                writeNuRelativeRateBlock(writer, prefix, parameter);
            }
        } else {
            if (model.hasCodonPartitions()) {
                writeParameter(num, dr.oldevomodelxml.sitemodel.GammaSiteModelParser.RELATIVE_RATE, "mu", model, writer);
            } else {
                writeParameter(dr.oldevomodelxml.sitemodel.GammaSiteModelParser.RELATIVE_RATE, "mu", model, writer);
            }
        }


        if (model.isGammaHetero()) {
            writer.writeOpenTag(GammaSiteModelParser.GAMMA_SHAPE, new Attribute.Default<String>(
                    GammaSiteModelParser.GAMMA_CATEGORIES, "" + model.getGammaCategories()));
            if (num == -1 || model.isUnlinkedHeterogeneityModel()) {
//                writeParameter(prefix + "alpha", model, writer);
                writeParameter(num, "alpha", model, writer);
            } else {
                // multiple partitions but linked heterogeneity
                if (num == 1) {
//                    writeParameter(prefix2 + "alpha", model, writer);
                    writeParameter("alpha", model, writer);
                } else {
                    writer.writeIDref(ParameterParser.PARAMETER, prefix2 + "alpha");
                }
            }
            writer.writeCloseTag(GammaSiteModelParser.GAMMA_SHAPE);
        }

        if (model.isInvarHetero()) {
            writer.writeOpenTag(GammaSiteModelParser.PROPORTION_INVARIANT);
            if (num == -1 || model.isUnlinkedHeterogeneityModel()) {
//                writeParameter(prefix + "pInv", model, writer);
                writeParameter(num, "pInv", model, writer);
            } else {
                // multiple partitions but linked heterogeneity
                if (num == 1) {
//                    writeParameter(prefix2 + "pInv", model, writer);
                    writeParameter("pInv", model, writer);
                } else {
                    writer.writeIDref(ParameterParser.PARAMETER, prefix2 + "pInv");
                }
            }
            writer.writeCloseTag(GammaSiteModelParser.PROPORTION_INVARIANT);
        }

        writer.writeCloseTag(GammaSiteModelParser.SITE_MODEL);

        if (!options.classicOperatorsAndPriors && options.NEW_RELATIVE_RATE_PARAMETERIZATION) {
            writeMuStatistic(writer, prefix, GammaSiteModelParser.SITE_MODEL);
        }

    }

    /**
     * Write the two states site model XML block.
     *
     * @param writer the writer
     * @param model  the partition model to write in BEAST XML
     */
    private void writeTwoStateSiteModel(XMLWriter writer, PartitionSubstitutionModel model) {

        String prefix = model.getPrefix();

        writer.writeComment("site model");
        writer.writeOpenTag(GammaSiteModelParser.SITE_MODEL,
                new Attribute[]{new Attribute.Default<String>(XMLParser.ID, prefix + GammaSiteModelParser.SITE_MODEL)});


        writer.writeOpenTag(GammaSiteModelParser.SUBSTITUTION_MODEL);

        switch (model.getBinarySubstitutionModel()) {
            case BIN_SIMPLE:
                //writer.writeIDref(dr.evomodel.substmodel.GeneralSubstitutionModel.GENERAL_SUBSTITUTION_MODEL, "bsimple");
                writer.writeIDref(BinarySubstitutionModelParser.BINARY_SUBSTITUTION_MODEL, prefix + "bsimple");
                break;
            case BIN_COVARION:
                writer.writeIDref(BinaryCovarionModelParser.COVARION_MODEL, prefix + "bcov");
                break;
            case BIN_DOLLO:
                break; // Handled by component
            default:
                throw new IllegalArgumentException("Unknown substitution model.");
        }

        writer.writeCloseTag(GammaSiteModelParser.SUBSTITUTION_MODEL);

        if (!options.classicOperatorsAndPriors && options.NEW_RELATIVE_RATE_PARAMETERIZATION) {
            Parameter parameter = model.getParameter("nu");
            String prefix1 = options.getPrefix();
            if (parameter.getSubParameters().size() > 0) {
                writeNuRelativeRateBlock(writer, prefix1, parameter);
            }
        } else {
            writeParameter(GammaSiteModelParser.RELATIVE_RATE, "mu", model, writer);
        }

        if (model.isGammaHetero()) {
            writer.writeOpenTag(GammaSiteModelParser.GAMMA_SHAPE,
                    new Attribute.Default<String>(GammaSiteModelParser.GAMMA_CATEGORIES, "" + model.getGammaCategories()));
            writeParameter(prefix + "alpha", model, writer);
            writer.writeCloseTag(GammaSiteModelParser.GAMMA_SHAPE);
        }

        if (model.isInvarHetero()) {
            writeParameter(GammaSiteModelParser.PROPORTION_INVARIANT, "pInv", model, writer);
        }

        writer.writeCloseTag(GammaSiteModelParser.SITE_MODEL);

        if (!options.classicOperatorsAndPriors && options.NEW_RELATIVE_RATE_PARAMETERIZATION) {
            writeMuStatistic(writer, prefix, GammaSiteModelParser.SITE_MODEL);
        }

    }

    /**
     * Write the AA site model XML block.
     *
     * @param writer the writer
     * @param model  the partition model to write in BEAST XML
     */
    private void writeAASiteModel(XMLWriter writer, PartitionSubstitutionModel model) {

        String prefix = model.getPrefix();

        writer.writeComment("site model");
        writer.writeOpenTag(GammaSiteModelParser.SITE_MODEL, new Attribute[]{
                new Attribute.Default<String>(XMLParser.ID, prefix + GammaSiteModelParser.SITE_MODEL)});


        writer.writeOpenTag(GammaSiteModelParser.SUBSTITUTION_MODEL);
        writer.writeIDref(EmpiricalAminoAcidModelParser.EMPIRICAL_AMINO_ACID_MODEL, prefix + "aa");
        writer.writeCloseTag(GammaSiteModelParser.SUBSTITUTION_MODEL);

        if (!options.classicOperatorsAndPriors && options.NEW_RELATIVE_RATE_PARAMETERIZATION) {
            Parameter parameter = model.getParameter("nu");
            String prefix1 = options.getPrefix();
            if (parameter.getSubParameters().size() > 0) {
                writeNuRelativeRateBlock(writer, prefix1, parameter);
            }

        } else {
            writeParameter(GammaSiteModelParser.RELATIVE_RATE, "mu", model, writer);
        }

        if (model.isGammaHetero()) {
            writer.writeOpenTag(GammaSiteModelParser.GAMMA_SHAPE,
                    new Attribute.Default<String>(
                            GammaSiteModelParser.GAMMA_CATEGORIES, "" + model.getGammaCategories()));
            writeParameter("alpha", model, writer);
            writer.writeCloseTag(GammaSiteModelParser.GAMMA_SHAPE);
        }

        if (model.isInvarHetero()) {
            writeParameter(GammaSiteModelParser.PROPORTION_INVARIANT, "pInv", model, writer);
        }

        writer.writeCloseTag(GammaSiteModelParser.SITE_MODEL);

        if (!options.classicOperatorsAndPriors && options.NEW_RELATIVE_RATE_PARAMETERIZATION) {
            writeMuStatistic(writer, prefix, GammaSiteModelParser.SITE_MODEL);
        }


    }

    /**
     * Write the relative rate block for site model XML block.
     *
     * @param writer the writer
     */
    private void writeNuRelativeRateBlock(XMLWriter writer, String prefix, Parameter parameter) {
        double weight = ((double) parameter.getParent().getDimensionWeight()) / parameter.getDimensionWeight();
        writer.writeOpenTag(GammaSiteModelParser.RELATIVE_RATE,
                new Attribute.Default<String>(GammaSiteModelParser.WEIGHT, "" + weight));
        // Initial values must sum to 1.0
        double initial = 1.0 / parameter.getParent().getSubParameters().size();
        writeParameter(prefix + "nu", 1, initial, 0.0, 1.0, writer);
        writer.writeCloseTag(GammaSiteModelParser.RELATIVE_RATE);
    }

    /**
     * Write a statistic for mu from a site model (when parameterised using nu).
     *
     * @param writer the writer
     */
    private void writeMuStatistic(XMLWriter writer, String prefix, String siteModelTag) {
        writer.writeComment("");
        writer.writeOpenTag(StatisticParser.STATISTIC, new Attribute[]{
                new Attribute.Default<String>(XMLParser.ID, prefix + "mu"),
                new Attribute.Default<String>(StatisticParser.NAME, "mu")});
        writer.writeIDref(siteModelTag, prefix + siteModelTag);
        writer.writeCloseTag(StatisticParser.STATISTIC);
    }

    private void writeMicrosatSubstModel(PartitionSubstitutionModel model, XMLWriter writer) {

        writer.writeOpenTag(AsymmetricQuadraticModel.ASYMQUAD_MODEL, new Attribute[]{
                new Attribute.Default<String>(XMLParser.ID, model.getPrefix() + AsymmetricQuadraticModel.ASYMQUAD_MODEL),
                new Attribute.Default<Boolean>(AsymQuadModelParser.IS_SUBMODEL,
                        !(model.getMutationBias() == MicroSatModelType.MutationalBias.UNBIASED
                                && model.getPhase() == MicroSatModelType.Phase.ONE_PHASE)), // ?U1 is false
        });
        writer.writeIDref(MicrosatelliteParser.MICROSAT, model.getMicrosatellite().getName());

//        if (model.getRatePorportion() == MicroSatModelType.RateProportionality.EQUAL_RATE) {
//            // no xml
//        } else
        if (model.getRatePorportion() == MicroSatModelType.RateProportionality.PROPORTIONAL_RATE) {
            writeParameter(AsymQuadModelParser.EXPANSION_LIN, "propLinear", model, writer);
            writeParameterRef(AsymQuadModelParser.CONTRACTION_LIN, model.getPrefix() + "propLinear", writer);
        } else if (model.getRatePorportion() == MicroSatModelType.RateProportionality.ASYM_QUAD) {

        }
        writer.writeCloseTag(AsymmetricQuadraticModel.ASYMQUAD_MODEL);

        //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        if (model.getMutationBias() != MicroSatModelType.MutationalBias.UNBIASED) {
            writer.writeOpenTag(LinearBiasModel.LINEAR_BIAS_MODEL, new Attribute[]{
                    new Attribute.Default<String>(XMLParser.ID, model.getPrefix() + LinearBiasModel.LINEAR_BIAS_MODEL),
                    new Attribute.Default<Boolean>(LinearBiasModelParser.LOGISTICS, true),
                    new Attribute.Default<Boolean>(LinearBiasModelParser.ESTIMATE_SUBMODEL_PARAMS,
                            model.getMutationBias() == MicroSatModelType.MutationalBias.LINEAR_BIAS),
                    new Attribute.Default<Boolean>(LinearBiasModelParser.IS_SUBMODEL,
                            model.getPhase() != MicroSatModelType.Phase.ONE_PHASE),
            });
            writer.writeIDref(MicrosatelliteParser.MICROSAT, model.getMicrosatellite().getName());

//            if (model.getMutationBias() == MicroSatModelType.MutationalBias.CONSTANT_BIAS)
            writeParameterRef(LinearBiasModelParser.SUBMODEL, model.getPrefix() + AsymmetricQuadraticModel.ASYMQUAD_MODEL, writer);

            writeParameter(LinearBiasModelParser.BIAS_CONSTANT, "biasConst", model, writer);

            if (model.getMutationBias() == MicroSatModelType.MutationalBias.LINEAR_BIAS) {
                writeParameter(LinearBiasModelParser.BIAS_LINEAR, "biasLinear", model, writer);
            }
            writer.writeCloseTag(LinearBiasModel.LINEAR_BIAS_MODEL);
        }

        //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        if (model.getPhase() != MicroSatModelType.Phase.ONE_PHASE) {
            writer.writeOpenTag(TwoPhaseModel.TWO_PHASE_MODEL, new Attribute[]{
                    new Attribute.Default<String>(XMLParser.ID, model.getPrefix() + TwoPhaseModel.TWO_PHASE_MODEL),
                    new Attribute.Default<Boolean>(TwoPhaseModelParser.ESTIMATE_SUBMODEL_PARAMS, true),
            });
            writer.writeIDref(MicrosatelliteParser.MICROSAT, model.getMicrosatellite().getName());

            if (model.getMutationBias() == MicroSatModelType.MutationalBias.UNBIASED) {
                writeParameterRef(TwoPhaseModelParser.SUBMODEL, model.getPrefix() + AsymmetricQuadraticModel.ASYMQUAD_MODEL, writer);
            } else {
                writeParameterRef(TwoPhaseModelParser.SUBMODEL, model.getPrefix() + LinearBiasModel.LINEAR_BIAS_MODEL, writer);
            }

            if (model.getPhase() == MicroSatModelType.Phase.TWO_PHASE) {
                writeParameter(TwoPhaseModelParser.GEO_PARAM, "geomDist", model, writer);
                writer.writeOpenTag(TwoPhaseModelParser.ONEPHASEPR_PARAM);
                writeParameter(model.getPrefix() + "onePhaseProb", 1, 0.0, Double.NaN, Double.NaN, writer);
                writer.writeCloseTag(TwoPhaseModelParser.ONEPHASEPR_PARAM);
            } else if (model.getPhase() == MicroSatModelType.Phase.TWO_PHASE_STAR) {
                writeParameter(TwoPhaseModelParser.GEO_PARAM, "geomDist", model, writer);
                writeParameter(TwoPhaseModelParser.ONEPHASEPR_PARAM, "onePhaseProb", model, writer);
            }
            writer.writeCloseTag(TwoPhaseModel.TWO_PHASE_MODEL);
        }
    }

}