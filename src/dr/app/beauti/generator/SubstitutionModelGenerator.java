package dr.app.beauti.generator;

import dr.app.beauti.components.ComponentFactory;
import dr.app.beauti.options.AbstractPartitionData;
import dr.app.beauti.options.BeautiOptions;
import dr.app.beauti.options.PartitionData;
import dr.app.beauti.options.PartitionSubstitutionModel;
import dr.app.beauti.types.DiscreteSubstModelType;
import dr.app.beauti.types.FrequencyPolicyType;
import dr.app.beauti.types.MicroSatModelType;
import dr.app.beauti.util.XMLWriter;
import dr.evolution.alignment.Alignment;
import dr.evolution.datatype.DataType;
import dr.evolution.datatype.Nucleotides;
import dr.evomodel.sitemodel.GammaSiteModel;
import dr.evomodel.sitemodel.SiteModel;
import dr.evomodel.substmodel.*;
import dr.evomodelxml.sitemodel.GammaSiteModelParser;
import dr.evomodelxml.substmodel.*;
import dr.evoxml.AlignmentParser;
import dr.evoxml.GeneralDataTypeParser;
import dr.evoxml.MicrosatelliteParser;
import dr.inference.model.ParameterParser;
import dr.inferencexml.loggers.ColumnsParser;
import dr.inferencexml.model.CompoundParameterParser;
import dr.inferencexml.model.ProductStatisticParser;
import dr.inferencexml.model.SumStatisticParser;
import dr.util.Attribute;
import dr.xml.XMLParser;

import java.util.List;

/**
 * @author Alexei Drummond
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
                // Jukes-Cantor model
                if (model.getNucSubstitutionModel() == NucModelType.JC) {
                    String prefix = model.getPrefix();
                    writer.writeComment("The JC substitution model (Jukes & Cantor, 1969)");
                    writer.writeOpenTag(NucModelType.HKY.getXMLName(),
                            new Attribute[]{new Attribute.Default<String>(XMLParser.ID, prefix + "jc")}
                    );
                    writer.writeOpenTag(HKYParser.FREQUENCIES);
                    writer.writeOpenTag(
                            FrequencyModelParser.FREQUENCY_MODEL,
                            new Attribute[]{
                                    new Attribute.Default<String>("dataType", dataTypeDescription)
                            }
                    );
                    writer.writeOpenTag(FrequencyModelParser.FREQUENCIES);
                    writer.writeTag(
                            ParameterParser.PARAMETER,
                            new Attribute[]{
                                    new Attribute.Default<String>(XMLParser.ID, prefix + "frequencies"),
                                    new Attribute.Default<String>(ParameterParser.VALUE, "0.25 0.25 0.25 0.25")
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

                    throw new IllegalArgumentException("AR: Need to check that kappa = 1 for JC (I have feeling it should be 0.5)");

                } else {
                    // Hasegawa Kishino and Yano 85 model
                    if (model.getNucSubstitutionModel() == NucModelType.HKY) {
                        if (model.isUnlinkedSubstitutionModel()) {
                            for (int i = 1; i <= model.getCodonPartitionCount(); i++) {
                                writeHKYModel(i, writer, model);
                            }
                        } else {
                            writeHKYModel(-1, writer, model);
                        }

                    } else if (model.getNucSubstitutionModel() == NucModelType.TN93) {
                        if (model.isUnlinkedSubstitutionModel()) {
                            for (int i = 1; i <= model.getCodonPartitionCount(); i++) {
                                writeTN93Model(i, writer, model);
                            }
                        } else {
                            writeTN93Model(-1, writer, model);
                        }

                    } else {
                        // General time reversible model
                        if (model.getNucSubstitutionModel() == NucModelType.GTR) {
                            if (model.isUnlinkedSubstitutionModel()) {
                                for (int i = 1; i <= model.getCodonPartitionCount(); i++) {
                                    writeGTRModel(i, writer, model);
                                }
                            } else {
                                writeGTRModel(-1, writer, model);
                            }
                        }
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
                writeDiscreteTraitsSubstModel(model, writer);
                //****************** Site Model *****************
                writeDiscreteTraitsSiteModel(model, writer);
                break;

            case DataType.MICRO_SAT:
                writeMicrosatSubstModel(model, writer);
                break;

            default:
                throw new IllegalArgumentException("Unknown data type");
        }
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

        writeParameter(num, GTRParser.A_TO_C, PartitionSubstitutionModel.GTR_RATE_NAMES[0], model, writer);
        writeParameter(num, GTRParser.A_TO_G, PartitionSubstitutionModel.GTR_RATE_NAMES[1], model, writer);
        writeParameter(num, GTRParser.A_TO_T, PartitionSubstitutionModel.GTR_RATE_NAMES[2], model, writer);
        writeParameter(num, GTRParser.C_TO_G, PartitionSubstitutionModel.GTR_RATE_NAMES[3], model, writer);
        writeParameter(num, GTRParser.G_TO_T, PartitionSubstitutionModel.GTR_RATE_NAMES[4], model, writer);
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
            if (model.getDataType() == Nucleotides.INSTANCE && model.getCodonPartitionCount() > 1 && model.isUnlinkedSubstitutionModel()) {
                for (AbstractPartitionData partition : options.getAllPartitionData(model)) { //?
                    if (num >= 0)
                        writeCodonPatternsRef(prefix + partition.getPrefix(), num, model.getCodonPartitionCount(), writer);
                }
            } else {
                for (AbstractPartitionData partition : options.getAllPartitionData(model)) { //?
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
            List<AbstractPartitionData> partitions = options.getAllPartitionData(model);
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

    /**
     * Discrete Traits Subst Model
     *
     * @param model  PartitionSubstitutionModel
     * @param writer XMLWriter
     */
    private void writeDiscreteTraitsSubstModel(PartitionSubstitutionModel model, XMLWriter writer) {

        int stateCount = options.getStatesForDiscreteModel(model).size();

        if (model.getDiscreteSubstType() == DiscreteSubstModelType.SYM_SUBST) {
            writer.writeComment("symmetric CTMC model for discrete state reconstructions");

            writer.writeOpenTag(GeneralSubstitutionModelParser.GENERAL_SUBSTITUTION_MODEL, new Attribute[]{
                    new Attribute.Default<String>(XMLParser.ID, model.getPrefix() + AbstractSubstitutionModel.MODEL)});

            writer.writeIDref(GeneralDataTypeParser.GENERAL_DATA_TYPE, model.getPrefix() + DiscreteTraitGenerator.DATA_TYPE);

            writer.writeOpenTag(GeneralSubstitutionModelParser.FREQUENCIES);

            writeDiscreteFrequencyModel(model, stateCount, true, writer);

            writer.writeCloseTag(GeneralSubstitutionModelParser.FREQUENCIES);

            //---------------- rates and indicators -----------------

            // AR - we are trying to unify this setup. The only difference between BSSVS and not is the presence of indicators...
//            if (!model.isActivateBSSVS()) {
//                writer.writeComment("Rates parameter in " + GeneralSubstitutionModelParser.GENERAL_SUBSTITUTION_MODEL
//                        + " element should have (" + (numOfStates * (numOfStates - 1) / 2) + " - 1) dimensions");
//                writeRatesAndIndicators(model, (numOfStates * (numOfStates - 1) / 2) - 1, null, writer);
//            } else {
//                writeRatesAndIndicators(model, numOfStates * (numOfStates - 1) / 2, null, writer);
//            }
            writeRatesAndIndicators(model, stateCount * (stateCount - 1) / 2, null, writer);
            writer.writeCloseTag(GeneralSubstitutionModelParser.GENERAL_SUBSTITUTION_MODEL);

        } else if (model.getDiscreteSubstType() == DiscreteSubstModelType.ASYM_SUBST) {
            writer.writeComment("asymmetric CTMC model for discrete state reconstructions");

            writer.writeOpenTag(GeneralSubstitutionModelParser.GENERAL_SUBSTITUTION_MODEL, new Attribute[]{
                    new Attribute.Default<String>(XMLParser.ID, model.getPrefix() + AbstractSubstitutionModel.MODEL),
                    new Attribute.Default<Boolean>(ComplexSubstitutionModelParser.RANDOMIZE, false)});

            writer.writeIDref(GeneralDataTypeParser.GENERAL_DATA_TYPE, model.getPrefix() + DiscreteTraitGenerator.DATA_TYPE);

            writer.writeOpenTag(GeneralSubstitutionModelParser.FREQUENCIES);

            writeDiscreteFrequencyModel(model, stateCount, true, writer);

            writer.writeCloseTag(GeneralSubstitutionModelParser.FREQUENCIES);

            //---------------- rates and indicators -----------------
            writeRatesAndIndicators(model, stateCount * (stateCount - 1), null, writer);

            writer.writeCloseTag(GeneralSubstitutionModelParser.GENERAL_SUBSTITUTION_MODEL);

        } else {

        }

        if (model.isActivateBSSVS()) // If "BSSVS" is not activated, rateIndicator should not be there.
            writeStatisticModel(model, writer);
    }

    private void writeDiscreteFrequencyModel(PartitionSubstitutionModel model, int stateCount, Boolean normalize, XMLWriter writer) {
        if (normalize == null) {
            writer.writeOpenTag(FrequencyModelParser.FREQUENCY_MODEL, new Attribute[]{
                    new Attribute.Default<String>(XMLParser.ID, model.getPrefix() + FrequencyModelParser.FREQUENCY_MODEL)});
        } else {
            writer.writeOpenTag(FrequencyModelParser.FREQUENCY_MODEL, new Attribute[]{
                    new Attribute.Default<String>(XMLParser.ID, model.getPrefix() + FrequencyModelParser.FREQUENCY_MODEL),
                    new Attribute.Default<Boolean>(FrequencyModelParser.NORMALIZE, normalize)});
        }

        writer.writeIDref(GeneralDataTypeParser.GENERAL_DATA_TYPE, model.getPrefix() + DiscreteTraitGenerator.DATA_TYPE);

        writer.writeOpenTag(FrequencyModelParser.FREQUENCIES);
        writeParameter(model.getPrefix() + "trait.frequencies", stateCount, Double.NaN, Double.NaN, Double.NaN, writer);
        writer.writeCloseTag(FrequencyModelParser.FREQUENCIES);

        writer.writeCloseTag(FrequencyModelParser.FREQUENCY_MODEL);
    }

    private void writeRatesAndIndicators(PartitionSubstitutionModel model, int dimension, Integer relativeTo, XMLWriter writer) {
        writer.writeComment("rates and indicators");

        if (relativeTo == null) {
            writer.writeOpenTag(GeneralSubstitutionModelParser.RATES);
        } else {
            writer.writeOpenTag(GeneralSubstitutionModelParser.RATES, new Attribute[]{
                    new Attribute.Default<Integer>(GeneralSubstitutionModelParser.RELATIVE_TO, relativeTo)});
        }

        model.getParameter("trait.rates").isFixed = true;
        writeParameter(model.getParameter("trait.rates"), dimension, writer);

        writer.writeCloseTag(GeneralSubstitutionModelParser.RATES);

        if (model.isActivateBSSVS()) { //If "BSSVS" is not activated, rateIndicator should not be there.
            writer.writeOpenTag(GeneralSubstitutionModelParser.INDICATOR);
            model.getParameter("trait.indicators").isFixed = true;
            writeParameter(model.getParameter("trait.indicators"), dimension, writer);
            writer.writeCloseTag(GeneralSubstitutionModelParser.INDICATOR);
        }

    }

    private void writeStatisticModel(PartitionSubstitutionModel model, XMLWriter writer) {
        writer.writeOpenTag(SumStatisticParser.SUM_STATISTIC, new Attribute[]{
                new Attribute.Default<String>(XMLParser.ID, model.getPrefix() + "trait.nonZeroRates"),
                new Attribute.Default<Boolean>(SumStatisticParser.ELEMENTWISE, true)});
        writer.writeIDref(ParameterParser.PARAMETER, model.getPrefix() + "trait.indicators");
        writer.writeCloseTag(SumStatisticParser.SUM_STATISTIC);

        writer.writeOpenTag(ProductStatisticParser.PRODUCT_STATISTIC, new Attribute[]{
                new Attribute.Default<String>(XMLParser.ID, model.getPrefix() + "actualRates"),
                new Attribute.Default<Boolean>(SumStatisticParser.ELEMENTWISE, false)});
        writer.writeIDref(ParameterParser.PARAMETER, model.getPrefix() + "trait.indicators");
        writer.writeIDref(ParameterParser.PARAMETER, model.getPrefix() + "trait.rates");
        writer.writeCloseTag(ProductStatisticParser.PRODUCT_STATISTIC);
    }

    /**
     * Write the site model XML block.
     *
     * @param model            the partition model to write in BEAST XML
     * @param writer           the writer
     */
//    public void writeSiteModel(PartitionSubstitutionModel model, XMLWriter writer) {
//
//        switch (model.getDataType().getType()) {
//            case DataType.NUCLEOTIDES:
//                if (model.getCodonPartitionCount() > 1) { //model.getCodonHeteroPattern() != null) {
//                    for (int i = 1; i <= model.getCodonPartitionCount(); i++) {
//                        writeNucSiteModel(i, writer, model);
//                    }
//                    writer.println();
//                } else {
//                    writeNucSiteModel(-1, writer, model);
//                }
//                break;
//
//            case DataType.AMINO_ACIDS:
//                writeAASiteModel(writer, model);
//                break;
//
//            case DataType.TWO_STATES:
//            case DataType.COVARION:
//                writeTwoStateSiteModel(writer, model);
//                break;
//
//            default:
//                throw new IllegalArgumentException("Unknown data type");
//        }
//    }

    /**
     * Write the allMus for each partition model.
     *
     * @param model  PartitionSubstitutionModel
     * @param writer XMLWriter
     */
    public void writeAllMus(PartitionSubstitutionModel model, XMLWriter writer) {
        if (model.hasCodon()) { // write allMus for codon model
            // allMus is global for each gene
            writer.writeOpenTag(CompoundParameterParser.COMPOUND_PARAMETER,
                    new Attribute[]{new Attribute.Default<String>(XMLParser.ID, model.getPrefix() + "allMus")});

            writeMuParameterRefs(model, writer);

            writer.writeCloseTag(CompoundParameterParser.COMPOUND_PARAMETER);
        }
    }

    /**
     * Write the all the mu parameters for this partition model.
     *
     * @param writer the writer
     * @param model  the partition model to write in BEAST XML
     */
    public void writeMuParameterRefs(PartitionSubstitutionModel model, XMLWriter writer) {

        if (model.getDataType().getType() == DataType.NUCLEOTIDES && model.getCodonHeteroPattern() != null) {
            for (int i = 1; i <= model.getCodonPartitionCount(); i++) {
                writer.writeIDref(ParameterParser.PARAMETER, model.getPrefix(i) + "mu");
            }
        } else {
            writer.writeIDref(ParameterParser.PARAMETER, model.getPrefix() + "mu");
        }

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
                                for (String rateName : PartitionSubstitutionModel.GTR_RATE_NAMES) {
                                    writer.writeIDref(ParameterParser.PARAMETER, model.getPrefix(i) + rateName);
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
                //TODO
                break;

            case DataType.MICRO_SAT:
                writeMicrosatSubstModelRef(model, writer);
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

    public void writeRateLog(PartitionSubstitutionModel model, XMLWriter writer) {
        writer.writeIDref(ParameterParser.PARAMETER, model.getPrefix() + "trait.rates");

        if (model.isActivateBSSVS()) { //If "BSSVS" is not activated, rateIndicator should not be there.
            writer.writeIDref(ParameterParser.PARAMETER, model.getPrefix() + "trait.indicators");
            writer.writeIDref(SumStatisticParser.SUM_STATISTIC, model.getPrefix() + "trait.nonZeroRates");
        }
    }

    public void writeStatisticLog(PartitionSubstitutionModel model, XMLWriter writer) {
        if (model.isActivateBSSVS()) { //If "BSSVS" is not activated, rateIndicator should not be there.
            writer.writeOpenTag(ColumnsParser.COLUMN,
                    new Attribute[]{
                            new Attribute.Default<String>(ColumnsParser.LABEL, model.getPrefix() + "nonZeroRates"),
                            new Attribute.Default<String>(ColumnsParser.SIGNIFICANT_FIGURES, "6"),
                            new Attribute.Default<String>(ColumnsParser.WIDTH, "12")
                    }
            );

            writer.writeIDref(SumStatisticParser.SUM_STATISTIC, model.getPrefix() + "trait.nonZeroRates");

            writer.writeCloseTag(ColumnsParser.COLUMN);
        }
    }

    public void writeMicrosatSubstModelRef(PartitionSubstitutionModel model, XMLWriter writer) {
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
        writer.writeOpenTag(GammaSiteModel.SITE_MODEL,
                new Attribute[]{new Attribute.Default<String>(XMLParser.ID, prefix + SiteModel.SITE_MODEL)});


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

        if (model.hasCodon()) {
            writeParameter(num, GammaSiteModelParser.RELATIVE_RATE, "mu", model, writer);
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

        writer.writeCloseTag(GammaSiteModel.SITE_MODEL);
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
        writer.writeOpenTag(GammaSiteModel.SITE_MODEL,
                new Attribute[]{new Attribute.Default<String>(XMLParser.ID, prefix + SiteModel.SITE_MODEL)});


        writer.writeOpenTag(GammaSiteModelParser.SUBSTITUTION_MODEL);

        switch (model.getBinarySubstitutionModel()) {
            case BIN_SIMPLE:
                //writer.writeIDref(dr.evomodel.substmodel.GeneralSubstitutionModel.GENERAL_SUBSTITUTION_MODEL, "bsimple");
                writer.writeIDref(BinarySubstitutionModelParser.BINARY_SUBSTITUTION_MODEL, prefix + "bsimple");
                break;
            case BIN_COVARION:
                writer.writeIDref(BinaryCovarionModelParser.COVARION_MODEL, prefix + "bcov");
                break;
            default:
                throw new IllegalArgumentException("Unknown substitution model.");
        }

        writer.writeCloseTag(GammaSiteModelParser.SUBSTITUTION_MODEL);

        if (model.hasCodon()) {
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

        writer.writeCloseTag(GammaSiteModel.SITE_MODEL);
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
        writer.writeOpenTag(GammaSiteModel.SITE_MODEL, new Attribute[]{
                new Attribute.Default<String>(XMLParser.ID, prefix + SiteModel.SITE_MODEL)});


        writer.writeOpenTag(GammaSiteModelParser.SUBSTITUTION_MODEL);
        writer.writeIDref(EmpiricalAminoAcidModelParser.EMPIRICAL_AMINO_ACID_MODEL, prefix + "aa");
        writer.writeCloseTag(GammaSiteModelParser.SUBSTITUTION_MODEL);

        if (model.hasCodon()) {
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

        writer.writeCloseTag(GammaSiteModel.SITE_MODEL);
    }

    private void writeDiscreteTraitsSiteModel(PartitionSubstitutionModel model, XMLWriter writer) {
        writer.writeOpenTag(SiteModel.SITE_MODEL, new Attribute[]{
                new Attribute.Default<String>(XMLParser.ID, model.getPrefix() + SiteModel.SITE_MODEL)});

        writer.writeOpenTag(GammaSiteModelParser.SUBSTITUTION_MODEL);
        writer.writeIDref(GeneralSubstitutionModelParser.GENERAL_SUBSTITUTION_MODEL, model.getPrefix() + AbstractSubstitutionModel.MODEL);
        writer.writeCloseTag(GammaSiteModelParser.SUBSTITUTION_MODEL);

//        writer.writeOpenTag(GammaSiteModelParser.MUTATION_RATE);
//        writeParameter(model.getParameter("trait.mu"), -1, writer);
//        writer.writeCloseTag(GammaSiteModelParser.MUTATION_RATE);

        writer.writeCloseTag(SiteModel.SITE_MODEL);
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