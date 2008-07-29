package dr.app.beauti.generator;

import dr.app.beauti.XMLWriter;
import dr.app.beauti.options.BeautiOptions;
import dr.app.beauti.options.ModelOptions;
import dr.app.beauti.options.NucModelType;
import dr.app.beauti.options.PartitionModel;
import dr.evolution.datatype.DataType;
import dr.evomodel.sitemodel.GammaSiteModel;
import dr.evomodel.substmodel.EmpiricalAminoAcidModel;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodelxml.BinarySubstitutionModelParser;
import dr.evomodelxml.HKYParser;
import dr.inference.model.CompoundParameter;
import dr.inference.model.ParameterParser;
import dr.util.Attribute;

/**
 * @author Alexei Drummond
 */
public class PartitionModelGenerator extends Generator {

    public PartitionModelGenerator(BeautiOptions options) {
        super(options);
    }

    /**
     * Writes the substitution model to XML.
     *
     * @param writer the writer
     * @param model  the partition model to write in BEAST XML
     */
    public void writeSubstitutionModel(PartitionModel model, XMLWriter writer) {

        DataType dataType = model.dataType;
        String dataTypeDescription = dataType.getDescription();

        switch (dataType.getType()) {
            case DataType.NUCLEOTIDES:
                // Jukes-Cantor model
                if (model.getNucSubstitutionModel() == NucModelType.JC) {
                    writer.writeComment("The JC substitution model (Jukes & Cantor, 1969)");
                    writer.writeOpenTag(NucModelType.HKY.getXMLName(),
                            new Attribute[]{new Attribute.Default<String>("id", model.getName() + ".jc")}
                    );
                    writer.writeOpenTag(HKYParser.FREQUENCIES);
                    writer.writeOpenTag(
                            FrequencyModel.FREQUENCY_MODEL,
                            new Attribute[]{
                                    new Attribute.Default<String>("dataType", dataTypeDescription)
                            }
                    );
                    writer.writeOpenTag(FrequencyModel.FREQUENCIES);
                    writer.writeTag(
                            ParameterParser.PARAMETER,
                            new Attribute[]{
                                    new Attribute.Default<String>("id", model.getName() + ".jc.frequencies"),
                                    new Attribute.Default<String>("value", "0.25 0.25 0.25 0.25")
                            },
                            true
                    );
                    writer.writeCloseTag(FrequencyModel.FREQUENCIES);

                    writer.writeCloseTag(FrequencyModel.FREQUENCY_MODEL);
                    writer.writeCloseTag(HKYParser.FREQUENCIES);

                    writer.writeOpenTag(HKYParser.KAPPA);
                    writeParameter("jc.kappa", 1, 1.0, Double.NaN, Double.NaN, writer);
                    writer.writeCloseTag(HKYParser.KAPPA);
                    writer.writeCloseTag(NucModelType.HKY.getXMLName());

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
                break;

            case DataType.AMINO_ACIDS:
                // Amino Acid model
                String aaModel = "";

                switch (model.getAaSubstitutionModel()) {
                    case 0:
                        aaModel = EmpiricalAminoAcidModel.BLOSUM_62;
                        break;
                    case 1:
                        aaModel = EmpiricalAminoAcidModel.DAYHOFF;
                        break;
                    case 2:
                        aaModel = EmpiricalAminoAcidModel.JTT;
                        break;
                    case 3:
                        aaModel = EmpiricalAminoAcidModel.MT_REV_24;
                        break;
                    case 4:
                        aaModel = EmpiricalAminoAcidModel.CP_REV_45;
                        break;
                    case 5:
                        aaModel = EmpiricalAminoAcidModel.WAG;
                        break;
                }

                writer.writeComment("The " + aaModel + " substitution model");
                writer.writeTag(
                        EmpiricalAminoAcidModel.EMPIRICAL_AMINO_ACID_MODEL,
                        new Attribute[]{new Attribute.Default<String>("id", model.getName() + ".aa"),
                                new Attribute.Default<String>("type", aaModel)}, true
                );

                break;

            case DataType.TWO_STATES:
            case DataType.COVARION:

                switch (model.getBinarySubstitutionModel()) {
                    case ModelOptions.BIN_SIMPLE:
                        writeBinarySimpleModel(writer, model);
                        break;
                    case ModelOptions.BIN_COVARION:
                        writeBinaryCovarionModel(writer, model);
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
     * @param model  the partition model to write in BEAST XML
     */
    public void writeHKYModel(int num, XMLWriter writer, PartitionModel model) {
        String dataTypeDescription = model.dataType.getDescription();

        String id = model.getName() + ".hky";
        if (num > 0) {
            id += Integer.toString(num);
        }
        // Hasegawa Kishino and Yano 85 model
        writer.writeComment("The HKY substitution model (Hasegawa, Kishino & Yano, 1985)");
        writer.writeOpenTag(NucModelType.HKY.getXMLName(),
                new Attribute[]{new Attribute.Default<String>("id", id)}
        );
        writer.writeOpenTag(HKYParser.FREQUENCIES);
        writer.writeOpenTag(
                FrequencyModel.FREQUENCY_MODEL,
                new Attribute[]{
                        new Attribute.Default<String>("dataType", dataTypeDescription)
                }
        );

        writer.writeTag("patterns", new Attribute[]{new Attribute.Default<String>("idref", model.getName() + ".patterns" + (num > 0 ? num : ""))}, true);
        writer.writeOpenTag(FrequencyModel.FREQUENCIES);
        if (model.getFrequencyPolicy() == ModelOptions.ALLEQUAL)
            writeParameter(id + ".frequencies", 4, writer);
        else
            writeParameter(id + ".frequencies", 4, Double.NaN, Double.NaN, Double.NaN, writer);
        writer.writeCloseTag(FrequencyModel.FREQUENCIES);
        writer.writeCloseTag(FrequencyModel.FREQUENCY_MODEL);
        writer.writeCloseTag(HKYParser.FREQUENCIES);

        writeParameter(HKYParser.KAPPA, id + ".kappa", writer, model);
        writer.writeCloseTag(NucModelType.HKY.getXMLName());
    }

    /**
     * Write the GTR model XML block.
     *
     * @param num    the model number
     * @param writer the writer
     * @param model  the partition model to write in BEAST XML
     */
    public void writeGTRModel(int num, XMLWriter writer, PartitionModel model) {
        String dataTypeDescription = model.dataType.getDescription();

        String id = model.getName() + ".gtr";
        if (num > 0) {
            id += Integer.toString(num);
        }

        writer.writeComment("The general time reversible (GTR) substitution model");
        writer.writeOpenTag(
                dr.evomodel.substmodel.GTR.GTR_MODEL,
                new Attribute[]{new Attribute.Default<String>("id", id)}
        );
        writer.writeOpenTag(dr.evomodel.substmodel.GTR.FREQUENCIES);
        writer.writeOpenTag(
                FrequencyModel.FREQUENCY_MODEL,
                new Attribute[]{
                        new Attribute.Default<String>("dataType", dataTypeDescription)
                }
        );
        writer.writeTag("patterns", new Attribute[]{new Attribute.Default<String>("idref", model.getName() + ".patterns" + (num > 0 ? num : ""))}, true);
        writer.writeOpenTag(FrequencyModel.FREQUENCIES);
        if (model.getFrequencyPolicy() == ModelOptions.ALLEQUAL)
            writeParameter(id + ".frequencies", 4, writer);
        else
            writeParameter(id + ".frequencies", 4, Double.NaN, Double.NaN, Double.NaN, writer);
        writer.writeCloseTag(FrequencyModel.FREQUENCIES);
        writer.writeCloseTag(FrequencyModel.FREQUENCY_MODEL);
        writer.writeCloseTag(dr.evomodel.substmodel.GTR.FREQUENCIES);

        writeParameter(dr.evomodel.substmodel.GTR.A_TO_C, id + ".ac", writer, model);
        writeParameter(dr.evomodel.substmodel.GTR.A_TO_G, id + ".ag", writer, model);
        writeParameter(dr.evomodel.substmodel.GTR.A_TO_T, id + ".at", writer, model);
        writeParameter(dr.evomodel.substmodel.GTR.C_TO_G, id + ".cg", writer, model);
        writeParameter(dr.evomodel.substmodel.GTR.G_TO_T, id + ".gt", writer, model);
    }


    /**
     * Write the Binary  simple model XML block.
     *
     * @param writer the writer
     * @param model  the partition model to write in BEAST XML
     */
    public void writeBinarySimpleModel(XMLWriter writer, PartitionModel model) {
        String dataTypeDescription = model.dataType.getDescription();

        final String id = model.getName() + ".bsimple";

        writer.writeComment("The Binary simple model (based on the general substitution model)");
        writer.writeOpenTag(
                BinarySubstitutionModelParser.BINARY_SUBSTITUTION_MODEL,
                new Attribute[]{new Attribute.Default<String>("id", id)}
        );
        writer.writeOpenTag(dr.evomodel.substmodel.GeneralSubstitutionModel.FREQUENCIES);
        writer.writeOpenTag(
                FrequencyModel.FREQUENCY_MODEL,
                new Attribute[]{
                        new Attribute.Default<String>("dataType", dataTypeDescription)
                }
        );
        // TODO what alignment should I actually use here?
        writer.writeTag("alignment", new Attribute[]{new Attribute.Default<String>("idref", "alignment")}, true);
        writer.writeOpenTag(FrequencyModel.FREQUENCIES);
        writeParameter(id + ".frequencies", 2, Double.NaN, Double.NaN, Double.NaN, writer);
        writer.writeCloseTag(FrequencyModel.FREQUENCIES);
        writer.writeCloseTag(FrequencyModel.FREQUENCY_MODEL);
        writer.writeCloseTag(dr.evomodel.substmodel.GeneralSubstitutionModel.FREQUENCIES);

        writer.writeCloseTag(BinarySubstitutionModelParser.BINARY_SUBSTITUTION_MODEL);
    }

    /**
     * Write the Binary covarion model XML block
     *
     * @param writer the writer
     */
    public void writeBinaryCovarionModel(XMLWriter writer, PartitionModel model) {
        String id = model.getName() + "bcov";

        writer.writeComment("The Binary covarion model");
        writer.writeOpenTag(
                dr.evomodel.substmodel.BinaryCovarionModel.COVARION_MODEL,
                new Attribute[]{new Attribute.Default<String>("id", id)}
        );

        writeParameter(dr.evomodel.substmodel.BinaryCovarionModel.FREQUENCIES,
                id + ".frequencies", 2, 0.5, 0.0, 1.0, writer);

        writeParameter(dr.evomodel.substmodel.BinaryCovarionModel.HIDDEN_FREQUENCIES,
                id + ".hfrequencies", 2, 0.5, 0.0, 1.0, writer);

        writeParameter(dr.evomodel.substmodel.BinaryCovarionModel.ALPHA, id + ".alpha", writer, options);
        writeParameter(dr.evomodel.substmodel.BinaryCovarionModel.SWITCHING_RATE, id + ".s", writer, options);

        writer.writeCloseTag(dr.evomodel.substmodel.BinaryCovarionModel.COVARION_MODEL);
    }

    /**
     * Write the site model XML block.
     *
     * @param writer the writer
     * @param model  the partition model to write in BEAST XML
     */
    public void writeSiteModel(PartitionModel model, boolean writeMuParameter, XMLWriter writer) {

        switch (model.dataType.getType()) {
            case DataType.NUCLEOTIDES:
                if (model.getCodonHeteroPattern() != null) {
                    for (int i = 1; i <= model.getCodonPartitionCount(); i++) {
                        writeNucSiteModel(i, writeMuParameter, writer, model);
                    }
                    writer.println();
                } else {
                    writeNucSiteModel(-1, writeMuParameter, writer, model);
                }
                break;

            case DataType.AMINO_ACIDS:
                writeAASiteModel(writer, writeMuParameter, model);
                break;

            case DataType.TWO_STATES:
            case DataType.COVARION:
                writeTwoStateSiteModel(writer, writeMuParameter, model);
                break;

            default:
                throw new IllegalArgumentException("Unknown data type");
        }
    }

    /**
     * Write the all the mu parameters for this partition model.
     *
     * @param writer the writer
     * @param model  the partition model to write in BEAST XML
     */
    public void writeMuParameterRefs(PartitionModel model, XMLWriter writer) {

        if (model.dataType.getType() == DataType.NUCLEOTIDES && model.getCodonHeteroPattern() != null) {
            for (int i = 1; i <= model.getCodonPartitionCount(); i++) {
                writer.writeTag(ParameterParser.PARAMETER,
                        new Attribute[]{new Attribute.Default<String>("idref", model.getName() + ".siteModel" + i + ".mu")}, true);
            }
        } else {
            writer.writeTag(ParameterParser.PARAMETER,
                    new Attribute[]{new Attribute.Default<String>("idref", model.getName() + ".siteModel.mu")}, true);
        }

    }

    public void writeLog(XMLWriter writer, PartitionModel model) {

        int codonPartitionCount = model.getCodonPartitionCount();
        String id = model.getName() + ".";

        switch (model.dataType.getType()) {
            case DataType.NUCLEOTIDES:
                if (codonPartitionCount > 1) {
                    for (int i = 1; i <= codonPartitionCount; i++) {
                        writer.writeTag(ParameterParser.PARAMETER,
                                new Attribute.Default<String>("idref", id + "siteModel" + i + ".mu"), true);
                    }
                }
                switch (model.getNucSubstitutionModel()) {
                    case HKY:
                        if (codonPartitionCount > 1 && model.isUnlinkedSubstitutionModel()) {
                            for (int i = 1; i <= codonPartitionCount; i++) {
                                writer.writeTag(ParameterParser.PARAMETER,
                                        new Attribute.Default<String>("idref", id + "hky" + i + ".kappa"), true);
                            }
                        } else {
                            writer.writeTag(ParameterParser.PARAMETER,
                                    new Attribute.Default<String>("idref", id + "hky.kappa"), true);
                        }
                        break;

                    case GTR:
                        id = model.getName() + ".gtr";
                        if (codonPartitionCount > 1 && model.isUnlinkedSubstitutionModel()) {

                            for (int i = 1; i <= codonPartitionCount; i++) {
                                writer.writeTag(ParameterParser.PARAMETER,
                                        new Attribute.Default<String>("idref", id + i + ".ac"), true);
                                writer.writeTag(ParameterParser.PARAMETER,
                                        new Attribute.Default<String>("idref", id + i + ".ag"), true);
                                writer.writeTag(ParameterParser.PARAMETER,
                                        new Attribute.Default<String>("idref", id + i + ".at"), true);
                                writer.writeTag(ParameterParser.PARAMETER,
                                        new Attribute.Default<String>("idref", id + i + ".cg"), true);
                                writer.writeTag(ParameterParser.PARAMETER,
                                        new Attribute.Default<String>("idref", id + i + ".gt"), true);
                            }
                        } else {
                            writer.writeTag(ParameterParser.PARAMETER,
                                    new Attribute.Default<String>("idref", id + ".ac"), true);
                            writer.writeTag(ParameterParser.PARAMETER,
                                    new Attribute.Default<String>("idref", id + ".ag"), true);
                            writer.writeTag(ParameterParser.PARAMETER,
                                    new Attribute.Default<String>("idref", id + ".at"), true);
                            writer.writeTag(ParameterParser.PARAMETER,
                                    new Attribute.Default<String>("idref", id + ".cg"), true);
                            writer.writeTag(ParameterParser.PARAMETER,
                                    new Attribute.Default<String>("idref", id + ".gt"), true);
                        }
                        break;
                }
                break;//NUCLEOTIDES

            case DataType.AMINO_ACIDS:
                break;//AMINO_ACIDS

            case DataType.TWO_STATES:
            case DataType.COVARION:

                switch (model.getBinarySubstitutionModel()) {
                    case ModelOptions.BIN_SIMPLE:
                        break;
                    case ModelOptions.BIN_COVARION:
                        writer.writeTag(ParameterParser.PARAMETER,
                                new Attribute.Default<String>("idref", id + "bcov.alpha"), true);
                        writer.writeTag(ParameterParser.PARAMETER,
                                new Attribute.Default<String>("idref", id + "bcov.s"), true);
                        writer.writeTag(ParameterParser.PARAMETER,
                                new Attribute.Default<String>("idref", id + "bcov.frequencies"), true);
                        writer.writeTag(ParameterParser.PARAMETER,
                                new Attribute.Default<String>("idref", id + "bcov.hfrequencies"), true);
                        break;

                }
                break;//BINARY
        }

        if (model.isGammaHetero()) {
            if (codonPartitionCount > 1 && model.isUnlinkedHeterogeneityModel()) {
                for (int i = 1; i <= codonPartitionCount; i++) {
                    writer.writeTag(ParameterParser.PARAMETER,
                            new Attribute.Default<String>("idref", id + "siteModel" + i + ".alpha"), true);
                }
            } else {
                writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default<String>("idref", id + "siteModel.alpha"), true);
            }
        }

        if (model.isInvarHetero()) {
            if (codonPartitionCount > 1 && model.isUnlinkedHeterogeneityModel()) {
                for (int i = 1; i <= codonPartitionCount; i++) {
                    writer.writeTag(ParameterParser.PARAMETER,
                            new Attribute.Default<String>("idref", id + "siteModel" + i + ".pInv"), true);
                }
            } else {
                writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default<String>("idref", id + "siteModel.pInv"), true);
            }
        }
    }

    /**
     * Write the nucleotide site model XML block.
     *
     * @param num    the model number
     * @param writeMuParameter
     * @param writer the writer
     * @param model  the partition model to write in BEAST XML
     */
    private void writeNucSiteModel(int num, boolean writeMuParameter, XMLWriter writer, PartitionModel model) {

        String id = model.getName() + ".siteModel";
        if (num > 0) {
            id += Integer.toString(num);
        }

        writer.writeComment("site model");
        writer.writeOpenTag(GammaSiteModel.SITE_MODEL, new Attribute[]{new Attribute.Default<String>("id", id)});


        writer.writeOpenTag(GammaSiteModel.SUBSTITUTION_MODEL);

        if (model.isUnlinkedSubstitutionModel()) {
            switch (model.getNucSubstitutionModel()) {
                // JC cannot be unlinked because it has no parameters
                case JC:
                    writer.writeTag(NucModelType.HKY.getXMLName(), new Attribute.Default<String>("idref", "jc"), true);
                    break;
                case HKY:
                    writer.writeTag(NucModelType.HKY.getXMLName(),
                            new Attribute.Default<String>("idref", model.getName() + ".hky" + num), true);
                    break;
                case GTR:
                    writer.writeTag(dr.evomodel.substmodel.GTR.GTR_MODEL,
                            new Attribute.Default<String>("idref", model.getName() + ".gtr" + num), true);
                    break;
                default:
                    throw new IllegalArgumentException("Unknown substitution model.");
            }
        } else {
            switch (model.getNucSubstitutionModel()) {
                case JC:
                    writer.writeTag(NucModelType.HKY.getXMLName(),
                            new Attribute.Default<String>("idref", model.getName() + ".jc"), true);
                    break;
                case HKY:
                    writer.writeTag(NucModelType.HKY.getXMLName(),
                            new Attribute.Default<String>("idref", model.getName() + ".hky"), true);
                    break;
                case GTR:
                    writer.writeTag(dr.evomodel.substmodel.GTR.GTR_MODEL,
                            new Attribute.Default<String>("idref", model.getName() + ".gtr"), true);
                    break;
                default:
                    throw new IllegalArgumentException("Unknown substitution model.");
            }
        }
        writer.writeCloseTag(GammaSiteModel.SUBSTITUTION_MODEL);

        if (writeMuParameter) {
            writeParameter(GammaSiteModel.RELATIVE_RATE, id + ".mu", writer, model);
        }

        if (model.isGammaHetero()) {
            writer.writeOpenTag(GammaSiteModel.GAMMA_SHAPE, new Attribute.Default<String>(GammaSiteModel.GAMMA_CATEGORIES, "" + model.getGammaCategories()));
            if (num == -1 || model.isUnlinkedHeterogeneityModel()) {
                writeParameter(id + ".alpha", writer, model);
            } else {
                // multiple partitions but linked heterogeneity
                if (num == 1) {
                    writeParameter(id + ".alpha", writer, model);
                } else {
                    writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default<String>("idref", id + ".alpha"), true);
                }
            }
            writer.writeCloseTag(GammaSiteModel.GAMMA_SHAPE);
        }

        if (model.isInvarHetero()) {
            writer.writeOpenTag(GammaSiteModel.PROPORTION_INVARIANT);
            if (num == -1 || model.isUnlinkedHeterogeneityModel()) {
                writeParameter(id + ".pInv", writer, model);
            } else {
                // multiple partitions but linked heterogeneity
                if (num == 1) {
                    writeParameter(id + ".pInv", writer, model);
                } else {
                    writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default<String>("idref", id + ".pInv"), true);
                }
            }
            writer.writeCloseTag(GammaSiteModel.PROPORTION_INVARIANT);
        }

        writer.writeCloseTag(GammaSiteModel.SITE_MODEL);
    }

    /**
     * Write the two states site model XML block.
     *
     * @param writer the writer
     * @param writeMuParameter
     * @param model  the partition model to write in BEAST XML
     */
    private void writeTwoStateSiteModel(XMLWriter writer, boolean writeMuParameter, PartitionModel model) {

        String id = model.getName() + ".siteModel";

        writer.writeComment("site model");
        writer.writeOpenTag(GammaSiteModel.SITE_MODEL, new Attribute[]{new Attribute.Default<String>("id", id)});


        writer.writeOpenTag(GammaSiteModel.SUBSTITUTION_MODEL);

        switch (model.getBinarySubstitutionModel()) {
            case ModelOptions.BIN_SIMPLE:
                //writer.writeTag(dr.evomodel.substmodel.GeneralSubstitutionModel.GENERAL_SUBSTITUTION_MODEL, new Attribute.Default<String>("idref", "bsimple"), true);
                writer.writeTag(BinarySubstitutionModelParser.BINARY_SUBSTITUTION_MODEL, new Attribute.Default<String>("idref", "bsimple"), true);
                break;
            case ModelOptions.BIN_COVARION:
                writer.writeTag(dr.evomodel.substmodel.BinaryCovarionModel.COVARION_MODEL, new Attribute.Default<String>("idref", "bcov"), true);
                break;
            default:
                throw new IllegalArgumentException("Unknown substitution model.");
        }

        writer.writeCloseTag(GammaSiteModel.SUBSTITUTION_MODEL);

        if (writeMuParameter) {
            writeParameter(GammaSiteModel.RELATIVE_RATE, id + ".mu", writer, model);
        }

        if (model.isGammaHetero()) {
            writer.writeOpenTag(GammaSiteModel.GAMMA_SHAPE,
                    new Attribute.Default<String>(GammaSiteModel.GAMMA_CATEGORIES, "" + model.getGammaCategories()));
            writeParameter(id + ".alpha", writer, model);
            writer.writeCloseTag(GammaSiteModel.GAMMA_SHAPE);
        }

        if (model.isInvarHetero()) {
            writeParameter(GammaSiteModel.PROPORTION_INVARIANT, id + ".pInv", writer, model);
        }

        writer.writeCloseTag(GammaSiteModel.SITE_MODEL);
    }

    /**
     * Write the AA site model XML block.
     *
     * @param writer the writer
     * @param writeMuParameter
     * @param model  the partition model to write in BEAST XML
     */
    private void writeAASiteModel(XMLWriter writer, boolean writeMuParameter, PartitionModel model) {

        String id = model.getName() + ".siteModel";

        writer.writeComment("site model");
        writer.writeOpenTag(GammaSiteModel.SITE_MODEL, new Attribute[]{
                new Attribute.Default<String>("id", id)});


        writer.writeOpenTag(GammaSiteModel.SUBSTITUTION_MODEL);
        writer.writeTag(EmpiricalAminoAcidModel.EMPIRICAL_AMINO_ACID_MODEL,
                new Attribute.Default<String>("idref", "aa"), true);
        writer.writeCloseTag(GammaSiteModel.SUBSTITUTION_MODEL);

        if (writeMuParameter) {
            writeParameter(GammaSiteModel.RELATIVE_RATE, id + ".mu", writer, model);
        }


        if (model.isGammaHetero()) {
            writer.writeOpenTag(GammaSiteModel.GAMMA_SHAPE,
                    new Attribute.Default<String>(
                            GammaSiteModel.GAMMA_CATEGORIES, "" + model.getGammaCategories()));
            writeParameter(id + ".alpha", writer, model);
            writer.writeCloseTag(GammaSiteModel.GAMMA_SHAPE);
        }

        if (model.isInvarHetero()) {
            writeParameter(GammaSiteModel.PROPORTION_INVARIANT, id + ".pInv", writer, model);
        }

        writer.writeCloseTag(GammaSiteModel.SITE_MODEL);
    }

}