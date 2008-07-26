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
                if (model.nucSubstitutionModel == NucModelType.JC) {
                    writer.writeComment("The JC substitution model (Jukes & Cantor, 1969)");
                    writer.writeOpenTag(NucModelType.HKY.getXMLName(),
                            new Attribute[]{new Attribute.Default<String>("id", "jc")}
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
                                    new Attribute.Default<String>("id", "jc.frequencies"),
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
                    if (model.nucSubstitutionModel == NucModelType.HKY) {
                        if (model.unlinkedSubstitutionModel) {
                            for (int i = 1; i <= model.codonPartitionCount; i++) {
                                writeHKYModel(i, writer, model);
                            }
                        } else {
                            writeHKYModel(-1, writer, model);
                        }
                    } else {
                        // General time reversible model
                        if (model.nucSubstitutionModel == NucModelType.GTR) {
                            if (model.unlinkedSubstitutionModel) {
                                for (int i = 1; i <= model.codonPartitionCount; i++) {
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

                switch (model.aaSubstitutionModel) {
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
                        new Attribute[]{new Attribute.Default<String>("id", "aa"),
                                new Attribute.Default<String>("type", aaModel)}, true
                );

                break;

            case DataType.TWO_STATES:
            case DataType.COVARION:

                switch (model.binarySubstitutionModel) {
                    case ModelOptions.BIN_SIMPLE:
                        writeBinarySimpleModel(writer, model);
                        break;
                    case ModelOptions.BIN_COVARION:
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
     * @param model  the partition model to write in BEAST XML
     */
    public void writeHKYModel(int num, XMLWriter writer, PartitionModel model) {
        String dataTypeDescription = model.dataType.getDescription();

        String id = "hky";
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
        writer.writeTag("alignment", new Attribute[]{new Attribute.Default<String>("idref", "alignment")}, true);
        writer.writeOpenTag(FrequencyModel.FREQUENCIES);
        if (model.frequencyPolicy == ModelOptions.ALLEQUAL)
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

        String id = "gtr";
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
        writer.writeTag("alignment", new Attribute[]{new Attribute.Default<String>("idref", "alignment")}, true);
        writer.writeOpenTag(FrequencyModel.FREQUENCIES);
        if (model.frequencyPolicy == ModelOptions.ALLEQUAL)
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

        final String id = "bsimple";

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
    public void writeBinaryCovarionModel(XMLWriter writer) {
        String id = "bcov";

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
    public void writeSiteModel(PartitionModel model, XMLWriter writer) {

        switch (model.dataType.getType()) {
            case DataType.NUCLEOTIDES:
                if (model.codonHeteroPattern != null) {
                    for (int i = 1; i <= model.codonPartitionCount; i++) {
                        writeNucSiteModel(i, writer, model);
                    }
                    writer.println();
                    writer.writeOpenTag(CompoundParameter.COMPOUND_PARAMETER, new Attribute[]{new Attribute.Default<String>("id", "allMus")});
                    for (int i = 1; i <= model.codonPartitionCount; i++) {
                        writer.writeTag(ParameterParser.PARAMETER,
                                new Attribute[]{new Attribute.Default<String>("idref", "siteModel" + i + ".mu")}, true);
                    }
                    writer.writeCloseTag(CompoundParameter.COMPOUND_PARAMETER);
                } else {
                    writeNucSiteModel(-1, writer, model);
                }
                break;

            case DataType.AMINO_ACIDS:
                writeAASiteModel(writer, model);
                break;

            case DataType.TWO_STATES:
            case DataType.COVARION:
                writeTwoStateSiteModel(writer, model);
                break;

            default:
                throw new IllegalArgumentException("Unknown data type");
        }
    }

    public void writeLog(XMLWriter writer, PartitionModel model) {

        int codonPartitionCount = model.codonPartitionCount;

        switch (model.dataType.getType()) {
            case DataType.NUCLEOTIDES:
                if (codonPartitionCount > 1) {
                    for (int i = 1; i <= codonPartitionCount; i++) {
                        writer.writeTag(ParameterParser.PARAMETER,
                                new Attribute.Default<String>("idref", "siteModel" + i + ".mu"), true);
                    }
                }
                switch (model.nucSubstitutionModel) {
                    case HKY:
                        if (codonPartitionCount > 1 && model.unlinkedSubstitutionModel) {
                            for (int i = 1; i <= codonPartitionCount; i++) {
                                writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default<String>("idref", "hky" + i + ".kappa"), true);
                            }
                        } else {
                            writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default<String>("idref", "hky.kappa"), true);
                        }
                        break;

                    case GTR:
                        if (codonPartitionCount > 1 && model.unlinkedSubstitutionModel) {
                            for (int i = 1; i <= codonPartitionCount; i++) {
                                writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default<String>("idref", "gtr" + i + ".ac"), true);
                                writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default<String>("idref", "gtr" + i + ".ag"), true);
                                writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default<String>("idref", "gtr" + i + ".at"), true);
                                writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default<String>("idref", "gtr" + i + ".cg"), true);
                                writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default<String>("idref", "gtr" + i + ".gt"), true);
                            }
                        } else {
                            writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default<String>("idref", "gtr.ac"), true);
                            writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default<String>("idref", "gtr.ag"), true);
                            writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default<String>("idref", "gtr.at"), true);
                            writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default<String>("idref", "gtr.cg"), true);
                            writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default<String>("idref", "gtr.gt"), true);
                        }
                        break;
                }
                break;//NUCLEOTIDES

            case DataType.AMINO_ACIDS:
                break;//AMINO_ACIDS

            case DataType.TWO_STATES:
            case DataType.COVARION:

                switch (model.binarySubstitutionModel) {
                    case ModelOptions.BIN_SIMPLE:
                        break;
                    case ModelOptions.BIN_COVARION:
                        writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default<String>("idref", "bcov.alpha"), true);
                        writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default<String>("idref", "bcov.s"), true);
                        writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default<String>("idref", "bcov.frequencies"), true);
                        writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default<String>("idref", "bcov.hfrequencies"), true);
                        break;

                }
                break;//BINARY
        }

        if (model.gammaHetero) {
            if (codonPartitionCount > 1 && model.unlinkedHeterogeneityModel) {
                for (int i = 1; i <= codonPartitionCount; i++) {
                    writer.writeTag(ParameterParser.PARAMETER,
                            new Attribute.Default<String>("idref", "siteModel" + i + ".alpha"), true);
                }
            } else {
                writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default<String>("idref", "siteModel.alpha"), true);
            }
        }

        if (model.invarHetero) {
            if (codonPartitionCount > 1 && model.unlinkedHeterogeneityModel) {
                for (int i = 1; i <= codonPartitionCount; i++) {
                    writer.writeTag(ParameterParser.PARAMETER,
                            new Attribute.Default<String>("idref", "siteModel" + i + ".pInv"), true);
                }
            } else {
                writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default<String>("idref", "siteModel.pInv"), true);
            }
        }
    }

    /**
     * Write the nucleotide site model XML block.
     *
     * @param num    the model number
     * @param writer the writer
     * @param model  the partition model to write in BEAST XML
     */
    private void writeNucSiteModel(int num, XMLWriter writer, PartitionModel model) {

        String id = "siteModel";
        if (num > 0) {
            id += Integer.toString(num);
        }

        writer.writeComment("site model");
        writer.writeOpenTag(GammaSiteModel.SITE_MODEL, new Attribute[]{new Attribute.Default<String>("id", id)});


        writer.writeOpenTag(GammaSiteModel.SUBSTITUTION_MODEL);

        if (model.unlinkedSubstitutionModel) {
            switch (model.nucSubstitutionModel) {
                // JC cannot be unlinked because it has no parameters
                case JC:
                    writer.writeTag(NucModelType.HKY.getXMLName(), new Attribute.Default<String>("idref", "jc"), true);
                    break;
                case HKY:
                    writer.writeTag(NucModelType.HKY.getXMLName(),
                            new Attribute.Default<String>("idref", "hky" + num), true);
                    break;
                case GTR:
                    writer.writeTag(dr.evomodel.substmodel.GTR.GTR_MODEL,
                            new Attribute.Default<String>("idref", "gtr" + num), true);
                    break;
                default:
                    throw new IllegalArgumentException("Unknown substitution model.");
            }
        } else {
            switch (model.nucSubstitutionModel) {
                case JC:
                    writer.writeTag(NucModelType.HKY.getXMLName(), new Attribute.Default<String>("idref", "jc"), true);
                    break;
                case HKY:
                    writer.writeTag(NucModelType.HKY.getXMLName(), new Attribute.Default<String>("idref", "hky"), true);
                    break;
                case GTR:
                    writer.writeTag(dr.evomodel.substmodel.GTR.GTR_MODEL, new Attribute.Default<String>("idref", "gtr"), true);
                    break;
                default:
                    throw new IllegalArgumentException("Unknown substitution model.");
            }
        }
        writer.writeCloseTag(GammaSiteModel.SUBSTITUTION_MODEL);

        if (num != -1) {
            writeParameter(GammaSiteModel.RELATIVE_RATE, id + ".mu", writer, model);
        } else {
//            The actual mutation rate is now in the BranchRateModel so relativeRate can be missing
        }

        if (model.gammaHetero) {
            writer.writeOpenTag(GammaSiteModel.GAMMA_SHAPE, new Attribute.Default<String>(GammaSiteModel.GAMMA_CATEGORIES, "" + model.gammaCategories));
            if (num == -1 || model.unlinkedHeterogeneityModel) {
                writeParameter(id + ".alpha", writer, model);
            } else {
                // multiple partitions but linked heterogeneity
                if (num == 1) {
                    writeParameter("siteModel.alpha", writer, model);
                } else {
                    writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default<String>("idref", "siteModel.alpha"), true);
                }
            }
            writer.writeCloseTag(GammaSiteModel.GAMMA_SHAPE);
        }

        if (model.invarHetero) {
            writer.writeOpenTag(GammaSiteModel.PROPORTION_INVARIANT);
            if (num == -1 || model.unlinkedHeterogeneityModel) {
                writeParameter(id + ".pInv", writer, model);
            } else {
                // multiple partitions but linked heterogeneity
                if (num == 1) {
                    writeParameter("siteModel.pInv", writer, model);
                } else {
                    writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default<String>("idref", "siteModel.pInv"), true);
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
     * @param model  the partition model to write in BEAST XML
     */
    private void writeTwoStateSiteModel(XMLWriter writer, PartitionModel model) {

        String id = "siteModel";

        writer.writeComment("site model");
        writer.writeOpenTag(GammaSiteModel.SITE_MODEL, new Attribute[]{new Attribute.Default<String>("id", id)});


        writer.writeOpenTag(GammaSiteModel.SUBSTITUTION_MODEL);

        switch (model.binarySubstitutionModel) {
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

        if (model.gammaHetero) {
            writer.writeOpenTag(GammaSiteModel.GAMMA_SHAPE,
                    new Attribute.Default<String>(GammaSiteModel.GAMMA_CATEGORIES, "" + model.gammaCategories));
            writeParameter(id + ".alpha", writer, model);
            writer.writeCloseTag(GammaSiteModel.GAMMA_SHAPE);
        }

        if (model.invarHetero) {
            writeParameter(GammaSiteModel.PROPORTION_INVARIANT, id + ".pInv", writer, model);
        }

        writer.writeCloseTag(GammaSiteModel.SITE_MODEL);
    }

    /**
     * Write the AA site model XML block.
     *
     * @param writer the writer
     * @param model  the partition model to write in BEAST XML
     */
    private void writeAASiteModel(XMLWriter writer, PartitionModel model) {

        writer.writeComment("site model");
        writer.writeOpenTag(GammaSiteModel.SITE_MODEL, new Attribute[]{
                new Attribute.Default<String>("id", "siteModel")});


        writer.writeOpenTag(GammaSiteModel.SUBSTITUTION_MODEL);
        writer.writeTag(EmpiricalAminoAcidModel.EMPIRICAL_AMINO_ACID_MODEL,
                new Attribute.Default<String>("idref", "aa"), true);
        writer.writeCloseTag(GammaSiteModel.SUBSTITUTION_MODEL);

        //The actual mutation rate is now in the BranchRateModel so relativeRate can be missing

        if (model.gammaHetero) {
            writer.writeOpenTag(GammaSiteModel.GAMMA_SHAPE,
                    new Attribute.Default<String>(
                            GammaSiteModel.GAMMA_CATEGORIES, "" + model.gammaCategories));
            writeParameter("siteModel.alpha", writer, model);
            writer.writeCloseTag(GammaSiteModel.GAMMA_SHAPE);
        }

        if (model.invarHetero) {
            writeParameter(GammaSiteModel.PROPORTION_INVARIANT, "siteModel.pInv", writer, model);
        }

        writer.writeCloseTag(GammaSiteModel.SITE_MODEL);
    }
}