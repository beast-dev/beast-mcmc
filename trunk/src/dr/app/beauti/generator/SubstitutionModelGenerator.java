package dr.app.beauti.generator;

import dr.app.beauti.util.XMLWriter;
import dr.app.beauti.components.ComponentFactory;
import dr.app.beauti.enumTypes.FrequencyPolicyType;
import dr.app.beauti.enumTypes.NucModelType;
import dr.app.beauti.options.*;
import dr.evolution.datatype.DataType;
import dr.evolution.datatype.Nucleotides;
import dr.evomodel.sitemodel.GammaSiteModel;
import dr.evomodel.sitemodel.SiteModel;
import dr.evomodel.substmodel.BinaryCovarionModel;
import dr.evomodel.substmodel.EmpiricalAminoAcidModel;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodelxml.BinarySubstitutionModelParser;
import dr.evomodelxml.HKYParser;
import dr.evoxml.AlignmentParser;
import dr.evoxml.MergePatternsParser;
import dr.evoxml.SitePatternsParser;
import dr.inference.model.ParameterParser;
import dr.util.Attribute;
import dr.xml.AttributeParser;
import dr.xml.XMLParser;

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
    public void writeSubstitutionModel(PartitionSubstitutionModel model, XMLWriter writer) {

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
                            FrequencyModel.FREQUENCY_MODEL,
                            new Attribute[]{
                                    new Attribute.Default<String>("dataType", dataTypeDescription)
                            }
                    );
                    writer.writeOpenTag(FrequencyModel.FREQUENCIES);
                    writer.writeTag(
                            ParameterParser.PARAMETER,
                            new Attribute[]{
                                    new Attribute.Default<String>(XMLParser.ID, prefix + "frequencies"),
                                    new Attribute.Default<String>(ParameterParser.VALUE, "0.25 0.25 0.25 0.25")
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
                String aaModel = model.getAaSubstitutionModel().getXMLName();

                writer.writeComment("The " + aaModel + " substitution model");
                writer.writeTag(
                        EmpiricalAminoAcidModel.EMPIRICAL_AMINO_ACID_MODEL,
                        new Attribute[]{new Attribute.Default<String>(XMLParser.ID, model.getPrefix() + "aa"),
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
    public void writeHKYModel(int num, XMLWriter writer, PartitionSubstitutionModel model) {

        String prefix = model.getPrefix(num);

        // Hasegawa Kishino and Yano 85 model
        writer.writeComment("The HKY substitution model (Hasegawa, Kishino & Yano, 1985)");
        writer.writeOpenTag(NucModelType.HKY.getXMLName(),
                new Attribute[]{new Attribute.Default<String>(XMLParser.ID, prefix + "hky")}
        );
        writer.writeOpenTag(HKYParser.FREQUENCIES);
        writeFrequencyModel(writer, model, num);
        writer.writeCloseTag(HKYParser.FREQUENCIES);

        writeParameter(num, HKYParser.KAPPA, "kappa", model, writer);
        writer.writeCloseTag(NucModelType.HKY.getXMLName());
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
        writer.writeOpenTag(
                dr.evomodel.substmodel.GTR.GTR_MODEL,
                new Attribute[]{new Attribute.Default<String>(XMLParser.ID, prefix + "gtr")}
        );
        writer.writeOpenTag(dr.evomodel.substmodel.GTR.FREQUENCIES);
        writeFrequencyModel(writer, model, num);
        writer.writeCloseTag(dr.evomodel.substmodel.GTR.FREQUENCIES);

        writeParameter(num, dr.evomodel.substmodel.GTR.A_TO_C, PartitionSubstitutionModel.GTR_RATE_NAMES[0], model, writer);
        writeParameter(num, dr.evomodel.substmodel.GTR.A_TO_G, PartitionSubstitutionModel.GTR_RATE_NAMES[1], model, writer);
        writeParameter(num, dr.evomodel.substmodel.GTR.A_TO_T, PartitionSubstitutionModel.GTR_RATE_NAMES[2], model, writer);
        writeParameter(num, dr.evomodel.substmodel.GTR.C_TO_G, PartitionSubstitutionModel.GTR_RATE_NAMES[3], model, writer);
        writeParameter(num, dr.evomodel.substmodel.GTR.G_TO_T, PartitionSubstitutionModel.GTR_RATE_NAMES[4], model, writer);
        writer.writeCloseTag(dr.evomodel.substmodel.GTR.GTR_MODEL);
    }

    private void writeFrequencyModel(XMLWriter writer, PartitionSubstitutionModel model, int num) {
        String dataTypeDescription = model.getDataType().getDescription();
        String prefix = model.getPrefix(num);
        writer.writeOpenTag(
                FrequencyModel.FREQUENCY_MODEL,
                new Attribute[]{
                        new Attribute.Default<String>("dataType", dataTypeDescription)
                }
        );

        if (model.getFrequencyPolicy() == FrequencyPolicyType.EMPIRICAL) {
        	if (model.getDataType() == Nucleotides.INSTANCE && model.getCodonHeteroPattern() != null && model.getCodonPartitionCount() > 1) {
        		for (PartitionData partition : model.getAllPartitionData()) { //?
        			writer.writeIDref(MergePatternsParser.MERGE_PATTERNS, prefix + partition.getName() + "." + SitePatternsParser.PATTERNS);    	    			
        		}   		
        	} else { 
        		for (PartitionData partition : model.getAllPartitionData()) { //?
        			writer.writeIDref(AlignmentParser.ALIGNMENT, partition.getAlignment().getId());    
        		}
        	}
        }

        writer.writeOpenTag(FrequencyModel.FREQUENCIES);
        switch (model.getFrequencyPolicy()) {
            case ALLEQUAL:
            case ESTIMATED:
                writer.writeTag(
                        ParameterParser.PARAMETER,
                        new Attribute[]{
                                new Attribute.Default<String>(XMLParser.ID, prefix + "frequencies"),
                                new Attribute.Default<String>(ParameterParser.VALUE, "0.25 0.25 0.25 0.25")
                        }, true);
                break;
            case EMPIRICAL:
                writeParameter(prefix + "frequencies", 4, Double.NaN, Double.NaN, Double.NaN, writer);
                break;
        }
        writer.writeCloseTag(FrequencyModel.FREQUENCIES);
        writer.writeCloseTag(FrequencyModel.FREQUENCY_MODEL);
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
        writer.writeOpenTag(dr.evomodel.substmodel.GeneralSubstitutionModel.FREQUENCIES);
        writer.writeOpenTag(
                FrequencyModel.FREQUENCY_MODEL,
                new Attribute[]{
                        new Attribute.Default<String>("dataType", dataTypeDescription)
                }
        );
        
        if (model.getFrequencyPolicy() == FrequencyPolicyType.EMPIRICAL) {
        	if (model.getDataType() == Nucleotides.INSTANCE && model.getCodonHeteroPattern() != null && model.getCodonPartitionCount() > 1) {
        		for (PartitionData partition : model.getAllPartitionData()) { //?
        			writer.writeIDref(MergePatternsParser.MERGE_PATTERNS, prefix + partition.getName() + "." + SitePatternsParser.PATTERNS);    	    			
        		}   		
        	} else { 
        		for (PartitionData partition : model.getAllPartitionData()) { //?
        			writer.writeIDref(AlignmentParser.ALIGNMENT, partition.getAlignment().getId());    
        		}
        	}
        }
        
        writer.writeOpenTag(FrequencyModel.FREQUENCIES);
        writeParameter(prefix + "frequencies", 2, Double.NaN, Double.NaN, Double.NaN, writer);
        writer.writeCloseTag(FrequencyModel.FREQUENCIES);
        writer.writeCloseTag(FrequencyModel.FREQUENCY_MODEL);
        writer.writeCloseTag(dr.evomodel.substmodel.GeneralSubstitutionModel.FREQUENCIES);

        writer.writeCloseTag(BinarySubstitutionModelParser.BINARY_SUBSTITUTION_MODEL);
    }

    /**
     * Write the Binary covarion model XML block
     *
     * @param writer the writer
     * @param model  the partition model to write
     */
    public void writeBinaryCovarionModel(XMLWriter writer, PartitionSubstitutionModel model) {
        String prefix = model.getPrefix();

        writer.writeComment("The Binary covarion model");
        writer.writeOpenTag(
                dr.evomodel.substmodel.BinaryCovarionModel.COVARION_MODEL,
                new Attribute[]{new Attribute.Default<String>(XMLParser.ID, prefix + "bcov")}
        );

        writeParameter(dr.evomodel.substmodel.BinaryCovarionModel.FREQUENCIES,
                prefix + "frequencies", 2, 0.5, 0.0, 1.0, writer);

        writeParameter(dr.evomodel.substmodel.BinaryCovarionModel.HIDDEN_FREQUENCIES,
                prefix + "hfrequencies", 2, 0.5, 0.0, 1.0, writer);

        writeParameter(dr.evomodel.substmodel.BinaryCovarionModel.ALPHA, "alpha", options, writer);
        writeParameter(dr.evomodel.substmodel.BinaryCovarionModel.SWITCHING_RATE, "bcov.s", options, writer);

        writer.writeCloseTag(dr.evomodel.substmodel.BinaryCovarionModel.COVARION_MODEL);
    }

    /**
     * Write the site model XML block.
     *
     * @param model            the partition model to write in BEAST XML
     * @param writeMuParameter the relative rate parameter for this site model
     * @param writer           the writer
     */
    public void writeSiteModel(PartitionSubstitutionModel model, boolean writeMuParameter, XMLWriter writer) {

        switch (model.getDataType().getType()) {
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
    public void writeMuParameterRefs(PartitionSubstitutionModel model, XMLWriter writer) {

        if (model.getDataType().getType() == DataType.NUCLEOTIDES && model.getCodonHeteroPattern() != null) {
            for (int i = 1; i <= model.getCodonPartitionCount(); i++) {
                writer.writeIDref(ParameterParser.PARAMETER, model.getPrefix(i) + "mu");
            }
        } else {
            writer.writeIDref(ParameterParser.PARAMETER, model.getPrefix() + "mu");
        }

    }

    public void writeLog(XMLWriter writer, PartitionSubstitutionModel model) {

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
                if (codonPartitionCount > 1 && model.isUnlinkedSubstitutionModel()) {
                	for (int i = 1; i <= codonPartitionCount; i++) {
                		writer.writeIDref(ParameterParser.PARAMETER, model.getPrefix(i) + "frequencies");
                    }
                } else {
                	writer.writeIDref(ParameterParser.PARAMETER, model.getPrefix() + "frequencies");
                }
                break;//NUCLEOTIDES

            case DataType.AMINO_ACIDS:
                break;//AMINO_ACIDS

            case DataType.TWO_STATES:
            case DataType.COVARION:

                String prefix = model.getPrefix();
                switch (model.getBinarySubstitutionModel()) {
                    case ModelOptions.BIN_SIMPLE:
                        break;
                    case ModelOptions.BIN_COVARION:
                        writer.writeIDref(ParameterParser.PARAMETER, prefix + "alpha");
                        writer.writeIDref(ParameterParser.PARAMETER, prefix + "bcov.s");
                        writer.writeIDref(ParameterParser.PARAMETER, prefix + "frequencies");
                        writer.writeIDref(ParameterParser.PARAMETER, prefix + "hfrequencies");
                        break;

                }
                break;//BINARY
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

    /**
     * Write the nucleotide site model XML block.
     *
     * @param num              the model number
     * @param writeMuParameter the relative rate parameter for this site model
     * @param writer           the writer
     * @param model            the partition model to write in BEAST XML
     */
    private void writeNucSiteModel(int num, boolean writeMuParameter, XMLWriter writer, PartitionSubstitutionModel model) {

        String prefix = model.getPrefix(num);
        String prefix2 = model.getPrefix();

        writer.writeComment("site model");
        writer.writeOpenTag(GammaSiteModel.SITE_MODEL,
                new Attribute[]{new Attribute.Default<String>(XMLParser.ID, prefix + SiteModel.SITE_MODEL)});


        writer.writeOpenTag(GammaSiteModel.SUBSTITUTION_MODEL);

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
                    writer.writeIDref(dr.evomodel.substmodel.GTR.GTR_MODEL, prefix + "gtr");
                    break;
                default:
                    throw new IllegalArgumentException("Unknown substitution model.");
            }
            

            if (writeMuParameter) {
                writeParameter(num, GammaSiteModel.RELATIVE_RATE, "mu", model, writer);
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
                    writer.writeIDref(dr.evomodel.substmodel.GTR.GTR_MODEL, prefix2 + "gtr");
                    break;
                default:
                    throw new IllegalArgumentException("Unknown substitution model.");
            }
            

            if (writeMuParameter) {
                writeParameter(GammaSiteModel.RELATIVE_RATE, "mu", model, writer);
            }
        }
        
        writer.writeCloseTag(GammaSiteModel.SUBSTITUTION_MODEL);

        if (model.isGammaHetero()) {
            writer.writeOpenTag(GammaSiteModel.GAMMA_SHAPE, new Attribute.Default<String>(GammaSiteModel.GAMMA_CATEGORIES, "" + model.getGammaCategories()));
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
            writer.writeCloseTag(GammaSiteModel.GAMMA_SHAPE);
        }

        if (model.isInvarHetero()) {
            writer.writeOpenTag(GammaSiteModel.PROPORTION_INVARIANT);
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
            writer.writeCloseTag(GammaSiteModel.PROPORTION_INVARIANT);
        }

        writer.writeCloseTag(GammaSiteModel.SITE_MODEL);
    }

    /**
     * Write the two states site model XML block.
     *
     * @param writer           the writer
     * @param writeMuParameter the relative rate parameter for this site model
     * @param model            the partition model to write in BEAST XML
     */
    private void writeTwoStateSiteModel(XMLWriter writer, boolean writeMuParameter, PartitionSubstitutionModel model) {

        String prefix = model.getPrefix();

        writer.writeComment("site model");
        writer.writeOpenTag(GammaSiteModel.SITE_MODEL,
                new Attribute[]{new Attribute.Default<String>(XMLParser.ID, prefix + SiteModel.SITE_MODEL)});


        writer.writeOpenTag(GammaSiteModel.SUBSTITUTION_MODEL);

        switch (model.getBinarySubstitutionModel()) {
            case ModelOptions.BIN_SIMPLE:
                //writer.writeIDref(dr.evomodel.substmodel.GeneralSubstitutionModel.GENERAL_SUBSTITUTION_MODEL, "bsimple");
                writer.writeIDref(BinarySubstitutionModelParser.BINARY_SUBSTITUTION_MODEL, "bsimple");
                break;
            case ModelOptions.BIN_COVARION:
                writer.writeIDref(BinaryCovarionModel.COVARION_MODEL, "bcov");
                break;
            default:
                throw new IllegalArgumentException("Unknown substitution model.");
        }

        writer.writeCloseTag(GammaSiteModel.SUBSTITUTION_MODEL);

        if (writeMuParameter) {
            writeParameter(GammaSiteModel.RELATIVE_RATE, "mu", model, writer);
        }

        if (model.isGammaHetero()) {
            writer.writeOpenTag(GammaSiteModel.GAMMA_SHAPE,
                    new Attribute.Default<String>(GammaSiteModel.GAMMA_CATEGORIES, "" + model.getGammaCategories()));
            writeParameter(prefix + "alpha", model, writer);
            writer.writeCloseTag(GammaSiteModel.GAMMA_SHAPE);
        }

        if (model.isInvarHetero()) {
            writeParameter(GammaSiteModel.PROPORTION_INVARIANT, "pInv", model, writer);
        }

        writer.writeCloseTag(GammaSiteModel.SITE_MODEL);
    }

    /**
     * Write the AA site model XML block.
     *
     * @param writer           the writer
     * @param writeMuParameter the relative rate parameter for this site model
     * @param model            the partition model to write in BEAST XML
     */
    private void writeAASiteModel(XMLWriter writer, boolean writeMuParameter, PartitionSubstitutionModel model) {

        String prefix = model.getPrefix();

        writer.writeComment("site model");
        writer.writeOpenTag(GammaSiteModel.SITE_MODEL, new Attribute[]{
                new Attribute.Default<String>(XMLParser.ID, prefix + SiteModel.SITE_MODEL)});


        writer.writeOpenTag(GammaSiteModel.SUBSTITUTION_MODEL);
        writer.writeIDref(EmpiricalAminoAcidModel.EMPIRICAL_AMINO_ACID_MODEL, prefix + "aa");
        writer.writeCloseTag(GammaSiteModel.SUBSTITUTION_MODEL);

        if (writeMuParameter) {
            writeParameter(GammaSiteModel.RELATIVE_RATE, "mu", model, writer);
        }


        if (model.isGammaHetero()) {
            writer.writeOpenTag(GammaSiteModel.GAMMA_SHAPE,
                    new Attribute.Default<String>(
                            GammaSiteModel.GAMMA_CATEGORIES, "" + model.getGammaCategories()));
            writeParameter("alpha", model, writer);
            writer.writeCloseTag(GammaSiteModel.GAMMA_SHAPE);
        }

        if (model.isInvarHetero()) {
            writeParameter(GammaSiteModel.PROPORTION_INVARIANT, "pInv", model, writer);
        }

        writer.writeCloseTag(GammaSiteModel.SITE_MODEL);
    }
}