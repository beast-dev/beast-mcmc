package dr.app.beauti.components.discrete;

import dr.app.beauti.options.*;
import dr.app.beauti.util.XMLWriter;
import dr.app.beauti.generator.BaseComponentGenerator;
import dr.evolution.datatype.GeneralDataType;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.sitemodel.SiteModel;
import dr.evomodel.substmodel.AbstractSubstitutionModel;
import dr.evomodel.tree.TreeModel;
import dr.evomodelxml.branchratemodel.*;
import dr.evomodelxml.clock.ACLikelihoodParser;
import dr.evomodelxml.sitemodel.GammaSiteModelParser;
import dr.evomodelxml.substmodel.*;
import dr.evomodelxml.treelikelihood.AncestralStateTreeLikelihoodParser;
import dr.evomodelxml.treelikelihood.TreeLikelihoodParser;
import dr.evoxml.*;
import dr.inference.model.ParameterParser;
import dr.inferencexml.loggers.ColumnsParser;
import dr.inferencexml.loggers.LoggerParser;
import dr.inferencexml.model.ProductStatisticParser;
import dr.inferencexml.model.SumStatisticParser;
import dr.util.Attribute;
import dr.xml.XMLParser;

import java.util.Set;

/**
 * @author Andrew Rambaut
 * @version $Id$
 */
public class DiscreteTraitsComponentGenerator extends BaseComponentGenerator {

    public DiscreteTraitsComponentGenerator(final BeautiOptions options) {
        super(options);
    }

    public boolean usesInsertionPoint(final InsertionPoint point) {
        if (options.getAllPartitionData(GeneralDataType.INSTANCE).size() == 0) {
            // Empty, so do nothing
            return false;
        }


        switch (point) {
            case AFTER_PATTERNS:
            case AFTER_TREE_LIKELIHOOD:
            case IN_MCMC_LIKELIHOOD:
            case IN_SCREEN_LOG:
            case IN_FILE_LOG_PARAMETERS:
            case IN_FILE_LOG_LIKELIHOODS:
            case AFTER_FILE_LOG:
            case IN_TREES_LOG:
                return true;
            default:
                return false;
        }
    }

    protected void generate(final InsertionPoint point, final Object item, final XMLWriter writer) {
        DiscreteTraitsComponentOptions comp = (DiscreteTraitsComponentOptions)options.getComponentOptions(DiscreteTraitsComponentOptions.class);

        switch (point) {
            case AFTER_PATTERNS:
                writeDiscreteTraitPatterns(writer, comp);
                break;

            case AFTER_TREE_LIKELIHOOD:
                writeAncestralTreeLikelihoods(writer, comp);
                break;

            case IN_MCMC_LIKELIHOOD:
                writeAncestralTreeLikelihoodReferences(writer);
                break;

            case IN_SCREEN_LOG:
                writeScreenLogEntries(writer);
                break;

            case IN_FILE_LOG_PARAMETERS:
                writeFileLogEntries(writer);
                break;

            case IN_FILE_LOG_LIKELIHOODS:
                writeAncestralTreeLikelihoodReferences(writer);
                break;

            case AFTER_FILE_LOG:
                writeDiscreteTraitFileLoggers(writer);

            case IN_TREES_LOG:
                writeTreeLogEntries((PartitionTreeModel)item, writer);
                break;
            default:
                throw new IllegalArgumentException("This insertion point is not implemented for " + this.getClass().getName());
        }

    }

    protected String getCommentLabel() {
        return "Discrete Traits Model";
    }

    private void writeDiscreteTraitPatterns(XMLWriter writer,
                                            DiscreteTraitsComponentOptions component) {

        for (PartitionSubstitutionModel model : options.getPartitionSubstitutionModels(GeneralDataType.INSTANCE)) {
            writeGeneralDataType(model, writer);
            writer.writeBlankLine();
        }

        // now create an attribute pattern for each trait that uses it
        for (AbstractPartitionData partition : options.getAllPartitionData(GeneralDataType.INSTANCE)) {
            if (partition.getTraits() != null) {
                writeAttributePatterns(partition, writer);
                writer.writeText("");

            }
        }

    }

    private void writeGeneralDataType(PartitionSubstitutionModel model, XMLWriter writer) {

        writer.writeComment("general data type for discrete trait model, '" + model.getName() + "'");

        Set<String> states = options.getStatesForDiscreteModel(model);

        // <generalDataType>
        writer.writeOpenTag(GeneralDataTypeParser.GENERAL_DATA_TYPE, new Attribute[]{
                new Attribute.Default<String>(XMLParser.ID, model.getPrefix() + "dataType")});

        int numOfStates = states.size();
        writer.writeComment("Number Of States = " + numOfStates);

        for (String eachGD : states) {
            writer.writeTag(GeneralDataTypeParser.STATE, new Attribute[]{
                    new Attribute.Default<String>(GeneralDataTypeParser.CODE, eachGD)}, true);
        }

        writer.writeCloseTag(GeneralDataTypeParser.GENERAL_DATA_TYPE);
    }

    /**
     * write <attributePatterns>
     *
     * @param partition PartitionData
     * @param writer    XMLWriter
     */
    private void writeAttributePatterns(AbstractPartitionData partition, XMLWriter writer) {
        writer.writeComment("Data pattern for discrete trait, '" + partition.getTraits().get(0).getName() + "'");

        // <attributePatterns>
        writer.writeOpenTag(AttributePatternsParser.ATTRIBUTE_PATTERNS, new Attribute[]{
                new Attribute.Default<String>(XMLParser.ID, partition.getPrefix() + AttributePatternsParser.ATTRIBUTE_PATTERNS),
                new Attribute.Default<String>(AttributePatternsParser.ATTRIBUTE, partition.getTraits().get(0).getName())});
        writer.writeIDref(TaxaParser.TAXA, TaxaParser.TAXA);
        writer.writeIDref(GeneralDataTypeParser.GENERAL_DATA_TYPE, partition.getPartitionSubstitutionModel().getPrefix() + "dataType");
        writer.writeCloseTag(AttributePatternsParser.ATTRIBUTE_PATTERNS);
    }

    private void writeDiscreteTraitsModels(XMLWriter writer, DiscreteTraitsComponentOptions componentOptions) {
        for (PartitionSubstitutionModel model : options.getPartitionSubstitutionModels(GeneralDataType.INSTANCE)) {
            writeDiscreteTraitsSubstitutionModel(model, writer);
            writeDiscreteTraitsSiteModel(model, writer);
        }
    }

    private void writeDiscreteTraitsSubstitutionModel(PartitionSubstitutionModel model, XMLWriter writer) {

        int stateCount = options.getStatesForDiscreteModel(model).size();

        if (model.getDiscreteSubstType() == DiscreteSubstModelType.SYM_SUBST) {
            writer.writeComment("symmetric CTMC model for discrete state reconstructions");

            writer.writeOpenTag(GeneralSubstitutionModelParser.GENERAL_SUBSTITUTION_MODEL, new Attribute[]{
                    new Attribute.Default<String>(XMLParser.ID, model.getPrefix() + AbstractSubstitutionModel.MODEL)});

            writer.writeIDref(GeneralDataTypeParser.GENERAL_DATA_TYPE, model.getPrefix() +  "dataType");

            writer.writeOpenTag(GeneralSubstitutionModelParser.FREQUENCIES);

            writeDiscreteFrequencyModel(model, stateCount, true, writer);

            writer.writeCloseTag(GeneralSubstitutionModelParser.FREQUENCIES);

            //---------------- rates and indicators -----------------

            writeRatesAndIndicators(model, stateCount * (stateCount - 1) / 2, null, writer);
            writer.writeCloseTag(GeneralSubstitutionModelParser.GENERAL_SUBSTITUTION_MODEL);

        } else if (model.getDiscreteSubstType() == DiscreteSubstModelType.ASYM_SUBST) {
            writer.writeComment("asymmetric CTMC model for discrete state reconstructions");

            writer.writeOpenTag(GeneralSubstitutionModelParser.GENERAL_SUBSTITUTION_MODEL, new Attribute[]{
                    new Attribute.Default<String>(XMLParser.ID, model.getPrefix() + AbstractSubstitutionModel.MODEL),
                    new Attribute.Default<Boolean>(ComplexSubstitutionModelParser.RANDOMIZE, false)});

            writer.writeIDref(GeneralDataTypeParser.GENERAL_DATA_TYPE, model.getPrefix() + "dataType");

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

        writer.writeIDref(GeneralDataTypeParser.GENERAL_DATA_TYPE, model.getPrefix() + "dataType");

        writer.writeOpenTag(FrequencyModelParser.FREQUENCIES);
        writeParameter(model.getPrefix() + "trait.frequencies", stateCount, Double.NaN, Double.NaN, Double.NaN, writer);
        writer.writeCloseTag(FrequencyModelParser.FREQUENCIES);

        writer.writeCloseTag(FrequencyModelParser.FREQUENCY_MODEL);
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

    private void writeRatesAndIndicators(PartitionSubstitutionModel model, int dimension, Integer relativeTo, XMLWriter writer) {
        writer.writeComment("rates and indicators");

        String prefix = model.getName() + ".";

        if (relativeTo == null) {
            writer.writeOpenTag(GeneralSubstitutionModelParser.RATES);
        } else {
            writer.writeOpenTag(GeneralSubstitutionModelParser.RATES, new Attribute[]{
                    new Attribute.Default<Integer>(GeneralSubstitutionModelParser.RELATIVE_TO, relativeTo)});
        }

        model.getParameter(prefix + "rates").isFixed = true;
        writeParameter(model.getParameter(prefix + "rates"), dimension, writer);

        writer.writeCloseTag(GeneralSubstitutionModelParser.RATES);

        if (model.isActivateBSSVS()) { //If "BSSVS" is not activated, rateIndicator should not be there.
            writer.writeOpenTag(GeneralSubstitutionModelParser.INDICATOR);
            model.getParameter(prefix + "indicators").isFixed = true;
            writeParameter(model.getParameter(prefix + "indicators"), dimension, writer);
            writer.writeCloseTag(GeneralSubstitutionModelParser.INDICATOR);
        }

    }

    private void writeStatisticModel(PartitionSubstitutionModel model, XMLWriter writer) {
        String prefix = model.getName() + ".";

        writer.writeOpenTag(SumStatisticParser.SUM_STATISTIC, new Attribute[]{
                new Attribute.Default<String>(XMLParser.ID, model.getPrefix() + prefix + "nonZeroRates"),
                new Attribute.Default<Boolean>(SumStatisticParser.ELEMENTWISE, true)});
        writer.writeIDref(ParameterParser.PARAMETER, model.getPrefix() + prefix + "indicators");
        writer.writeCloseTag(SumStatisticParser.SUM_STATISTIC);

        writer.writeOpenTag(ProductStatisticParser.PRODUCT_STATISTIC, new Attribute[]{
                new Attribute.Default<String>(XMLParser.ID, prefix + "actualRates"),
                new Attribute.Default<Boolean>(SumStatisticParser.ELEMENTWISE, false)});
        writer.writeIDref(ParameterParser.PARAMETER, prefix + "indicators");
        writer.writeIDref(ParameterParser.PARAMETER, prefix + "rates");
        writer.writeCloseTag(ProductStatisticParser.PRODUCT_STATISTIC);
    }

    private void writeAncestralTreeLikelihoods(XMLWriter writer,
                                               DiscreteTraitsComponentOptions component) {
        // generate tree likelihoods for discrete trait partitions
        if (options.hasDiscreteTraitPartition()) {
            writer.writeComment("Likelihood for tree given discrete trait data");
        }
        for (AbstractPartitionData partition : options.dataPartitions) {
            if (partition.getTraits() != null) {
                TraitData trait = partition.getTraits().get(0);
                if (trait.getTraitType() == TraitData.TraitType.DISCRETE) {
                    writeAncestralTreeLikelihood(partition, writer);
                }
            }
        }
    }

    /**
     * Ancestral Tree Likelihood
     *
     * @param partition PartitionData
     * @param writer    XMLWriter
     */
    public void writeAncestralTreeLikelihood(AbstractPartitionData partition, XMLWriter writer) {
        String prefix = partition.getName() + ".";

        PartitionSubstitutionModel substModel = partition.getPartitionSubstitutionModel();
        PartitionTreeModel treeModel = partition.getPartitionTreeModel();
        PartitionClockModel clockModel = partition.getPartitionClockModel();

        writer.writeOpenTag(AncestralStateTreeLikelihoodParser.RECONSTRUCTING_TREE_LIKELIHOOD, new Attribute[]{
                new Attribute.Default<String>(XMLParser.ID, prefix + TreeLikelihoodParser.TREE_LIKELIHOOD),
                new Attribute.Default<String>(AncestralStateTreeLikelihoodParser.TAG_NAME, prefix + AncestralStateTreeLikelihoodParser.RECONSTRUCTION_TAG),
        });

        writer.writeIDref(AttributePatternsParser.ATTRIBUTE_PATTERNS, prefix + AttributePatternsParser.ATTRIBUTE_PATTERNS);
        writer.writeIDref(TreeModel.TREE_MODEL, treeModel.getPrefix() + TreeModel.TREE_MODEL);
        writer.writeIDref(SiteModel.SITE_MODEL, substModel.getName() + "." + SiteModel.SITE_MODEL);
        writer.writeIDref(GeneralSubstitutionModelParser.GENERAL_SUBSTITUTION_MODEL, substModel.getName() + "." + AbstractSubstitutionModel.MODEL);

        switch (clockModel.getClockType()) {
            case STRICT_CLOCK:
                writer.writeIDref(StrictClockBranchRatesParser.STRICT_CLOCK_BRANCH_RATES,
                        clockModel.getPrefix() + BranchRateModel.BRANCH_RATES);
                break;
            case UNCORRELATED:
                writer.writeIDref(DiscretizedBranchRatesParser.DISCRETIZED_BRANCH_RATES,
                        clockModel.getPrefix() + BranchRateModel.BRANCH_RATES);
                break;
            case RANDOM_LOCAL_CLOCK:
                writer.writeIDref(RandomLocalClockModelParser.LOCAL_BRANCH_RATES,
                        clockModel.getPrefix() + BranchRateModel.BRANCH_RATES);
                break;
            case AUTOCORRELATED:
                writer.writeIDref(ACLikelihoodParser.AC_LIKELIHOOD,
                        clockModel.getPrefix() + BranchRateModel.BRANCH_RATES);
                break;

            default:
                throw new IllegalArgumentException("Unknown clock model");
        }

        writer.writeCloseTag(AncestralStateTreeLikelihoodParser.RECONSTRUCTING_TREE_LIKELIHOOD);
    }

    private void writeAncestralTreeLikelihoodReferences(XMLWriter writer) {
        for (AbstractPartitionData partition : options.dataPartitions) {
            if (partition.getTraits() != null) {
                TraitData trait = partition.getTraits().get(0);

                if (trait.getTraitType() == TraitData.TraitType.DISCRETE) {
                    writer.writeIDref(AncestralStateTreeLikelihoodParser.RECONSTRUCTING_TREE_LIKELIHOOD,
                            partition.getPrefix() + TreeLikelihoodParser.TREE_LIKELIHOOD);
                }
            }
        }
    }

    public void writeScreenLogEntries(XMLWriter writer) {
        for (PartitionSubstitutionModel model : options.getPartitionSubstitutionModels(GeneralDataType.INSTANCE)) {
            String prefix = model.getName() + ".";

            if (model.isActivateBSSVS()) { //If "BSSVS" is not activated, rateIndicator should not be there.
                writer.writeOpenTag(ColumnsParser.COLUMN,
                        new Attribute[]{
                                new Attribute.Default<String>(ColumnsParser.LABEL, model.getPrefix() + "nonZeroRates"),
                                new Attribute.Default<String>(ColumnsParser.SIGNIFICANT_FIGURES, "6"),
                                new Attribute.Default<String>(ColumnsParser.WIDTH, "12")
                        }
                );

                writer.writeIDref(SumStatisticParser.SUM_STATISTIC, prefix + "nonZeroRates");

                writer.writeCloseTag(ColumnsParser.COLUMN);
            }
        }
    }

    private void writeFileLogEntries(XMLWriter writer) {
        for (PartitionSubstitutionModel model : options.getPartitionSubstitutionModels(GeneralDataType.INSTANCE)) {
            String prefix = model.getName() + ".";

            writer.writeIDref(ParameterParser.PARAMETER, prefix + "rates");

            if (model.isActivateBSSVS()) { //If "BSSVS" is not activated, rateIndicator should not be there.
                writer.writeIDref(ParameterParser.PARAMETER, prefix + "indicators");
                writer.writeIDref(SumStatisticParser.SUM_STATISTIC, prefix + "nonZeroRates");
            }
        }

    }

    private void writeDiscreteTraitFileLoggers(XMLWriter writer) {
        for (PartitionSubstitutionModel model : options.getPartitionSubstitutionModels(GeneralDataType.INSTANCE)) {
            writeDiscreteTraitFileLogger(writer, model);
        }
    }

    private void writeDiscreteTraitFileLogger(XMLWriter writer,
                                              PartitionSubstitutionModel model) {

        String prefix = model.getName() + ".";

        String fileName = options.logFileName.substring(0, options.logFileName.indexOf(".log")) + model.getName();
        fileName = (fileName.endsWith(".") ? "" : ".") + "rates.log";

        writer.writeOpenTag(LoggerParser.LOG, new Attribute[]{
                new Attribute.Default<String>(XMLParser.ID, prefix + "rateMatrixLog"),
                new Attribute.Default<String>(LoggerParser.LOG_EVERY, options.logEvery + ""),
                new Attribute.Default<String>(LoggerParser.FILE_NAME, fileName)});

        writeLogEntries(model, writer);

        writer.writeCloseTag(LoggerParser.LOG);
    }

    private void writeLogEntries(PartitionSubstitutionModel model, XMLWriter writer) {
        String prefix = model.getName() + ".";

           writer.writeIDref(ParameterParser.PARAMETER, prefix + "rates");

        if (model.isActivateBSSVS()) { //If "BSSVS" is not activated, rateIndicator should not be there.
            writer.writeIDref(ParameterParser.PARAMETER, prefix + "indicators");
            writer.writeIDref(SumStatisticParser.SUM_STATISTIC, prefix + "nonZeroRates");
        }
    }

    private void writeTreeLogEntries(PartitionTreeModel treeModel, XMLWriter writer) {
        for (AbstractPartitionData partitionData : options.getAllPartitionData(GeneralDataType.INSTANCE)) {
            if (partitionData.getPartitionTreeModel() == treeModel) {
                writer.writeIDref(AncestralStateTreeLikelihoodParser.RECONSTRUCTING_TREE_LIKELIHOOD,
                        partitionData.getPrefix() + TreeLikelihoodParser.TREE_LIKELIHOOD);
            }
        }
    }


}