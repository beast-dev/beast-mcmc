/*
 * DiscreteTraitsComponentGenerator.java
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

package dr.app.beauti.components.discrete;

import dr.app.beagle.evomodel.parsers.MarkovJumpsTreeLikelihoodParser;
import dr.app.beauti.components.ancestralstates.AncestralStatesComponentOptions;
import dr.app.beauti.generator.BaseComponentGenerator;
import dr.app.beauti.generator.ComponentGenerator;
import dr.app.beauti.options.*;
import dr.app.beauti.util.XMLWriter;
import dr.evolution.datatype.GeneralDataType;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.sitemodel.SiteModel;
import dr.evomodel.substmodel.AbstractSubstitutionModel;
import dr.evomodel.tree.TreeModel;
import dr.evomodelxml.branchratemodel.DiscretizedBranchRatesParser;
import dr.evomodelxml.branchratemodel.RandomLocalClockModelParser;
import dr.evomodelxml.branchratemodel.StrictClockBranchRatesParser;
import dr.evomodelxml.clock.ACLikelihoodParser;
import dr.evomodelxml.sitemodel.GammaSiteModelParser;
import dr.evomodelxml.substmodel.ComplexSubstitutionModelParser;
import dr.evomodelxml.substmodel.FrequencyModelParser;
import dr.evomodelxml.substmodel.GeneralSubstitutionModelParser;
import dr.evomodelxml.treelikelihood.AncestralStateTreeLikelihoodParser;
import dr.evomodelxml.treelikelihood.TreeLikelihoodParser;
import dr.evoxml.AttributePatternsParser;
import dr.evoxml.GeneralDataTypeParser;
import dr.evoxml.TaxaParser;
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
        if (options.getDataPartitions(GeneralDataType.INSTANCE).size() == 0) {
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
//            case IN_TREES_LOG:
                return true;
            case IN_MCMC_PRIOR:
                return hasBSSVS();
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
                writeDiscreteTraitsModels(writer);
                writeTreeLikelihoods(writer, comp);
                break;

            case IN_MCMC_PRIOR:
                writeDiscreteTraitsSubstitutionModelReferences(writer);
                break;

            case IN_MCMC_LIKELIHOOD:
                writeTreeLikelihoodReferences(writer);
                break;

            case IN_SCREEN_LOG:
                writeScreenLogEntries(writer);
                break;

            case IN_FILE_LOG_PARAMETERS:
                writeFileLogEntries(writer);
                break;

            case IN_FILE_LOG_LIKELIHOODS:
                writeTreeLikelihoodReferences(writer);
                break;

            case AFTER_FILE_LOG:
                writeDiscreteTraitFileLoggers(writer);

                // This was for ancestral state writing which is now handled by the ancestral state component
//            case IN_TREES_LOG:
//                writeTreeLogEntries((PartitionTreeModel)item, writer);
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

        boolean first = true;
        for (PartitionSubstitutionModel model : options.getPartitionSubstitutionModels(GeneralDataType.INSTANCE)) {
            if (!first) {
            writer.writeBlankLine();
            } else {
                first = false;
        }
            writeGeneralDataType(model, writer);
        }

        // now create an attribute pattern for each trait that uses it
        for (AbstractPartitionData partition : options.getDataPartitions(GeneralDataType.INSTANCE)) {
            if (partition.getTraits() != null) {
                writer.writeBlankLine();
                writeAttributePatterns(partition, writer);
            }
        }

    }

    private void writeGeneralDataType(PartitionSubstitutionModel model, XMLWriter writer) {

        writer.writeComment("general data type for discrete trait model, '" + model.getName() + "'");

        Set<String> states = options.getStatesForDiscreteModel(model);
        String prefix = model.getName() + ".";

        // <generalDataType>
        writer.writeOpenTag(GeneralDataTypeParser.GENERAL_DATA_TYPE, new Attribute[]{
                new Attribute.Default<String>(XMLParser.ID, prefix + "dataType")});

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

        String prefix = partition.getName() + ".";
        // <attributePatterns>
        writer.writeOpenTag(AttributePatternsParser.ATTRIBUTE_PATTERNS, new Attribute[]{
                new Attribute.Default<String>(XMLParser.ID, prefix + "pattern"),
                new Attribute.Default<String>(AttributePatternsParser.ATTRIBUTE, partition.getTraits().get(0).getName())});
        writer.writeIDref(TaxaParser.TAXA, TaxaParser.TAXA);
        writer.writeIDref(GeneralDataTypeParser.GENERAL_DATA_TYPE, prefix + "dataType");
        writer.writeCloseTag(AttributePatternsParser.ATTRIBUTE_PATTERNS);
    }

    private void writeDiscreteTraitsModels(XMLWriter writer) {
        for (PartitionSubstitutionModel model : options.getPartitionSubstitutionModels(GeneralDataType.INSTANCE)) {
            writeDiscreteTraitsSubstitutionModel(model, writer);
            writeDiscreteTraitsSiteModel(model, writer);
        }
    }

    private void writeDiscreteTraitsSubstitutionModel(PartitionSubstitutionModel model, XMLWriter writer) {

        int stateCount = options.getStatesForDiscreteModel(model).size();
        String prefix = model.getName() + ".";

        if (model.getDiscreteSubstType() == DiscreteSubstModelType.SYM_SUBST) {
            writer.writeComment("symmetric CTMC model for discrete state reconstructions");

            writer.writeOpenTag(GeneralSubstitutionModelParser.GENERAL_SUBSTITUTION_MODEL, new Attribute[]{
                    new Attribute.Default<String>(XMLParser.ID, prefix + AbstractSubstitutionModel.MODEL)});

            writer.writeIDref(GeneralDataTypeParser.GENERAL_DATA_TYPE, prefix +  "dataType");

            writer.writeOpenTag(GeneralSubstitutionModelParser.FREQUENCIES);

            writeDiscreteFrequencyModel(model, stateCount, true, writer);

            writer.writeCloseTag(GeneralSubstitutionModelParser.FREQUENCIES);

            //---------------- rates and indicators -----------------

            writeRatesAndIndicators(model, stateCount * (stateCount - 1) / 2, null, writer);
            writer.writeCloseTag(GeneralSubstitutionModelParser.GENERAL_SUBSTITUTION_MODEL);
        } else if (model.getDiscreteSubstType() == DiscreteSubstModelType.ASYM_SUBST) {
            writer.writeComment("asymmetric CTMC model for discrete state reconstructions");

            writer.writeOpenTag(GeneralSubstitutionModelParser.GENERAL_SUBSTITUTION_MODEL, new Attribute[]{
                    new Attribute.Default<String>(XMLParser.ID, prefix + AbstractSubstitutionModel.MODEL),
                    new Attribute.Default<Boolean>(ComplexSubstitutionModelParser.RANDOMIZE, false)});

            writer.writeIDref(GeneralDataTypeParser.GENERAL_DATA_TYPE, prefix + "dataType");

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

    private void writeDiscreteTraitsSubstitutionModelReferences(XMLWriter writer) {
        for (PartitionSubstitutionModel model : options.getPartitionSubstitutionModels(GeneralDataType.INSTANCE)) {
            writer.writeIDref(GeneralSubstitutionModelParser.GENERAL_SUBSTITUTION_MODEL, model.getName() + "." + AbstractSubstitutionModel.MODEL);
        }
    }

    private boolean hasBSSVS() {
        for (PartitionSubstitutionModel model : options.getPartitionSubstitutionModels(GeneralDataType.INSTANCE)) {
            if (model.isActivateBSSVS())
                return true;
        }
        return false;
    }

    private void writeDiscreteFrequencyModel(PartitionSubstitutionModel model, int stateCount, Boolean normalize, XMLWriter writer) {
        String prefix = model.getName() + ".";
        if (normalize == null) {
            writer.writeOpenTag(FrequencyModelParser.FREQUENCY_MODEL, new Attribute[]{
                    new Attribute.Default<String>(XMLParser.ID, prefix + FrequencyModelParser.FREQUENCY_MODEL)});
        } else {
            writer.writeOpenTag(FrequencyModelParser.FREQUENCY_MODEL, new Attribute[]{
                    new Attribute.Default<String>(XMLParser.ID, prefix + FrequencyModelParser.FREQUENCY_MODEL),
                    new Attribute.Default<Boolean>(FrequencyModelParser.NORMALIZE, normalize)});
        }

        writer.writeIDref(GeneralDataTypeParser.GENERAL_DATA_TYPE, prefix + "dataType");

        writer.writeOpenTag(FrequencyModelParser.FREQUENCIES);
        writeParameter(prefix + "frequencies", stateCount, Double.NaN, Double.NaN, Double.NaN, writer);
        writer.writeCloseTag(FrequencyModelParser.FREQUENCIES);

        writer.writeCloseTag(FrequencyModelParser.FREQUENCY_MODEL);
    }

    private void writeDiscreteTraitsSiteModel(PartitionSubstitutionModel model, XMLWriter writer) {
        String prefix = model.getName() + ".";
        writer.writeOpenTag(SiteModel.SITE_MODEL, new Attribute[]{
                new Attribute.Default<String>(XMLParser.ID, prefix + SiteModel.SITE_MODEL)});

        writer.writeOpenTag(GammaSiteModelParser.SUBSTITUTION_MODEL);
        writer.writeIDref(GeneralSubstitutionModelParser.GENERAL_SUBSTITUTION_MODEL, prefix + AbstractSubstitutionModel.MODEL);
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
        options.getParameter(prefix + "rates").isFixed = true;
        writeParameter(options.getParameter(prefix + "rates"), dimension, writer);

        writer.writeCloseTag(GeneralSubstitutionModelParser.RATES);

        if (model.isActivateBSSVS()) { //If "BSSVS" is not activated, rateIndicator should not be there.
            writer.writeOpenTag(GeneralSubstitutionModelParser.INDICATOR);
            options.getParameter(prefix + "indicators").isFixed = true;
            writeParameter(options.getParameter(prefix + "indicators"), dimension, writer);
            writer.writeCloseTag(GeneralSubstitutionModelParser.INDICATOR);
        }
    }

    private void writeStatisticModel(PartitionSubstitutionModel model, XMLWriter writer) {
        String prefix = model.getName() + ".";

        writer.writeOpenTag(SumStatisticParser.SUM_STATISTIC, new Attribute[]{
                new Attribute.Default<String>(XMLParser.ID, prefix + "nonZeroRates"),
                new Attribute.Default<Boolean>(SumStatisticParser.ELEMENTWISE, true)});
        writer.writeIDref(ParameterParser.PARAMETER, prefix + "indicators");
        writer.writeCloseTag(SumStatisticParser.SUM_STATISTIC);

        writer.writeOpenTag(ProductStatisticParser.PRODUCT_STATISTIC, new Attribute[]{
                new Attribute.Default<String>(XMLParser.ID, prefix + "actualRates"),
                new Attribute.Default<Boolean>(SumStatisticParser.ELEMENTWISE, false)});
        writer.writeIDref(ParameterParser.PARAMETER, prefix + "indicators");
        writer.writeIDref(ParameterParser.PARAMETER, prefix + "rates");
        writer.writeCloseTag(ProductStatisticParser.PRODUCT_STATISTIC);
    }

    private void writeTreeLikelihoods(XMLWriter writer,
                                               DiscreteTraitsComponentOptions component) {
        // generate tree likelihoods for discrete trait partitions
        if (options.hasDiscreteTraitPartition()) {
            writer.writeComment("Likelihood for tree given discrete trait data");
        }
        for (AbstractPartitionData partition : options.dataPartitions) {
            if (partition.getTraits() != null) {
                TraitData trait = partition.getTraits().get(0);
                if (trait.getTraitType() == TraitData.TraitType.DISCRETE) {
                    writeTreeLikelihood(partition, writer);
                }
            }
        }
    }

    /**
     * Write Tree Likelihood
     *
     * @param partition PartitionData
     * @param writer    XMLWriter
     */
    public void writeTreeLikelihood(AbstractPartitionData partition, XMLWriter writer) {
        String prefix = partition.getName() + ".";

        PartitionSubstitutionModel substModel = partition.getPartitionSubstitutionModel();
        PartitionTreeModel treeModel = partition.getPartitionTreeModel();
        PartitionClockModel clockModel = partition.getPartitionClockModel();

        AncestralStatesComponentOptions ancestralStatesOptions = (AncestralStatesComponentOptions)options.getComponentOptions(AncestralStatesComponentOptions.class);
        String treeLikelihoodTag = TreeLikelihoodParser.ANCESTRAL_TREE_LIKELIHOOD;
        if (ancestralStatesOptions.isCountingStates(partition)) {
            treeLikelihoodTag = MarkovJumpsTreeLikelihoodParser.MARKOV_JUMP_TREE_LIKELIHOOD;
        }

        writer.writeOpenTag(treeLikelihoodTag, new Attribute[]{
                new Attribute.Default<String>(XMLParser.ID, prefix + TreeLikelihoodParser.TREE_LIKELIHOOD),
                new Attribute.Default<String>(AncestralStateTreeLikelihoodParser.RECONSTRUCTION_TAG_NAME, prefix + AncestralStateTreeLikelihoodParser.RECONSTRUCTION_TAG),
        });

        writer.writeIDref(AttributePatternsParser.ATTRIBUTE_PATTERNS, prefix + "pattern");
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

        getCallingGenerator().generateInsertionPoint(ComponentGenerator.InsertionPoint.IN_TREE_LIKELIHOOD, partition, writer);

        writer.writeCloseTag(treeLikelihoodTag);
    }

    private void writeTreeLikelihoodReferences(XMLWriter writer) {
        for (AbstractPartitionData partition : options.dataPartitions) {
            if (partition.getTraits() != null) {
                AncestralStatesComponentOptions ancestralStatesOptions = (AncestralStatesComponentOptions)options.getComponentOptions(AncestralStatesComponentOptions.class);
                String treeLikelihoodTag = TreeLikelihoodParser.ANCESTRAL_TREE_LIKELIHOOD;
                if (ancestralStatesOptions.isCountingStates(partition)) {
                    treeLikelihoodTag = MarkovJumpsTreeLikelihoodParser.MARKOV_JUMP_TREE_LIKELIHOOD;
                }

                TraitData trait = partition.getTraits().get(0);
                String prefix = partition.getName() + ".";
                if (trait.getTraitType() == TraitData.TraitType.DISCRETE) {
                    writer.writeIDref(treeLikelihoodTag,
                            prefix + TreeLikelihoodParser.TREE_LIKELIHOOD);
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
                                new Attribute.Default<String>(ColumnsParser.LABEL, prefix + "nonZeroRates"),
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

        String prefix = options.fileNameStem + "." + model.getName();

        writer.writeOpenTag(LoggerParser.LOG, new Attribute[]{
                new Attribute.Default<String>(XMLParser.ID, prefix + "rateMatrixLog"),
                new Attribute.Default<String>(LoggerParser.LOG_EVERY, options.logEvery + ""),
                new Attribute.Default<String>(LoggerParser.FILE_NAME, prefix + ".rates.log")});

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

    // This was for ancestral state writing which is now handled by the ancestral state component
//    private void writeTreeLogEntries(PartitionTreeModel treeModel, XMLWriter writer) {
//        for (AbstractPartitionData partitionData : options.getDataPartitions(GeneralDataType.INSTANCE)) {
//            if (partitionData.getPartitionTreeModel() == treeModel) {
//                String prefix = partitionData.getName() + ".";
//                writer.writeIDref(AncestralStateTreeLikelihoodParser.RECONSTRUCTING_TREE_LIKELIHOOD,
//                        prefix + OldTreeLikelihoodParser.TREE_LIKELIHOOD);
//            }
//        }
//    }


}