/*
 * DiscreteTraitsComponentGenerator.java
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

package dr.app.beauti.components.discrete;

import dr.evomodelxml.treelikelihood.MarkovJumpsTreeLikelihoodParser;
import dr.app.beauti.components.ancestralstates.AncestralStatesComponentOptions;
import dr.app.beauti.generator.BaseComponentGenerator;
import dr.app.beauti.generator.ClockModelGenerator;
import dr.app.beauti.generator.ComponentGenerator;
import dr.app.beauti.options.*;
import dr.app.beauti.util.XMLWriter;
import dr.evolution.datatype.GeneralDataType;
import dr.inference.distribution.BinomialLikelihood;
import dr.inference.model.DesignMatrix;
import dr.inference.operators.MultivariateNormalOperator;
import dr.inferencexml.distribution.BinomialLikelihoodParser;
import dr.inferencexml.distribution.GeneralizedLinearModelParser;
import dr.oldevomodel.sitemodel.SiteModel;
import dr.oldevomodel.substmodel.AbstractSubstitutionModel;
import dr.evomodel.tree.TreeModel;
import dr.oldevomodelxml.sitemodel.GammaSiteModelParser;
import dr.oldevomodelxml.substmodel.ComplexSubstitutionModelParser;
import dr.oldevomodelxml.substmodel.FrequencyModelParser;
import dr.oldevomodelxml.substmodel.GLMSubstitutionModelParser;
import dr.oldevomodelxml.substmodel.GeneralSubstitutionModelParser;
import dr.oldevomodelxml.treelikelihood.AncestralStateTreeLikelihoodParser;
import dr.oldevomodelxml.treelikelihood.TreeLikelihoodParser;
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

import static dr.evomodelxml.substmodel.ComplexSubstitutionModelParser.ROOT_FREQUENCIES;

/**
 * @author Andrew Rambaut
 * @version $Id$
 */
public class DiscreteTraitsComponentGenerator extends BaseComponentGenerator {

    public DiscreteTraitsComponentGenerator(final BeautiOptions options) {
        super(options);
    }

    @Override
    public void checkOptions() throws GeneratorException {
        super.checkOptions();

        for (PartitionSubstitutionModel model : options.getPartitionSubstitutionModels(GeneralDataType.INSTANCE)) {
            if (model.getDiscreteSubstType() == DiscreteSubstModelType.GLM_SUBST) {
                if (model.getTraitData().getIncludedPredictors().size() < 1) {
                    throw new GeneratorException("The GLM model for trait, " + model.getTraitData().getName() + ", has no predictors included.");
                }
            }
        }
    }

    public boolean usesInsertionPoint(final InsertionPoint point) {
        if (options.getDataPartitions(GeneralDataType.INSTANCE).size() == 0) {
            // Empty, so do nothing
            return false;
        }

        boolean hasGLM = false;
        for (PartitionSubstitutionModel model : options.getPartitionSubstitutionModels(GeneralDataType.INSTANCE)) {
            if (model.getDiscreteSubstType() == DiscreteSubstModelType.GLM_SUBST) {
                hasGLM = true;
            }
        }

        switch (point) {
            case AFTER_PATTERNS:
            case AFTER_TREE_LIKELIHOOD:
            case IN_MCMC_LIKELIHOOD:
            case IN_SCREEN_LOG:
            case IN_FILE_LOG_PARAMETERS:
            case IN_FILE_LOG_LIKELIHOODS:
            case AFTER_FILE_LOG:
                return true;
            case IN_OPERATORS:
                return hasGLM;
            case IN_MCMC_PRIOR:
                return hasGLM || hasBSSVS();
            default:
                return false;
        }
    }

    protected void generate(final InsertionPoint point, final Object item, String prefix, final XMLWriter writer) {
        DiscreteTraitsComponentOptions comp = (DiscreteTraitsComponentOptions)options.getComponentOptions(DiscreteTraitsComponentOptions.class);

        switch (point) {
            case AFTER_PATTERNS:
                writeDiscreteTraitPatterns(writer, comp);
                break;

            case AFTER_TREE_LIKELIHOOD:
                writeDiscreteTraitsModels(writer);
                writeTreeLikelihoods(writer, comp);
                break;

            case IN_OPERATORS:
                for (PartitionSubstitutionModel model : options.getPartitionSubstitutionModels(GeneralDataType.INSTANCE)) {
                    if (model.getDiscreteSubstType() == DiscreteSubstModelType.GLM_SUBST) {
                        writeGLMCoefficientOperator(model, writer);
                    }
                }
                break;

            case IN_MCMC_PRIOR:
                for (PartitionSubstitutionModel model : options.getPartitionSubstitutionModels(GeneralDataType.INSTANCE)) {
                    if (model.getDiscreteSubstType() == DiscreteSubstModelType.GLM_SUBST) {
                        writeGLMBinomialLikelihood(model, writer);
                    }
                }
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
        String traitName = partition.getTraits().get(0).getName();

        writer.writeComment("Data pattern for discrete trait, '" + traitName + "'");

        String prefix = partition.getName() + ".";
        // <attributePatterns>
        writer.writeOpenTag(AttributePatternsParser.ATTRIBUTE_PATTERNS, new Attribute[]{
                new Attribute.Default<String>(XMLParser.ID, prefix + "pattern"),
                new Attribute.Default<String>(AttributePatternsParser.ATTRIBUTE, traitName)});
        String prefix2 = "";
        if (!options.hasIdenticalTaxa()) {
            prefix2 = partition.getPartitionTreeModel().getPrefix();
        }
        writer.writeIDref(TaxaParser.TAXA, prefix2 + TaxaParser.TAXA);
        writer.writeIDref(GeneralDataTypeParser.GENERAL_DATA_TYPE,
                partition.getPartitionSubstitutionModel().getName() + ".dataType");
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

            writeDiscreteFrequencyModel(prefix, prefix, stateCount, true, writer);

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

            writeDiscreteFrequencyModel(prefix, prefix, stateCount, true, writer);

            writer.writeCloseTag(GeneralSubstitutionModelParser.FREQUENCIES);

            //---------------- rates and indicators -----------------
            writeRatesAndIndicators(model, stateCount * (stateCount - 1), null, writer);

            writer.writeCloseTag(GeneralSubstitutionModelParser.GENERAL_SUBSTITUTION_MODEL);
        } else if (model.getDiscreteSubstType() == DiscreteSubstModelType.GLM_SUBST) {
            writer.writeComment("GLM substitution model");

            writer.writeOpenTag(GLMSubstitutionModelParser.GLM_SUBSTITUTION_MODEL, new Attribute[] {
                    new Attribute.Default<String>(XMLParser.ID, prefix + AbstractSubstitutionModel.MODEL)
            });

            writer.writeIDref(GeneralDataTypeParser.GENERAL_DATA_TYPE, prefix + "dataType");

            writer.writeOpenTag(ROOT_FREQUENCIES);

            writeDiscreteFrequencyModel(prefix, prefix, stateCount, true, writer);

            writer.writeCloseTag(ROOT_FREQUENCIES);

            writer.writeOpenTag(GeneralizedLinearModelParser.GLM_LIKELIHOOD, new Attribute[] {
                    new Attribute.Default<String>("family", "logLinear"),
                    new Attribute.Default<String>("checkIdentifiability", "true")
            });

            writer.writeOpenTag(GeneralizedLinearModelParser.INDEPENDENT_VARIABLES);

            writeParameter(options.getParameter(prefix + "coefficients"), 1, writer);

            writeParameter(GeneralizedLinearModelParser.INDICATOR, prefix + "coefIndicators", 1, writer);

            writer.writeOpenTag(DesignMatrix.DESIGN_MATRIX, new Attribute[] {
                    new Attribute.Default<String>("id", prefix + DesignMatrix.DESIGN_MATRIX)
            } );

            for (Predictor predictor : model.getTraitData().getIncludedPredictors()) {
                if (predictor.getType() == Predictor.Type.ORIGIN_VECTOR || predictor.getType() == Predictor.Type.BOTH_VECTOR) {
                    writer.writeTag(ParameterParser.PARAMETER, new Attribute[]{
                            new Attribute.Default<String>(XMLParser.ID, prefix + predictor + "_origin"),
                            new Attribute.Default<String>(ParameterParser.VALUE, predictor.getValueString(Predictor.Type.ORIGIN_VECTOR))
                    }, true);
                }
                if (predictor.getType() == Predictor.Type.DESTINATION_VECTOR || predictor.getType() == Predictor.Type.BOTH_VECTOR) {
                    writer.writeTag(ParameterParser.PARAMETER, new Attribute[]{
                            new Attribute.Default<String>(XMLParser.ID, prefix + predictor + "_destination"),
                            new Attribute.Default<String>(ParameterParser.VALUE, predictor.getValueString(Predictor.Type.DESTINATION_VECTOR))
                    }, true);
                }
                if (predictor.getType() == Predictor.Type.MATRIX) {
                    writer.writeTag(ParameterParser.PARAMETER, new Attribute[]{
                            new Attribute.Default<String>(XMLParser.ID, prefix + predictor),
                            new Attribute.Default<String>(ParameterParser.VALUE, predictor.getValueString(Predictor.Type.MATRIX))
                    }, true);
                }
            }
            writer.writeCloseTag(DesignMatrix.DESIGN_MATRIX);

            writer.writeCloseTag(GeneralizedLinearModelParser.INDEPENDENT_VARIABLES);

            writer.writeCloseTag(GeneralizedLinearModelParser.GLM_LIKELIHOOD);

            writer.writeCloseTag(GLMSubstitutionModelParser.GLM_SUBSTITUTION_MODEL);

            writer.writeComment("GLM: statistic that returns the product of the coefficients and the respective indicators for the predictors");

            writer.writeOpenTag(ProductStatisticParser.PRODUCT_STATISTIC, new Attribute[] {
                    new Attribute.Default<String>(XMLParser.ID, prefix + "coefficientsTimesIndicators")
            });
            writer.writeTag(ParameterParser.PARAMETER, new Attribute[]{
                    new Attribute.Default<String>(XMLParser.IDREF, prefix + "coefficients")
            }, true);
            writer.writeTag(ParameterParser.PARAMETER, new Attribute[]{
                    new Attribute.Default<String>(XMLParser.IDREF, prefix + "coefIndicators")
            }, true);
            writer.writeCloseTag(ProductStatisticParser.PRODUCT_STATISTIC);

            writer.writeComment("GLM: the number of predictors included in the model (non-zero indicators)");
            writer.writeOpenTag(SumStatisticParser.SUM_STATISTIC, new Attribute[]{
                    new Attribute.Default<String>(XMLParser.ID, prefix + "includedPredictors")});
            writer.writeIDref(ParameterParser.PARAMETER, prefix + "coefIndicators");
            writer.writeCloseTag(SumStatisticParser.SUM_STATISTIC);

        } else {
            throw new IllegalArgumentException("Unknown discreteSubstType");
        }

        if (model.isActivateBSSVS()) {
            // If "BSSVS" is not activated, rateIndicator should not be there.
            writeStatisticModel(model, writer);
        }
    }

    private void writeDiscreteTraitsSubstitutionModelReferences(XMLWriter writer) {
        for (PartitionSubstitutionModel model : options.getPartitionSubstitutionModels(GeneralDataType.INSTANCE)) {
            if (model.getDiscreteSubstType() == DiscreteSubstModelType.GLM_SUBST) {
                // not strictly necessary but makes the XML consistent
                writer.writeIDref(GLMSubstitutionModelParser.GLM_SUBSTITUTION_MODEL, model.getName() + "." + AbstractSubstitutionModel.MODEL);
            } else {
                writer.writeIDref(GeneralSubstitutionModelParser.GENERAL_SUBSTITUTION_MODEL, model.getName() + "." + AbstractSubstitutionModel.MODEL);
            }
        }
    }

    private boolean hasBSSVS() {
        for (PartitionSubstitutionModel model : options.getPartitionSubstitutionModels(GeneralDataType.INSTANCE)) {
            if (model.isActivateBSSVS())
                return true;
        }
        return false;
    }

    private void writeDiscreteFrequencyModel(String prefix, String dataTypePrefix, int stateCount, Boolean normalize, XMLWriter writer) {
        if (normalize == null) {
            writer.writeOpenTag(FrequencyModelParser.FREQUENCY_MODEL, new Attribute[]{
                    new Attribute.Default<String>(XMLParser.ID, prefix + FrequencyModelParser.FREQUENCY_MODEL)});
        } else {
            writer.writeOpenTag(FrequencyModelParser.FREQUENCY_MODEL, new Attribute[]{
                    new Attribute.Default<String>(XMLParser.ID, prefix + FrequencyModelParser.FREQUENCY_MODEL),
                    new Attribute.Default<Boolean>(FrequencyModelParser.NORMALIZE, normalize)});
        }

        writer.writeIDref(GeneralDataTypeParser.GENERAL_DATA_TYPE, dataTypePrefix + "dataType");

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
        if (model.getDiscreteSubstType() == DiscreteSubstModelType.GLM_SUBST) {
            // not strictly necessary but makes the XML consistent
            writer.writeIDref(GLMSubstitutionModelParser.GLM_SUBSTITUTION_MODEL, prefix + AbstractSubstitutionModel.MODEL);
        } else {
            writer.writeIDref(GeneralSubstitutionModelParser.GENERAL_SUBSTITUTION_MODEL, prefix + AbstractSubstitutionModel.MODEL);
        }
        writer.writeCloseTag(GammaSiteModelParser.SUBSTITUTION_MODEL);

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
        writeParameter(options.getParameter(prefix + "rates"), dimension, writer);

        writer.writeCloseTag(GeneralSubstitutionModelParser.RATES);

        if (model.isActivateBSSVS()) { //If "BSSVS" is not activated, rateIndicator should not be there.
            writer.writeOpenTag(GeneralSubstitutionModelParser.INDICATOR);
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

        if (partition.getPartitionSubstitutionModel().getDiscreteSubstType() == DiscreteSubstModelType.GLM_SUBST) {
            // not strictly necessary but makes the XML consistent
            writer.writeIDref(GLMSubstitutionModelParser.GLM_SUBSTITUTION_MODEL, prefix + AbstractSubstitutionModel.MODEL);
        } else {
            writer.writeIDref(GeneralSubstitutionModelParser.GENERAL_SUBSTITUTION_MODEL, prefix + AbstractSubstitutionModel.MODEL);
        }

        writer.writeIDref(GeneralSubstitutionModelParser.GENERAL_SUBSTITUTION_MODEL, substModel.getName() + "." + AbstractSubstitutionModel.MODEL);

        ClockModelGenerator.writeBranchRatesModelRef(clockModel, writer);

        if (substModel.getDiscreteSubstType() == DiscreteSubstModelType.ASYM_SUBST) {
            int stateCount = options.getStatesForDiscreteModel(substModel).size();
            writer.writeComment("The root state frequencies");
            writeDiscreteFrequencyModel(partition.getPrefix() + "root.", substModel.getName() + ".", stateCount, true, writer);
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

    private void writeGLMCoefficientOperator(PartitionSubstitutionModel model, XMLWriter writer) {
        writer.writeOpenTag(MultivariateNormalOperator.MVN_OPERATOR,
                new Attribute[]{
                        new Attribute.Default<Double>(MultivariateNormalOperator.SCALE_FACTOR, 1.0),
                        new Attribute.Default<Double>(MultivariateNormalOperator.WEIGHT,5.0),
                        new Attribute.Default<String>(MultivariateNormalOperator.FORM_XTX, "true")
                }
        );
        String prefix = model.getName() + ".";

        writeParameterRef(prefix + "coefficients", writer);

        writer.writeOpenTag(MultivariateNormalOperator.VARIANCE_MATRIX);
        writeParameterRef(prefix + DesignMatrix.DESIGN_MATRIX, writer);
        writer.writeCloseTag(MultivariateNormalOperator.VARIANCE_MATRIX);

        writer.writeCloseTag(MultivariateNormalOperator.MVN_OPERATOR);
    }

    private void writeGLMBinomialLikelihood(PartitionSubstitutionModel model, XMLWriter writer) {
        double proportion = 1.0 - Math.exp(Math.log(0.5) / model.getTraitData().getIncludedPredictorCount());

        String prefix = model.getName() + ".";

        writer.writeComment("Using the binomialLikelihood we specify a 50% prior mass on no predictors being included.");
        writer.writeOpenTag(BinomialLikelihood.BINOMIAL_LIKELIHOOD);
        writer.writeOpenTag(BinomialLikelihoodParser.PROPORTION);
        writer.writeTag("parameter", new Attribute.Default<Double>("value", proportion), true);
        writer.writeCloseTag(BinomialLikelihoodParser.PROPORTION);
        // the dimension of this parameter will be set automatically to be the same as the counts.
        writer.writeOpenTag(BinomialLikelihoodParser.TRIALS);
        writer.writeTag("parameter", new Attribute[]{
                new Attribute.Default<Double>("value", 1.0)
        }, true);
        writer.writeCloseTag(BinomialLikelihoodParser.TRIALS);
        writer.writeOpenTag(BinomialLikelihoodParser.COUNTS);
        writeParameterRef(prefix + "coefIndicators", writer);
        writer.writeCloseTag(BinomialLikelihoodParser.COUNTS);

        writer.writeCloseTag(BinomialLikelihood.BINOMIAL_LIKELIHOOD);
    }

    public void writeScreenLogEntries(XMLWriter writer) {
        for (PartitionSubstitutionModel model : options.getPartitionSubstitutionModels(GeneralDataType.INSTANCE)) {
            String prefix = model.getName() + ".";

            if (model.isActivateBSSVS()) { //If "BSSVS" is not activated, rateIndicator should not be there.
                writer.writeOpenTag(ColumnsParser.COLUMN,
                        new Attribute[]{
                                new Attribute.Default<String>(ColumnsParser.LABEL, prefix + "nonZeroRates"),
                                new Attribute.Default<String>(ColumnsParser.DECIMAL_PLACES, "0"),
                                new Attribute.Default<String>(ColumnsParser.WIDTH, "6")
                        }
                );

                writer.writeIDref(SumStatisticParser.SUM_STATISTIC, prefix + "nonZeroRates");

                writer.writeCloseTag(ColumnsParser.COLUMN);
            }
            if (model.getDiscreteSubstType() == DiscreteSubstModelType.GLM_SUBST) {
                writer.writeOpenTag(ColumnsParser.COLUMN,
                        new Attribute[]{
                                new Attribute.Default<String>(ColumnsParser.LABEL, prefix + "incPredictors"),
                                new Attribute.Default<String>(ColumnsParser.DECIMAL_PLACES, "0"),
                                new Attribute.Default<String>(ColumnsParser.WIDTH, "6")
                        }
                );

                writer.writeIDref(SumStatisticParser.SUM_STATISTIC, prefix + "includedPredictors");

                writer.writeCloseTag(ColumnsParser.COLUMN);
            }
        }
    }

    private void writeFileLogEntries(XMLWriter writer) {
        for (PartitionSubstitutionModel model : options.getPartitionSubstitutionModels(GeneralDataType.INSTANCE)) {
            String prefix = model.getName() + ".";

            if (model.getDiscreteSubstType() == DiscreteSubstModelType.GLM_SUBST) {
                writer.writeIDref(SumStatisticParser.SUM_STATISTIC, prefix + "includedPredictors");
                writer.writeIDref(ProductStatisticParser.PRODUCT_STATISTIC, prefix + "coefficientsTimesIndicators");
            } else {
                writer.writeIDref(ParameterParser.PARAMETER, prefix + "rates");

                if (model.isActivateBSSVS()) { //If "BSSVS" is not activated, rateIndicator should not be there.
                    writer.writeIDref(ParameterParser.PARAMETER, prefix + "indicators");
                    writer.writeIDref(SumStatisticParser.SUM_STATISTIC, prefix + "nonZeroRates");
                }
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

        if (model.getDiscreteSubstType() == DiscreteSubstModelType.GLM_SUBST) {
            writer.writeOpenTag(LoggerParser.LOG, new Attribute[]{
                    new Attribute.Default<String>(XMLParser.ID, prefix + "glmLog"),
                    new Attribute.Default<String>(LoggerParser.LOG_EVERY, options.logEvery + ""),
                    new Attribute.Default<String>(LoggerParser.FILE_NAME, prefix + ".glm.log")});
        } else {
            writer.writeOpenTag(LoggerParser.LOG, new Attribute[]{
                    new Attribute.Default<String>(XMLParser.ID, prefix + "rateMatrixLog"),
                    new Attribute.Default<String>(LoggerParser.LOG_EVERY, options.logEvery + ""),
                    new Attribute.Default<String>(LoggerParser.FILE_NAME, prefix + ".rates.log")});
        }

        writeLogEntries(model, writer);

        writer.writeCloseTag(LoggerParser.LOG);
    }

    private void writeLogEntries(PartitionSubstitutionModel model, XMLWriter writer) {
        String prefix = model.getName() + ".";

        if (model.getDiscreteSubstType() == DiscreteSubstModelType.GLM_SUBST) {
            writer.writeIDref(ParameterParser.PARAMETER, prefix + "coefficients");
            writer.writeIDref(ParameterParser.PARAMETER, prefix + "coefIndicators");
            writer.writeIDref(ProductStatisticParser.PRODUCT_STATISTIC, prefix + "coefficientsTimesIndicators");
        } else {
            writer.writeIDref(ParameterParser.PARAMETER, prefix + "rates");

            if (model.isActivateBSSVS()) { //If "BSSVS" is not activated, rateIndicator should not be there.
                writer.writeIDref(ParameterParser.PARAMETER, prefix + "indicators");
                writer.writeIDref(SumStatisticParser.SUM_STATISTIC, prefix + "nonZeroRates");
            }
        }
    }

}