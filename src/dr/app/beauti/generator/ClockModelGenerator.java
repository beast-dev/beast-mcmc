/*
 * BranchRatesModelGenerator.java
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

import dr.app.beauti.components.ComponentFactory;
import dr.app.beauti.options.*;
import dr.app.beauti.types.ClockDistributionType;
import dr.app.beauti.types.ClockType;
import dr.app.beauti.types.OperatorSetType;
import dr.app.beauti.util.XMLWriter;
import dr.evolution.util.Taxa;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.inference.model.Statistic;
import dr.inference.model.StatisticParser;
import dr.oldevomodel.clock.RateEvolutionLikelihood;
import dr.evomodel.tree.TreeModel;
import dr.evomodelxml.branchratemodel.*;
import dr.oldevomodelxml.clock.ACLikelihoodParser;
import dr.evomodelxml.tree.RateCovarianceStatisticParser;
import dr.evomodelxml.tree.RateStatisticParser;
import dr.evomodelxml.tree.TreeModelParser;
import dr.evoxml.TaxaParser;
import dr.inference.distribution.ExponentialDistributionModel;
import dr.inference.distribution.GammaDistributionModel;
import dr.inference.model.ParameterParser;
import dr.inferencexml.distribution.LogNormalDistributionModelParser;
import dr.inferencexml.model.CompoundParameterParser;
import dr.inferencexml.model.SumStatisticParser;
import dr.util.Attribute;
import dr.xml.XMLParser;

import java.util.List;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 */
public class ClockModelGenerator extends Generator {

    public ClockModelGenerator(BeautiOptions options, ComponentFactory[] components) {
        super(options, components);
    }

    /**
     * Write the relaxed clock branch rates block.
     *
     * @param clockModel  PartitionClockModel
     * @param writer the writer
     */
    public void writeBranchRatesModel(PartitionClockModel clockModel, XMLWriter writer) {

        Attribute[] attributes;

        PartitionTreeModel treeModel = clockModel.getPartitionTreeModel();

        String treePrefix = treeModel.getPrefix();
        String prefix = clockModel.getPrefix();
        String tag;
//        List<PartitionTreeModel> activeTrees = options.getPartitionTreeModels(options.getDataPartitions(clockModel));

        switch (clockModel.getClockType()) {
            case STRICT_CLOCK:
                tag = StrictClockBranchRatesParser.STRICT_CLOCK_BRANCH_RATES;
                writer.writeComment("The strict clock (Uniform rates across branches)");

                writer.writeOpenTag(
                        tag,
                        new Attribute[]{new Attribute.Default<String>(XMLParser.ID, prefix + BranchRateModel.BRANCH_RATES)}
                );
                writeParameter("rate", "clock.rate", clockModel, writer);
                writer.writeCloseTag(tag);

                writeMeanRateStatistic(writer, tag, prefix, treePrefix);

                break;

            case UNCORRELATED:

                if (clockModel.performModelAveraging()) {
                    tag = MixtureModelBranchRatesParser.MIXTURE_MODEL_BRANCH_RATES;

                    writer.writeComment("Bayesian Model Averaging (BMA) of the available relaxed clock models as described by" +
                            " Li & Drummond (2012) MBE 29:751-61.");
                    //Bayesian Model Averaging uses the continuous quantile implementation by default
                    writer.writeComment("  Continuous quantile implementation (Li & Drummond (2012) Mol Biol Evol 29:751-61)");

                    attributes = new Attribute[]{
                            new Attribute.Default<String>(XMLParser.ID,
                                    prefix  + BranchRateModel.BRANCH_RATES)};
                    writer.writeOpenTag(tag, attributes);
                    // tree
                    writer.writeIDref(TreeModel.TREE_MODEL, treePrefix + TreeModel.TREE_MODEL);

                    writer.writeOpenTag("distribution");
                    writer.writeOpenTag(LogNormalDistributionModelParser.LOGNORMAL_DISTRIBUTION_MODEL,
                            new Attribute.Default<String>(LogNormalDistributionModelParser.MEAN_IN_REAL_SPACE, "true"));
                    writeParameter("mean", ClockType.UCLD_MEAN, clockModel, writer);
                    writeParameter("stdev", ClockType.UCLD_STDEV, clockModel, writer);
                    writer.writeCloseTag(LogNormalDistributionModelParser.LOGNORMAL_DISTRIBUTION_MODEL);
                    writer.writeCloseTag("distribution");

                    writer.writeOpenTag("distribution");
                    writer.writeOpenTag(GammaDistributionModel.GAMMA_DISTRIBUTION_MODEL);
                    writeParameter("mean", ClockType.UCGD_MEAN, clockModel, writer);
                    writeParameter("shape", ClockType.UCGD_SHAPE, clockModel, writer);
                    writer.writeCloseTag(GammaDistributionModel.GAMMA_DISTRIBUTION_MODEL);
                    writer.writeCloseTag("distribution");

                    writer.writeOpenTag("distribution");
                    writer.writeOpenTag(ExponentialDistributionModel.EXPONENTIAL_DISTRIBUTION_MODEL);
                    writeParameter("mean", ClockType.UCED_MEAN, clockModel, writer);
                    writer.writeCloseTag(ExponentialDistributionModel.EXPONENTIAL_DISTRIBUTION_MODEL);
                    writer.writeCloseTag("distribution");

                    writer.writeOpenTag(MixtureModelBranchRatesParser.DISTRIBUTION_INDEX);
                    writeParameter(clockModel.getParameter("branchRates.distributionIndex"), -1, writer);
                    writer.writeCloseTag(MixtureModelBranchRatesParser.DISTRIBUTION_INDEX);

                    writer.writeOpenTag(MixtureModelBranchRatesParser.RATE_CATEGORY_QUANTILES);
                    writeParameter(clockModel.getParameter("branchRates.quantiles"), -1, writer);
                    writer.writeCloseTag(MixtureModelBranchRatesParser.RATE_CATEGORY_QUANTILES);
                    writer.writeCloseTag(tag);

                    writeMeanRateStatistic(writer, tag, prefix, treePrefix);

                    writeCoefficientOfVariationStatistic(writer, tag, prefix, treePrefix);

                    writeCovarianceStatistic(writer, tag, prefix, treePrefix);

                } else {
                    tag = DiscretizedBranchRatesParser.DISCRETIZED_BRANCH_RATES;

                    writer.writeComment("The uncorrelated relaxed clock (Drummond, Ho, Phillips & Rambaut (2006) PLoS Biology 4, e88 )");

                    if (clockModel.isContinuousQuantile()) {
                        writer.writeComment("  Continuous quantile implementation (Li & Drummond (2012) Mol Biol Evol 29:751-61)");
                        tag = ContinuousBranchRatesParser.CONTINUOUS_BRANCH_RATES;
                    }

                    attributes = new Attribute[]{
                            new Attribute.Default<String>(XMLParser.ID,
                                    prefix  + BranchRateModel.BRANCH_RATES)};
                    writer.writeOpenTag(tag, attributes);
                    // tree
                    writer.writeIDref(TreeModel.TREE_MODEL, treePrefix + TreeModel.TREE_MODEL);

                    writer.writeOpenTag("distribution");

                    switch (clockModel.getClockDistributionType()) {

                        case LOGNORMAL:
                            writer.writeOpenTag(LogNormalDistributionModelParser.LOGNORMAL_DISTRIBUTION_MODEL,
                                    new Attribute.Default<String>(LogNormalDistributionModelParser.MEAN_IN_REAL_SPACE, "true"));

                            writeParameter("mean", ClockType.UCLD_MEAN, clockModel, writer);
                            writeParameter("stdev", ClockType.UCLD_STDEV, clockModel, writer);

                            writer.writeCloseTag(LogNormalDistributionModelParser.LOGNORMAL_DISTRIBUTION_MODEL);
                            break;
                        case GAMMA:
                            writer.writeOpenTag(GammaDistributionModel.GAMMA_DISTRIBUTION_MODEL);

                            writeParameter("mean", ClockType.UCGD_MEAN, clockModel, writer);
                            writeParameter("shape", ClockType.UCGD_SHAPE, clockModel, writer);

                            writer.writeCloseTag(GammaDistributionModel.GAMMA_DISTRIBUTION_MODEL);
                            break;
                        case CAUCHY:
                            throw new UnsupportedOperationException("Uncorrelated Cauchy model not implemented yet");
//                            break;
                        case EXPONENTIAL:
                            writer.writeOpenTag(ExponentialDistributionModel.EXPONENTIAL_DISTRIBUTION_MODEL);

                            writeParameter("mean", ClockType.UCED_MEAN, clockModel, writer);

                            writer.writeCloseTag(ExponentialDistributionModel.EXPONENTIAL_DISTRIBUTION_MODEL);
                            break;
                    }

                    writer.writeCloseTag("distribution");

                    if (clockModel.isContinuousQuantile()) {
                        writer.writeOpenTag(ContinuousBranchRatesParser.RATE_QUANTILES);
                        writeParameter(clockModel.getParameter("branchRates.quantiles"), -1, writer);
                        writer.writeCloseTag(ContinuousBranchRatesParser.RATE_QUANTILES);
                        writer.writeCloseTag(tag);
                    } else {
                        writer.writeOpenTag(DiscretizedBranchRatesParser.RATE_CATEGORIES);
                        writeParameter(clockModel.getParameter("branchRates.categories"), -1, writer);
                        writer.writeCloseTag(DiscretizedBranchRatesParser.RATE_CATEGORIES);
                        writer.writeCloseTag(tag);
                    }

                    writeMeanRateStatistic(writer, tag, prefix, treePrefix);

                    writeCoefficientOfVariationStatistic(writer, tag, prefix, treePrefix);

                    writeCovarianceStatistic(writer, tag, prefix, treePrefix);

                }

                break;

            case AUTOCORRELATED:
                writer.writeComment("The autocorrelated relaxed clock (Rannala & Yang, 2007)");

                tag = ACLikelihoodParser.AC_LIKELIHOOD;

                attributes = new Attribute[]{
                        new Attribute.Default<String>(XMLParser.ID, prefix + BranchRateModel.BRANCH_RATES),
                        new Attribute.Default<String>("episodic", "false"),
                        new Attribute.Default<String>("logspace", "true"),
                };

                writer.writeOpenTag(tag, attributes);
                writer.writeIDref(TreeModel.TREE_MODEL, treePrefix + TreeModel.TREE_MODEL);

//                    if (!model.isEstimatedRate()) { //TODO move to options or panel select method
//                        Parameter parameter = tree.getParameter(TreeModel.TREE_MODEL + "." + RateEvolutionLikelihood.ROOTRATE);//"treeModel.rootRate"
//                        parameter.isFixed = true;
//                        parameter.initial = model.getRate();
//                    }

                writer.writeOpenTag(RateEvolutionLikelihood.RATES,
                        new Attribute[]{
                                new Attribute.Default<String>(TreeModelParser.ROOT_NODE, "false"),
                                new Attribute.Default<String>(TreeModelParser.INTERNAL_NODES, "true"),
                                new Attribute.Default<String>(TreeModelParser.LEAF_NODES, "true")
                        });
                writer.writeTag(ParameterParser.PARAMETER,
                        new Attribute.Default<String>(XMLParser.ID, treePrefix + TreeModel.TREE_MODEL + "."
                                + TreeModelParser.NODE_RATES), true);
                writer.writeCloseTag(RateEvolutionLikelihood.RATES);

                writer.writeOpenTag(RateEvolutionLikelihood.ROOTRATE,
                        new Attribute[]{
                                new Attribute.Default<String>(TreeModelParser.ROOT_NODE, "true"),
                                new Attribute.Default<String>(TreeModelParser.INTERNAL_NODES, "false"),
                                new Attribute.Default<String>(TreeModelParser.LEAF_NODES, "false")
                        });
                writer.writeTag(ParameterParser.PARAMETER,
                        new Attribute.Default<String>(XMLParser.ID, treePrefix + TreeModel.TREE_MODEL + "."
                                + RateEvolutionLikelihood.ROOTRATE), true);
                writer.writeCloseTag(RateEvolutionLikelihood.ROOTRATE);
                //                writeParameterRef("rates", treePrefix + "treeModel.nodeRates", writer);
                //                writeParameterRef(RateEvolutionLikelihood.ROOTRATE, treePrefix + "treeModel.rootRate", writer);
                writeParameter("variance", "branchRates.var", treeModel, writer);

                writer.writeCloseTag(tag);

//                    if (model.isEstimatedRate()) {//TODO
                writer.writeText("");
                writer.writeOpenTag(CompoundParameterParser.COMPOUND_PARAMETER,
                        new Attribute[]{new Attribute.Default<String>(XMLParser.ID, treePrefix + TreeModel.TREE_MODEL
                                + "." + "allRates")});
                writer.writeIDref(ParameterParser.PARAMETER, treePrefix + TreeModel.TREE_MODEL + "."
                        + TreeModelParser.NODE_RATES);
                writer.writeIDref(ParameterParser.PARAMETER, treePrefix + TreeModel.TREE_MODEL + "."
                        + RateEvolutionLikelihood.ROOTRATE);
                writer.writeCloseTag(CompoundParameterParser.COMPOUND_PARAMETER);
//                    }

                writeMeanRateStatistic(writer, tag, prefix, treePrefix);

                writeCoefficientOfVariationStatistic(writer, tag, prefix, treePrefix);

                writeCovarianceStatistic(writer, tag, prefix, treePrefix);

            case RANDOM_LOCAL_CLOCK: // 1 random local clock CANNOT have different tree models
                writer.writeComment("The random local clock model (Drummond & Suchard, 2010)");

                tag = RandomLocalClockModelParser.LOCAL_BRANCH_RATES;

                writer.writeOpenTag(
                        tag,
                        new Attribute[]{
                                new Attribute.Default<String>(XMLParser.ID, prefix + BranchRateModel.BRANCH_RATES),
                                new Attribute.Default<String>("ratesAreMultipliers", "false")
                        }
                );
                writer.writeIDref(TreeModel.TREE_MODEL, treePrefix + TreeModel.TREE_MODEL);

                writer.writeOpenTag("rates");
                writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default<String>
                                (XMLParser.ID, prefix + ClockType.LOCAL_CLOCK + ".relativeRates")
                        , true);
                writer.writeCloseTag("rates");
                writer.writeOpenTag("rateIndicator");
                writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default<String>
                                (XMLParser.ID, prefix + ClockType.LOCAL_CLOCK + ".changes")
                        , true);
                writer.writeCloseTag("rateIndicator");

                writeParameter("clockRate", "clock.rate", clockModel, writer);

                writer.writeCloseTag(tag);

                writer.writeText("");
                writer.writeOpenTag(
                        SumStatisticParser.SUM_STATISTIC,
                        new Attribute[]{
                                new Attribute.Default<String>(XMLParser.ID, prefix + "rateChanges"),
                                new Attribute.Default<String>("name", prefix + "rateChangeCount"),
                                new Attribute.Default<String>("elementwise", "true"),
                        }
                );
                writer.writeIDref(ParameterParser.PARAMETER, prefix + ClockType.LOCAL_CLOCK + ".changes");
                writer.writeCloseTag(SumStatisticParser.SUM_STATISTIC);

                writeMeanRateStatistic(writer, tag, prefix, treePrefix);

                writeCoefficientOfVariationStatistic(writer, tag, prefix, treePrefix);

                writeCovarianceStatistic(writer, tag, prefix, treePrefix);

                break;

            case FIXED_LOCAL_CLOCK:
                writer.writeComment("The a priori local clock model (Yoder & Yang, 2000)");

                tag = LocalClockModelParser.LOCAL_CLOCK_MODEL;

                writer.writeOpenTag(
                        tag,
                        new Attribute[]{
                                new Attribute.Default<String>(XMLParser.ID, prefix + BranchRateModel.BRANCH_RATES)                        }
                );

                writer.writeIDref(TreeModel.TREE_MODEL, treePrefix + TreeModel.TREE_MODEL);

                writeParameter(LocalClockModelParser.RATE, "clock.rate", clockModel, writer);

                for (Taxa taxonSet : options.taxonSets) {
                    if (options.taxonSetsMono.get(taxonSet)) {
                        String parameterName = taxonSet.getId() + ".rate";
                        writer.writeOpenTag(
                                LocalClockModelParser.CLADE,
                                new Attribute[]{
                                        new Attribute.Default<String>("includeStem", options.taxonSetsIncludeStem.get(taxonSet).toString())
                                }
                        );
                        writeParameter(parameterName, clockModel, writer);
                        writer.writeIDref(TaxaParser.TAXA, taxonSet.getId());
                        writer.writeCloseTag(LocalClockModelParser.CLADE);
                    }
                }

                writer.writeCloseTag(tag);

                writeMeanRateStatistic(writer, tag, prefix, treePrefix);

                writeCoefficientOfVariationStatistic(writer, tag, prefix, treePrefix);

                writeCovarianceStatistic(writer, tag, prefix, treePrefix);
                break;
            default:
                throw new IllegalArgumentException("Unknown clock model");
        }

    }

    private void writeMeanRateStatistic(XMLWriter writer, String tag, String prefix, String treePrefix) {
        writer.writeText("");
        writer.writeOpenTag(
                RateStatisticParser.RATE_STATISTIC,
                new Attribute[]{
                        new Attribute.Default<String>(XMLParser.ID, prefix + "meanRate"),
                        new Attribute.Default<String>("name", prefix + "meanRate"),
                        new Attribute.Default<String>("mode", "mean"),
                        new Attribute.Default<String>("internal", "true"),
                        new Attribute.Default<String>("external", "true")
                }
        );
        writer.writeIDref(TreeModel.TREE_MODEL, treePrefix + TreeModel.TREE_MODEL);
        writer.writeIDref(tag, prefix
                + BranchRateModel.BRANCH_RATES);
        writer.writeCloseTag(RateStatisticParser.RATE_STATISTIC);
    }

    private void writeCovarianceStatistic(XMLWriter writer, String tag, String prefix, String treePrefix) {
        writer.writeText("");
        writer.writeOpenTag(
                RateCovarianceStatisticParser.RATE_COVARIANCE_STATISTIC,
                new Attribute[]{
                        new Attribute.Default<String>(XMLParser.ID, prefix + "covariance"),
                        new Attribute.Default<String>("name", prefix + "covariance")
                }
        );
        writer.writeIDref(TreeModel.TREE_MODEL, treePrefix + TreeModel.TREE_MODEL);
        writer.writeIDref(tag, prefix + BranchRateModel.BRANCH_RATES);
        writer.writeCloseTag(RateCovarianceStatisticParser.RATE_COVARIANCE_STATISTIC);
    }

    private void writeCoefficientOfVariationStatistic(XMLWriter writer, String tag, String prefix, String treePrefix) {
        writer.writeText("");
        writer.writeOpenTag(
                RateStatisticParser.RATE_STATISTIC,
                new Attribute[]{
                        new Attribute.Default<String>(XMLParser.ID, prefix + RateStatisticParser.COEFFICIENT_OF_VARIATION),
                        new Attribute.Default<String>("name", prefix + RateStatisticParser.COEFFICIENT_OF_VARIATION),
                        new Attribute.Default<String>("mode", RateStatisticParser.COEFFICIENT_OF_VARIATION),
                        new Attribute.Default<String>("internal", "true"),
                        new Attribute.Default<String>("external", "true")
                }
        );
        writer.writeIDref(TreeModel.TREE_MODEL, treePrefix + TreeModel.TREE_MODEL);
        writer.writeIDref(tag, prefix + BranchRateModel.BRANCH_RATES);
        writer.writeCloseTag(RateStatisticParser.RATE_STATISTIC);
    }

    /**
     * Write the branch rates model reference.
     *
     * @param model  PartitionClockModel
     * @param writer the writer
     */
    public static void writeBranchRatesModelRef(PartitionClockModel model, XMLWriter writer) {
        String tag = "";
        String id = "";

        switch (model.getClockType()) {
            case STRICT_CLOCK:
                tag = StrictClockBranchRatesParser.STRICT_CLOCK_BRANCH_RATES;
                id = model.getPrefix() + BranchRateModel.BRANCH_RATES;
                break;

            case UNCORRELATED:
                if (model.performModelAveraging()) {
                    tag = MixtureModelBranchRatesParser.MIXTURE_MODEL_BRANCH_RATES;
                } else {
                    tag = model.isContinuousQuantile() ?
                            ContinuousBranchRatesParser.CONTINUOUS_BRANCH_RATES :
                            DiscretizedBranchRatesParser.DISCRETIZED_BRANCH_RATES;
                }
                id = model.getPrefix() + BranchRateModel.BRANCH_RATES;
                break;

            case RANDOM_LOCAL_CLOCK:
                tag = RandomLocalClockModelParser.LOCAL_BRANCH_RATES;
                id = model.getPrefix() + BranchRateModel.BRANCH_RATES;
                break;

            case FIXED_LOCAL_CLOCK:
                tag = LocalClockModelParser.LOCAL_CLOCK_MODEL;
                id = model.getPrefix() + BranchRateModel.BRANCH_RATES;
                break;
            case AUTOCORRELATED:
                tag = ACLikelihoodParser.AC_LIKELIHOOD;
                throw new UnsupportedOperationException("Autocorrelated relaxed clock model not implemented yet");

            default:
                throw new IllegalArgumentException("Unknown clock model");
        }
        writer.writeIDref(tag, id);
    }

    /**
     * Write the allMus for each partition model.
     *
     * @param model  PartitionClockModel
     * @param writer XMLWriter
     */
    public void writeAllMus(PartitionClockModel model, XMLWriter writer) {
        String parameterName = !options.classicOperatorsAndPriors && BeautiOptions.NEW_RELATIVE_RATE_PARAMETERIZATION ? "allNus" : "allMus";

        Parameter allMus = model.getParameter(parameterName);
        if (allMus.getSubParameters().size() > 1) {
            writer.writeComment("Collecting together relative rates for partitions");

            // allMus is global for each gene
            writer.writeOpenTag(CompoundParameterParser.COMPOUND_PARAMETER,
                    new Attribute[]{new Attribute.Default<String>(XMLParser.ID, model.getPrefix() + parameterName)});

            for (Parameter parameter : allMus.getSubParameters()) {
                writer.writeIDref(ParameterParser.PARAMETER, parameter.getName());
            }

            writer.writeCloseTag(CompoundParameterParser.COMPOUND_PARAMETER);
            writer.writeText("");
        }
    }


    public void writeAllClockRateRefs(PartitionClockModel model, XMLWriter writer) {
        writer.writeIDref(ParameterParser.PARAMETER, getClockRateString(model));
    }


    public String getClockRateString(PartitionClockModel model) {
        String prefix = model.getPrefix();

        if (model.performModelAveraging()) {
            return prefix + "meanRate";
        }

        switch (model.getClockType()) {
            case STRICT_CLOCK:
            case RANDOM_LOCAL_CLOCK:
            case FIXED_LOCAL_CLOCK:
                return prefix + "clock.rate";

            case UNCORRELATED:
                switch (model.getClockDistributionType()) {

                    case LOGNORMAL:
                        return prefix + ClockType.UCLD_MEAN;
                    case GAMMA:
                        return prefix + ClockType.UCGD_MEAN;
                    case CAUCHY:
                        throw new UnsupportedOperationException("Uncorrelated Cauchy model not supported yet");
//                        return null;
                    case EXPONENTIAL:
                        return prefix + ClockType.UCED_MEAN;
                }

            case AUTOCORRELATED:
                //TODO
                throw new IllegalArgumentException("Autocorrelated Relaxed Clock, writeAllClockRateRefs(PartitionClockModel model, XMLWriter writer)");
//	        	break;

            default:
                throw new IllegalArgumentException("Unknown clock model");
        }
    }

    public void writeLog(PartitionClockModel model, XMLWriter writer) {
        String prefix = model.getPrefix();

        if (!options.classicOperatorsAndPriors && BeautiOptions.NEW_RELATIVE_RATE_PARAMETERIZATION) {
            Parameter allNus = model.getParameter("allNus");
            if (allNus.getSubParameters().size() > 1) {
                // The mu's are the more relevent parameter and allow comparisons with the old parameterization
                // It would be confusing to log the nus and mus.
//                writer.writeIDref(CompoundParameterParser.COMPOUND_PARAMETER, prefix + "allNus");

                for (Parameter parameter : allNus.getSubParameters()) {
                    String name = parameter.getName();
                    writer.writeIDref(StatisticParser.STATISTIC, name.substring(0, name.lastIndexOf(".")) + ".mu");
                }
            }

        } else {
            Parameter allMus = model.getParameter("allMus");
            if (allMus.getSubParameters().size() > 1) {
                writer.writeIDref(CompoundParameterParser.COMPOUND_PARAMETER, prefix + "allMus");
            }
        }

        switch (model.getClockType()) {
            case STRICT_CLOCK:
            case RANDOM_LOCAL_CLOCK:
                writer.writeIDref(ParameterParser.PARAMETER, prefix + "clock.rate");
                break;

            case FIXED_LOCAL_CLOCK:
                writer.writeIDref(ParameterParser.PARAMETER, prefix + "clock.rate");
                for (Taxa taxonSet : options.taxonSets) {
                    if (options.taxonSetsMono.get(taxonSet)) {
                        String parameterName = taxonSet.getId() + ".rate";
                        writer.writeIDref(ParameterParser.PARAMETER, model.getPrefix() + parameterName);
                    }
                }
                break;

            case UNCORRELATED:

                if (model.performModelAveraging()) {

                    writer.writeIDref(ParameterParser.PARAMETER, prefix + ClockType.UCLD_MEAN);
                    writer.writeIDref(ParameterParser.PARAMETER, prefix + ClockType.UCLD_STDEV);
                    writer.writeIDref(ParameterParser.PARAMETER, prefix + ClockType.UCGD_MEAN);
                    writer.writeIDref(ParameterParser.PARAMETER, prefix + ClockType.UCGD_SHAPE);
                    writer.writeIDref(ParameterParser.PARAMETER, prefix + ClockType.UCED_MEAN);

                    writer.writeIDref(ParameterParser.PARAMETER, "branchRates.distributionIndex");
                    writer.writeIDref(ParameterParser.PARAMETER, "branchRates.quantiles");

                } else {

                    switch (model.getClockDistributionType()) {
                        case LOGNORMAL:
                            writer.writeIDref(ParameterParser.PARAMETER, prefix + ClockType.UCLD_MEAN);
                            writer.writeIDref(ParameterParser.PARAMETER, prefix + ClockType.UCLD_STDEV);
                            break;
                        case GAMMA:
                            writer.writeIDref(ParameterParser.PARAMETER, prefix + ClockType.UCGD_MEAN);
                            writer.writeIDref(ParameterParser.PARAMETER, prefix + ClockType.UCGD_SHAPE);
                            break;
                        case CAUCHY:
                            throw new UnsupportedOperationException("Uncorrelated Couchy model not supported yet");
//                        break;
                        case EXPONENTIAL:
                            writer.writeIDref(ParameterParser.PARAMETER, prefix + ClockType.UCED_MEAN);
                            break;
                    }

                }

            case AUTOCORRELATED:
// TODO
                break;

            default:
                throw new IllegalArgumentException("Unknown clock model");
        }

    }

    public void writeLogStatistic(PartitionClockModel model, XMLWriter writer) {
        String prefix = model.getPrefix();

        switch (model.getClockType()) {
            case STRICT_CLOCK:
                writer.writeIDref(RateStatisticParser.RATE_STATISTIC, prefix + "meanRate");
                break;

            case UNCORRELATED:
            case FIXED_LOCAL_CLOCK:
                writer.writeIDref(RateStatisticParser.RATE_STATISTIC, prefix + "meanRate");
                writer.writeIDref(RateStatisticParser.RATE_STATISTIC, prefix + RateStatisticParser.COEFFICIENT_OF_VARIATION);
                writer.writeIDref(RateCovarianceStatisticParser.RATE_COVARIANCE_STATISTIC, prefix + "covariance");
                break;

            case AUTOCORRELATED:
                writer.writeIDref(RateStatisticParser.RATE_STATISTIC, prefix + "meanRate");
                writer.writeIDref(RateStatisticParser.RATE_STATISTIC, prefix + RateStatisticParser.COEFFICIENT_OF_VARIATION);
                writer.writeIDref(RateCovarianceStatisticParser.RATE_COVARIANCE_STATISTIC, prefix + "covariance");
                writer.writeIDref(ParameterParser.PARAMETER, prefix + "branchRates.var");
                writer.writeIDref(ParameterParser.PARAMETER, prefix + "treeModel.rootRate");
                break;

            case RANDOM_LOCAL_CLOCK:
                writer.writeIDref(SumStatisticParser.SUM_STATISTIC, model.getPrefix() + "rateChanges");
                writer.writeIDref(RateStatisticParser.RATE_STATISTIC, prefix + "meanRate");
                writer.writeIDref(RateStatisticParser.RATE_STATISTIC, prefix + RateStatisticParser.COEFFICIENT_OF_VARIATION);
                writer.writeIDref(RateCovarianceStatisticParser.RATE_COVARIANCE_STATISTIC, prefix + "covariance");
                break;

            default:
                throw new IllegalArgumentException("Unknown clock model");
        }

    }

    public void writeClockLikelihoodReferences(XMLWriter writer) {
        for (PartitionClockModel clockModel : options.getPartitionClockModels()) { // Each PD has one TreeLikelihood
            writeBranchRatesModelRef(clockModel, writer);
        }
    }

}
