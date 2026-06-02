/*
 * ClockModelGenerator.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
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
 *
 */

package dr.app.beauti.generator;

import dr.app.beauti.components.ComponentFactory;
import dr.app.beauti.options.*;
import dr.app.beauti.types.ClockType;
import dr.app.beauti.util.XMLWriter;
import dr.evolution.util.Taxa;
import dr.evomodel.branchratemodel.ArbitraryBranchRates;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.branchratemodel.BranchSpecificFixedEffects;
import dr.evomodel.tree.DefaultTreeModel;
import dr.evomodel.treedatalikelihood.discrete.LocationGradient;
import dr.evomodel.treedatalikelihood.discrete.ScaleGradient;
import dr.evomodelxml.branchmodel.BranchSpecificBranchModelParser;
import dr.evomodelxml.branchratemodel.*;
import dr.evomodelxml.continuous.hmc.LocationScaleGradientParser;
import dr.evomodelxml.tree.CTMCScalePriorParser;
import dr.evomodelxml.tree.RateCovarianceStatisticParser;
import dr.evomodelxml.tree.RateStatisticParser;
import dr.evomodelxml.tree.TreeModelParser;
import dr.evomodelxml.treedatalikelihood.TreeDataLikelihoodParser;
import dr.evoxml.TaxaParser;
import dr.inference.distribution.DistributionLikelihood;
import dr.inference.distribution.ExponentialDistributionModel;
import dr.inference.distribution.GammaDistributionModel;
import dr.inference.model.ParameterParser;
import dr.inference.model.StatisticParser;
import dr.inferencexml.SignTransformParser;
import dr.inferencexml.distribution.*;
import dr.inferencexml.distribution.shrinkage.BayesianBridgeDistributionModelParser;
import dr.inferencexml.hmc.CompoundGradientParser;
import dr.inferencexml.hmc.HessianWrapperParser;
import dr.inferencexml.hmc.JointGradientParser;
import dr.inferencexml.model.CompoundParameterParser;
import dr.inferencexml.model.SumStatisticParser;
import dr.oldevomodel.clock.RateEvolutionLikelihood;
import dr.oldevomodelxml.clock.ACLikelihoodParser;
import dr.util.Attribute;
import dr.xml.XMLParser;

import static dr.inference.model.ParameterParser.PARAMETER;
import static dr.inferencexml.distribution.PriorParsers.*;
import static dr.inferencexml.distribution.shrinkage.BayesianBridgeLikelihoodParser.*;

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
                    writer.writeIDref(DefaultTreeModel.TREE_MODEL, treePrefix + DefaultTreeModel.TREE_MODEL);

                    writer.writeOpenTag(DistributionLikelihoodParser.DISTRIBUTION);
                    writer.writeOpenTag(LogNormalDistributionModelParser.LOGNORMAL_DISTRIBUTION_MODEL,
                            new Attribute.Default<String>(LogNormalDistributionModelParser.MEAN_IN_REAL_SPACE, "true"));
                    writeParameter("mean", ClockType.UCLD_MEAN, clockModel, writer);
                    writeParameter("stdev", ClockType.UCLD_STDEV, clockModel, writer);
                    writer.writeCloseTag(LogNormalDistributionModelParser.LOGNORMAL_DISTRIBUTION_MODEL);
                    writer.writeCloseTag(DistributionLikelihoodParser.DISTRIBUTION);

                    writer.writeOpenTag(DistributionLikelihoodParser.DISTRIBUTION);
                    writer.writeOpenTag(GammaDistributionModel.GAMMA_DISTRIBUTION_MODEL);
                    writeParameter("mean", ClockType.UCGD_MEAN, clockModel, writer);
                    writeParameter("shape", ClockType.UCGD_SHAPE, clockModel, writer);
                    writer.writeCloseTag(GammaDistributionModel.GAMMA_DISTRIBUTION_MODEL);
                    writer.writeCloseTag(DistributionLikelihoodParser.DISTRIBUTION);

                    writer.writeOpenTag(DistributionLikelihoodParser.DISTRIBUTION);
                    writer.writeOpenTag(ExponentialDistributionModel.EXPONENTIAL_DISTRIBUTION_MODEL);
                    writeParameter("mean", ClockType.UCED_MEAN, clockModel, writer);
                    writer.writeCloseTag(ExponentialDistributionModel.EXPONENTIAL_DISTRIBUTION_MODEL);
                    writer.writeCloseTag(DistributionLikelihoodParser.DISTRIBUTION);

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
                    writer.writeIDref(DefaultTreeModel.TREE_MODEL, treePrefix + DefaultTreeModel.TREE_MODEL);

                    writer.writeOpenTag(DistributionLikelihoodParser.DISTRIBUTION);

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

                    writer.writeCloseTag(DistributionLikelihoodParser.DISTRIBUTION);

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
            case HMC_CLOCK:
                tag = ArbitraryBranchRatesParser.ARBITRARY_BRANCH_RATES;

                writer.writeComment("The Hamiltonian Monte Carlo relaxed clock");
                attributes = new Attribute[] {
                        new Attribute.Default<>(XMLParser.ID,
                                prefix  + BranchRateModel.BRANCH_RATES),
                        new Attribute.Default<>(ArbitraryBranchRatesParser.CENTER_AT_ONE, false)
                };
                writer.writeOpenTag(tag, attributes);
                // tree
                writer.writeIDref(DefaultTreeModel.TREE_MODEL, treePrefix + DefaultTreeModel.TREE_MODEL);

                writer.writeOpenTag(ArbitraryBranchRatesParser.RATES);
                writeParameter(clockModel.getParameter(ClockType.HMC_CLOCK_BRANCH_RATES), -1, writer);
                writer.writeCloseTag(ArbitraryBranchRatesParser.RATES);
                writer.writeOpenTag(ArbitraryBranchRatesParser.LOCATION);
                writeParameter(clockModel.getParameter(ClockType.HMC_CLOCK_LOCATION), -1, writer);
                writer.writeCloseTag(ArbitraryBranchRatesParser.LOCATION);
                writer.writeOpenTag(ArbitraryBranchRatesParser.SCALE);
                writeParameter(clockModel.getParameter(ClockType.HMCLN_SCALE), -1, writer);
                writer.writeCloseTag(ArbitraryBranchRatesParser.SCALE);
                writer.writeCloseTag(tag);

                //rates prior
                writer.writeOpenTag(DistributionLikelihood.DISTRIBUTION_LIKELIHOOD,
                        new Attribute.Default<>(XMLParser.ID,
                                prefix + BranchSpecificFixedEffects.RATES_PRIOR));

                writeParameterRef(MixedDistributionLikelihoodParser.DATA, prefix + ClockType.HMC_CLOCK_BRANCH_RATES, writer);

                writer.writeOpenTag(DistributionLikelihoodParser.DISTRIBUTION);

                switch (clockModel.getClockDistributionType()) {
                    case LOGNORMAL:
                        writer.writeOpenTag(LogNormalDistributionModelParser.LOGNORMAL_DISTRIBUTION_MODEL,
                                new Attribute.Default<String>(LogNormalDistributionModelParser.MEAN_IN_REAL_SPACE, "true"));

                        writer.writeOpenTag("mean");
                        writeParameter(null, 1, 1.0, 0.0, Double.NaN, writer);
                        writer.writeCloseTag("mean");
                        writer.writeOpenTag("stdev");
                        writeParameter(null, 1, 0.1, 0.0, Double.NaN, writer);
                        writer.writeCloseTag("stdev");

                        writer.writeCloseTag(LogNormalDistributionModelParser.LOGNORMAL_DISTRIBUTION_MODEL);
                        break;
                    default:
                        throw new UnsupportedOperationException("Only lognormal is supported for HMC relaxed clock");
                }
                writer.writeCloseTag(DistributionLikelihoodParser.DISTRIBUTION);

                writer.writeCloseTag(DistributionLikelihood.DISTRIBUTION_LIKELIHOOD);

                writeMeanRateStatistic(writer, tag, prefix, treePrefix);

                writeCoefficientOfVariationStatistic(writer, tag, prefix, treePrefix);

                writeCovarianceStatistic(writer, tag, prefix, treePrefix);

                boolean generateRatesGradient = false;
                boolean generateScaleGradient = false;

                for (Operator operator : options.selectOperators()) {
                    if (operator.getName().equals(ClockType.HMC_CLOCK_RATES_DESCRIPTION) && operator.isUsed()) {
                        generateRatesGradient = true;
                    }
                    if (operator.getName().equals(ClockType.HMC_CLOCK_LOCATION_SCALE_DESCRIPTION) && operator.isUsed()) {
                        generateScaleGradient = true;
                    }
                }

                if (generateRatesGradient) {

                    //scale prior
                    writer.writeOpenTag(DistributionLikelihood.DISTRIBUTION_LIKELIHOOD,
                            new Attribute.Default<>(XMLParser.ID,
                                    prefix + BranchSpecificFixedEffects.SCALE_PRIOR));
                    writeParameterRef(MixedDistributionLikelihoodParser.DATA, prefix + ClockType.HMCLN_SCALE, writer);
                    writer.writeOpenTag(DistributionLikelihoodParser.DISTRIBUTION);
                    writer.writeOpenTag(ExponentialDistributionModel.EXPONENTIAL_DISTRIBUTION_MODEL);
                    writer.writeOpenTag(ExponentialDistributionModelParser.MEAN);
                    writeParameter(null, 1, 1.0, 0.0, Double.NaN, writer);
                    writer.writeCloseTag(ExponentialDistributionModelParser.MEAN);
                    writer.writeCloseTag(ExponentialDistributionModel.EXPONENTIAL_DISTRIBUTION_MODEL);
                    writer.writeCloseTag(DistributionLikelihoodParser.DISTRIBUTION);
                    writer.writeCloseTag(DistributionLikelihood.DISTRIBUTION_LIKELIHOOD);

                    //compound parameter
                    writer.writeOpenTag(CompoundParameterParser.COMPOUND_PARAMETER, new Attribute.Default<>(XMLParser.ID, prefix + LocationScaleGradientParser.LOCATION_SCALE));
                    writer.writeIDref(ParameterParser.PARAMETER, prefix + ClockType.HMC_CLOCK_LOCATION);
                    writer.writeIDref(ParameterParser.PARAMETER, prefix + ClockType.HMCLN_SCALE);
                    writer.writeCloseTag(CompoundParameterParser.COMPOUND_PARAMETER);

                    //CTMC scale prior
                    writer.writeOpenTag(CTMCScalePriorParser.MODEL_NAME, new Attribute.Default<>(XMLParser.ID, prefix + BranchSpecificFixedEffects.LOCATION_PRIOR));
                    writer.writeOpenTag(CTMCScalePriorParser.SCALEPARAMETER);
                    writer.writeIDref(ParameterParser.PARAMETER, prefix + ClockType.HMC_CLOCK_LOCATION);
                    writer.writeCloseTag(CTMCScalePriorParser.SCALEPARAMETER);
                    writer.writeIDref(DefaultTreeModel.TREE_MODEL, treePrefix + DefaultTreeModel.TREE_MODEL);
                    writer.writeCloseTag(CTMCScalePriorParser.MODEL_NAME);

                }

                if (generateScaleGradient){
                    //location gradient
                    writer.writeOpenTag(LocationScaleGradientParser.NAME, new Attribute[]{
                            new Attribute.Default<>(XMLParser.ID, prefix + LocationGradient.LOCATION_GRADIENT),
                            new Attribute.Default<>("traitName", "Sequence"),
                            new Attribute.Default<>(LocationScaleGradientParser.USE_HESSIAN, "false")
                    });
                    writer.writeIDref(TreeDataLikelihoodParser.TREE_DATA_LIKELIHOOD, treePrefix + "treeLikelihood");
                    writer.writeOpenTag(LocationScaleGradientParser.LOCATION);
                    writer.writeIDref(ParameterParser.PARAMETER, prefix + ClockType.HMC_CLOCK_LOCATION);
                    writer.writeCloseTag(LocationScaleGradientParser.LOCATION);
                    writer.writeCloseTag(LocationScaleGradientParser.NAME);

                    //scale gradient
                    writer.writeOpenTag(LocationScaleGradientParser.NAME, new Attribute[]{
                            new Attribute.Default<>(XMLParser.ID, prefix + ScaleGradient.SCALE_GRADIENT),
                            new Attribute.Default<>("traitName", "Sequence"),
                            new Attribute.Default<>(LocationScaleGradientParser.USE_HESSIAN, "false")
                    });
                    writer.writeIDref(TreeDataLikelihoodParser.TREE_DATA_LIKELIHOOD, treePrefix + "treeLikelihood");
                    writer.writeOpenTag(LocationScaleGradientParser.SCALE);
                    writer.writeIDref(ParameterParser.PARAMETER, prefix + ClockType.HMCLN_SCALE);
                    writer.writeCloseTag(LocationScaleGradientParser.SCALE);
                    writer.writeCloseTag(LocationScaleGradientParser.NAME);

                    //location scale (compound) gradient
                    writer.writeOpenTag(CompoundGradientParser.COMPOUND_GRADIENT, new Attribute.Default<>(XMLParser.ID, prefix + LocationScaleGradientParser.NAME));
                    writer.writeIDref(LocationScaleGradientParser.NAME, prefix + LocationGradient.LOCATION_GRADIENT);
                    writer.writeIDref(LocationScaleGradientParser.NAME, prefix + ScaleGradient.SCALE_GRADIENT);
                    writer.writeCloseTag(CompoundGradientParser.COMPOUND_GRADIENT);

                    //location scale (compound) prior gradient
                    writer.writeOpenTag(CompoundGradientParser.COMPOUND_GRADIENT, new Attribute.Default<>(XMLParser.ID, prefix + LocationScaleGradientParser.LOCATION_SCALE_PRIOR_GRADIENT));
                    writer.writeIDref(CTMCScalePriorParser.MODEL_NAME, prefix + BranchSpecificFixedEffects.LOCATION_PRIOR);
                    writer.writeOpenTag(HessianWrapperParser.NAME);
                    writer.writeIDref(DistributionLikelihood.DISTRIBUTION_LIKELIHOOD, BranchSpecificFixedEffects.SCALE_PRIOR);
                    writer.writeIDref(ParameterParser.PARAMETER, prefix + ClockType.HMCLN_SCALE);
                    writer.writeCloseTag(HessianWrapperParser.NAME);
                    writer.writeCloseTag(CompoundGradientParser.COMPOUND_GRADIENT);

                    //location scale joint gradient
                    writer.writeOpenTag(JointGradientParser.JOINT_GRADIENT, new Attribute.Default<>(XMLParser.ID, prefix + LocationScaleGradientParser.LOCATION_SCALE_JOINT_GRADIENT));
                    writer.writeIDref(CompoundGradientParser.COMPOUND_GRADIENT, prefix + LocationScaleGradientParser.LOCATION_SCALE_PRIOR_GRADIENT);
                    writer.writeIDref(CompoundGradientParser.COMPOUND_GRADIENT, prefix + LocationScaleGradientParser.NAME);
                    writer.writeCloseTag(JointGradientParser.JOINT_GRADIENT);

                }

                break;

            case SHRINKAGE_LOCAL_CLOCK:
                writer.writeComment("The shrinkage local clock");

                tag = ArbitraryBranchRatesParser.ARBITRARY_BRANCH_RATES;

                attributes = new Attribute[] {
                        new Attribute.Default<>(XMLParser.ID,
                                prefix + "substBranchRates"),
                        new Attribute.Default<>(ArbitraryBranchRatesParser.CENTER_AT_ONE, false),
                        new Attribute.Default<>(ArbitraryBranchRatesParser.RANDOMIZE_RATES, true),
                        new Attribute.Default<>(ArbitraryBranchRatesParser.RANDOM_SCALE, "0.1")
                };
                writer.writeOpenTag(tag, attributes);
                // tree
                writer.writeIDref(DefaultTreeModel.TREE_MODEL, treePrefix + DefaultTreeModel.TREE_MODEL);

                writer.writeOpenTag(ArbitraryBranchRatesParser.RATES);
                writeParameter(clockModel.getParameter("substBranchRates.rates"), -1, writer);
                writer.writeCloseTag(ArbitraryBranchRatesParser.RATES);
                writer.writeCloseTag(tag);

                writer.writeOpenTag(BayesianBridgeDistributionModelParser.BAYESIAN_BRIDGE_DISTRIBUTION,
                        new Attribute.Default<>(XMLParser.ID,prefix  + "bbDistribution"));
                writer.writeOpenTag(BayesianBridgeDistributionModelParser.DIMENSION);
                writer.writeIDref(PARAMETER, prefix + "substBranchRates.rates");
                writer.writeCloseTag(BayesianBridgeDistributionModelParser.DIMENSION);
                writer.writeOpenTag(SLAB_WIDTH);
                writeParameter(prefix + SLAB_WIDTH, -1, 2.0, Double.NaN, Double.NaN, writer);
                writer.writeCloseTag(SLAB_WIDTH);
                writer.writeOpenTag(GLOBAL_SCALE);
                writeParameter(prefix + GLOBAL_SCALE, -1, 1.0, Double.NaN, Double.NaN, writer);
                writer.writeCloseTag(GLOBAL_SCALE);
                writer.writeOpenTag(EXPONENT);
                writeParameter(prefix + EXPONENT, -1, 0.25, Double.NaN, Double.NaN, writer);
                writer.writeCloseTag(EXPONENT);
                writer.writeOpenTag(LOCAL_SCALE);
                writeParameter(prefix + LOCAL_SCALE, -1, 10.0, Double.NaN, Double.NaN, writer);
                writer.writeCloseTag(LOCAL_SCALE);
                writer.writeCloseTag(BayesianBridgeDistributionModelParser.BAYESIAN_BRIDGE_DISTRIBUTION);

                writer.writeOpenTag(AutoCorrelatedBranchRatesDistributionParser.AUTO_CORRELATED_RATES,
                        new Attribute[] {
                                new Attribute.Default<>(XMLParser.ID,
                                        prefix  + "substBranchRatesPrior"),
                                new Attribute.Default<>(AutoCorrelatedBranchRatesDistributionParser.SCALING, "none"),
                                new Attribute.Default<>(AutoCorrelatedBranchRatesDistributionParser.LOG, true),
                                new Attribute.Default<>(AutoCorrelatedBranchRatesDistributionParser.OPERATE_ON_INCREMENTS, true),
                        });
                writer.writeIDref(ArbitraryBranchRatesParser.ARBITRARY_BRANCH_RATES, prefix  + "substBranchRates");
                writer.writeIDref(BayesianBridgeDistributionModelParser.BAYESIAN_BRIDGE_DISTRIBUTION, prefix  + "bbDistribution");
                writer.writeCloseTag(AutoCorrelatedBranchRatesDistributionParser.AUTO_CORRELATED_RATES);

                writer.writeOpenTag(GAMMA_PRIOR,
                        new Attribute[] {
                                new Attribute.Default<>(XMLParser.ID,
                                        prefix  + "globalScalePrior"),
                                new Attribute.Default<>(SHAPE, "1"),
                                new Attribute.Default<>(SCALE, "2")
                        });
                writer.writeIDref(PARAMETER, prefix + GLOBAL_SCALE);
                writer.writeCloseTag(GAMMA_PRIOR);

                writer.writeOpenTag(ScaledByTreeTimeBranchRateModelParser.TREE_TIME_BRANCH_RATES,
                        new Attribute.Default<>(XMLParser.ID, prefix  + BranchRateModel.BRANCH_RATES));
                writer.writeIDref(ArbitraryBranchRatesParser.ARBITRARY_BRANCH_RATES, prefix  + "substBranchRates");
                writer.writeIDref(DefaultTreeModel.TREE_MODEL, prefix  + DefaultTreeModel.TREE_MODEL);
                writeParameter(clockModel.getParameter(ClockType.SHRINKAGE_CLOCK_LOCATION), -1, writer);
                writer.writeCloseTag(ScaledByTreeTimeBranchRateModelParser.TREE_TIME_BRANCH_RATES);

//                writeMeanRateStatistic(writer, tag, prefix, treePrefix);
//
//                writeCoefficientOfVariationStatistic(writer, tag, prefix, treePrefix);
//
//                writeCovarianceStatistic(writer, tag, prefix, treePrefix);

                writer.writeComment("gradient of subst branch rates prior wrt increments");
                writer.writeOpenTag(AutoCorrelatedGradientWrtIncrementsParser.GRADIENT,
                        new Attribute.Default<>(XMLParser.ID, prefix  + "incrementGradient"));
                writer.writeIDref(AutoCorrelatedBranchRatesDistributionParser.AUTO_CORRELATED_RATES, prefix + "substBranchRatesPrior");
                writer.writeCloseTag(AutoCorrelatedGradientWrtIncrementsParser.GRADIENT);

//                writer.writeComment("log branch rates for tree");
//                writer.writeOpenTag(TransformedTreeTraitParser.NAME,
//                        new Attribute.Default<>(XMLParser.ID, prefix  + "logBranchRates"));
//                writer.writeTag("signTransform", true);
//                writer.writeIDref(ArbitraryBranchRatesParser.RATES, prefix  + "branchRates");
//                writer.writeCloseTag(TransformedTreeTraitParser.NAME);

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
                writer.writeIDref(DefaultTreeModel.TREE_MODEL, treePrefix + DefaultTreeModel.TREE_MODEL);

//                    if (!model.isEstimatedRate()) { //TODO move to options or panel select method
//                        Parameter parameter = tree.getParameter(DefaultTreeModel.TREE_MODEL + "." + RateEvolutionLikelihood.ROOTRATE);//"treeModel.rootRate"
//                        parameter.isFixed = true;
//                        parameter.initial = model.getRate();
//                    }

                writer.writeOpenTag(RateEvolutionLikelihood.RATES,
                        new Attribute[]{
                                new Attribute.Default<String>(TreeModelParser.ROOT_NODE, "false"),
                                new Attribute.Default<String>(TreeModelParser.INTERNAL_NODES, "true"),
                                new Attribute.Default<String>(TreeModelParser.LEAF_NODES, "true")
                        });
                writer.writeTag(PARAMETER,
                        new Attribute.Default<String>(XMLParser.ID, treePrefix + DefaultTreeModel.TREE_MODEL + "."
                                + TreeModelParser.NODE_RATES), true);
                writer.writeCloseTag(RateEvolutionLikelihood.RATES);

                writer.writeOpenTag(RateEvolutionLikelihood.ROOTRATE,
                        new Attribute[]{
                                new Attribute.Default<String>(TreeModelParser.ROOT_NODE, "true"),
                                new Attribute.Default<String>(TreeModelParser.INTERNAL_NODES, "false"),
                                new Attribute.Default<String>(TreeModelParser.LEAF_NODES, "false")
                        });
                writer.writeTag(PARAMETER,
                        new Attribute.Default<String>(XMLParser.ID, treePrefix + DefaultTreeModel.TREE_MODEL + "."
                                + RateEvolutionLikelihood.ROOTRATE), true);
                writer.writeCloseTag(RateEvolutionLikelihood.ROOTRATE);
                //                writeParameterRef("rates", treePrefix + "treeModel.nodeRates", writer);
                //                writeParameterRef(RateEvolutionLikelihood.ROOTRATE, treePrefix + "treeModel.rootRate", writer);
                writeParameter("variance", "branchRates.var", treeModel, writer);

                writer.writeCloseTag(tag);

//                    if (model.isEstimatedRate()) {//TODO
                writer.writeText("");
                writer.writeOpenTag(CompoundParameterParser.COMPOUND_PARAMETER,
                        new Attribute[]{new Attribute.Default<String>(XMLParser.ID, treePrefix + DefaultTreeModel.TREE_MODEL
                                + "." + "allRates")});
                writer.writeIDref(PARAMETER, treePrefix + DefaultTreeModel.TREE_MODEL + "."
                        + TreeModelParser.NODE_RATES);
                writer.writeIDref(PARAMETER, treePrefix + DefaultTreeModel.TREE_MODEL + "."
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
                writer.writeIDref(DefaultTreeModel.TREE_MODEL, treePrefix + DefaultTreeModel.TREE_MODEL);

                writer.writeOpenTag("rates");
                writer.writeTag(PARAMETER, new Attribute.Default<String>
                                (XMLParser.ID, prefix + ClockType.LOCAL_CLOCK + ".relativeRates")
                        , true);
                writer.writeCloseTag("rates");
                writer.writeOpenTag("rateIndicator");
                writer.writeTag(PARAMETER, new Attribute.Default<String>
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
                writer.writeIDref(PARAMETER, prefix + ClockType.LOCAL_CLOCK + ".changes");
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
                                new Attribute.Default<String>(XMLParser.ID, prefix + BranchRateModel.BRANCH_RATES)
                        }
                );

                writer.writeIDref(DefaultTreeModel.TREE_MODEL, treePrefix + DefaultTreeModel.TREE_MODEL);

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

            case MIXED_EFFECTS_CLOCK:
                writer.writeComment("The mixed effects clock model (Bletsa et al., Virus Evol., 2019)");

                tag = CompoundParameterParser.COMPOUND_PARAMETER;

                //first write the CompoundParameter XML bit
                writer.writeOpenTag(tag,
                        new Attribute[]{
                                new Attribute.Default<String>(XMLParser.ID, prefix + BranchSpecificFixedEffectsParser.FIXED_EFFECTS)
                        }
                );

                writer.writeTag(PARAMETER, new Attribute[]{
                        new Attribute.Default<String>(XMLParser.ID, prefix + BranchSpecificFixedEffectsParser.INTERCEPT),
                        new Attribute.Default<String>(ParameterParser.VALUE, "-0.01")}, true);
                int parameterNumber = 1;
                for (Taxa taxonSet : options.taxonSets) {
                    if (options.taxonSetsMono.get(taxonSet)) {
                        writer.writeTag(PARAMETER, new Attribute[]{
                                new Attribute.Default<String>(XMLParser.ID, prefix + BranchSpecificFixedEffectsParser.COEFFICIENT + parameterNumber),
                                new Attribute.Default<String>(ParameterParser.VALUE, "0.01")}, true);
                    }
                    parameterNumber++;
                }
                writer.writeCloseTag(tag);

                //continue with the fixedEffects XML block
                tag = BranchSpecificFixedEffectsParser.FIXED_EFFECTS;
                writer.writeOpenTag(tag,
                        new Attribute[]{
                                new Attribute.Default<String>(XMLParser.ID, prefix + BranchSpecificFixedEffectsParser.FIXED_EFFECTS_MODEL)
                        }
                );

                writer.writeIDref(DefaultTreeModel.TREE_MODEL, treePrefix + DefaultTreeModel.TREE_MODEL);
                writer.writeIDref(CompoundParameterParser.COMPOUND_PARAMETER, prefix + BranchSpecificFixedEffectsParser.FIXED_EFFECTS);

                int counter = 1;
                for (Taxa taxonSet : options.taxonSets) {
                    writer.writeOpenTag(BranchSpecificFixedEffectsParser.CATEGORY,
                            new Attribute[]{
                                    new Attribute.Default<String>(XMLParser.ID, "" + counter)
                            }
                    );
                    writer.writeOpenTag(BranchSpecificBranchModelParser.CLADE,
                            new Attribute[]{
                                    new Attribute.Default<String>(BranchSpecificFixedEffectsParser.CATEGORY, "" + (counter+1)),
                                    new Attribute.Default<String>(LocalClockModelParser.INCLUDE_STEM, Boolean.toString(options.taxonSetsIncludeStem.get(taxonSet))),
                                    new Attribute.Default<String>(LocalClockModelParser.EXCLUDE_CLADE, "false")
                            }
                    );
                    writer.writeIDref(TaxaParser.TAXA, taxonSet.getId());
                    writer.writeCloseTag(BranchSpecificBranchModelParser.CLADE);
                    writer.writeCloseTag(BranchSpecificFixedEffectsParser.CATEGORY);
                    counter++;
                }
                writer.writeTag(SignTransformParser.NAME, true);

                writer.writeCloseTag(tag);

                //and then the arbitraryBranchRates
                tag = ArbitraryBranchRatesParser.ARBITRARY_BRANCH_RATES;
                writer.writeOpenTag(tag,
                        new Attribute[]{
                                new Attribute.Default<String>(XMLParser.ID, prefix + BranchRateModel.BRANCH_RATES),
                                new Attribute.Default<String>(ArbitraryBranchRatesParser.RECIPROCAL, "false"),
                                new Attribute.Default<String>(ArbitraryBranchRatesParser.EXP, "false"),
                                new Attribute.Default<String>(ArbitraryBranchRatesParser.CENTER_AT_ONE, "false")
                        }
                );

                writer.writeOpenTag(ArbitraryBranchRatesParser.RATES);
                writer.writeTag(PARAMETER,
                        new Attribute[]{
                                new Attribute.Default<String>(XMLParser.ID, prefix + ClockType.ME_CLOCK_LOCATION),
                                new Attribute.Default<Double>(ParameterParser.VALUE, 0.2),
                                new Attribute.Default<Double>(ParameterParser.LOWER, 0.0)
                        }, true);
                writer.writeCloseTag(ArbitraryBranchRatesParser.RATES);

                writer.writeOpenTag(ArbitraryBranchRatesParser.SCALE);
                writer.writeTag(PARAMETER,
                        new Attribute[]{
                                new Attribute.Default<String>(XMLParser.ID, prefix + ClockType.ME_CLOCK_SCALE),
                                new Attribute.Default<Double>(ParameterParser.VALUE, 0.15),
                                new Attribute.Default<Double>(ParameterParser.LOWER, 0.0)
                        }, true);
                writer.writeCloseTag(ArbitraryBranchRatesParser.SCALE);

                writer.writeOpenTag(ArbitraryBranchRatesParser.LOCATION);
                writer.writeIDref(BranchSpecificFixedEffectsParser.FIXED_EFFECTS, BranchSpecificFixedEffectsParser.FIXED_EFFECTS_MODEL);
                writer.writeCloseTag(ArbitraryBranchRatesParser.LOCATION);

                writer.writeIDref(DefaultTreeModel.TREE_MODEL, treePrefix + DefaultTreeModel.TREE_MODEL);

                writer.writeCloseTag(tag);

                tag = DistributionLikelihood.DISTRIBUTION_LIKELIHOOD;

                // branch rates prior
                writer.writeOpenTag(tag,
                        new Attribute[]{
                                new Attribute.Default<String>(XMLParser.ID, BranchSpecificFixedEffects.RATES_PRIOR)
                        }
                        );
                writer.writeOpenTag(DistributionLikelihoodParser.DATA);
                writer.writeIDref(PARAMETER, ClockType.ME_CLOCK_LOCATION);
                writer.writeCloseTag(DistributionLikelihoodParser.DATA);
                writer.writeOpenTag(DistributionLikelihoodParser.DISTRIBUTION);
                writer.writeOpenTag(LogNormalDistributionModelParser.LOGNORMAL_DISTRIBUTION_MODEL,
                        new Attribute[]{
                                new Attribute.Default<String>(LogNormalDistributionModelParser.MEAN_IN_REAL_SPACE, "true")
                        }
                        );
                writer.writeOpenTag(LogNormalDistributionModelParser.MEAN);
                writer.writeTag(PARAMETER,
                        new Attribute[]{
                                new Attribute.Default<String>(ParameterParser.VALUE, "1.0"),
                                new Attribute.Default<String>(ParameterParser.LOWER, "0.0")
                        }, true);
                writer.writeCloseTag(LogNormalDistributionModelParser.MEAN);
                writer.writeOpenTag(LogNormalDistributionModelParser.STDEV);
                writer.writeTag(PARAMETER,
                        new Attribute[]{
                                new Attribute.Default<String>(ParameterParser.VALUE, "1.0"),
                                new Attribute.Default<String>(ParameterParser.LOWER, "0.0")
                        }, true);
                writer.writeCloseTag(LogNormalDistributionModelParser.STDEV);
                writer.writeCloseTag(LogNormalDistributionModelParser.LOGNORMAL_DISTRIBUTION_MODEL);
                writer.writeCloseTag(DistributionLikelihoodParser.DISTRIBUTION);
                writer.writeCloseTag(tag);

                // branch rates scale prior
                writer.writeOpenTag(tag,
                        new Attribute[]{
                                new Attribute.Default<String>(XMLParser.ID, BranchSpecificFixedEffects.SCALE_PRIOR)
                        }
                );
                writer.writeOpenTag(DistributionLikelihoodParser.DATA);
                writer.writeIDref(PARAMETER, ClockType.ME_CLOCK_SCALE);
                writer.writeCloseTag(DistributionLikelihoodParser.DATA);
                writer.writeOpenTag(DistributionLikelihoodParser.DISTRIBUTION);
                writer.writeOpenTag(ExponentialDistributionModel.EXPONENTIAL_DISTRIBUTION_MODEL);
                writer.writeOpenTag(ExponentialDistributionModelParser.MEAN);
                writer.writeTag(PARAMETER,
                        new Attribute[]{
                                new Attribute.Default<String>(ParameterParser.VALUE, "1.0"),
                                new Attribute.Default<String>(ParameterParser.LOWER, "0.0")
                        }, true);
                writer.writeCloseTag(ExponentialDistributionModelParser.MEAN);
                writer.writeCloseTag(ExponentialDistributionModel.EXPONENTIAL_DISTRIBUTION_MODEL);
                writer.writeCloseTag(DistributionLikelihoodParser.DISTRIBUTION);
                writer.writeCloseTag(tag);

                // intercept prior
                writer.writeOpenTag(tag,
                        new Attribute[]{
                                new Attribute.Default<String>(XMLParser.ID, BranchSpecificFixedEffects.INTERCEPT_PRIOR)
                        }
                );
                writer.writeOpenTag(DistributionLikelihoodParser.DATA);
                writer.writeIDref(PARAMETER, BranchSpecificFixedEffectsParser.INTERCEPT);
                writer.writeCloseTag(DistributionLikelihoodParser.DATA);
                writer.writeOpenTag(DistributionLikelihoodParser.DISTRIBUTION);
                writer.writeOpenTag(NormalDistributionModelParser.NORMAL_DISTRIBUTION_MODEL);
                writer.writeOpenTag(LogNormalDistributionModelParser.MEAN);
                writer.writeTag(PARAMETER,
                        new Attribute[]{
                                new Attribute.Default<String>(ParameterParser.VALUE, "0.0"),
                                new Attribute.Default<String>(ParameterParser.LOWER, "0.0")
                        }, true);
                writer.writeCloseTag(LogNormalDistributionModelParser.MEAN);
                writer.writeOpenTag(LogNormalDistributionModelParser.STDEV);
                writer.writeTag(PARAMETER,
                        new Attribute[]{
                                new Attribute.Default<String>(ParameterParser.VALUE, "100.0"),
                                new Attribute.Default<String>(ParameterParser.LOWER, "0.0")
                        }, true);
                writer.writeCloseTag(LogNormalDistributionModelParser.STDEV);
                writer.writeCloseTag(NormalDistributionModelParser.NORMAL_DISTRIBUTION_MODEL);
                writer.writeCloseTag(DistributionLikelihoodParser.DISTRIBUTION);
                writer.writeCloseTag(tag);

                // fixed effects priors
                counter = 1;
                for (Taxa taxonSet : options.taxonSets) {
                    writer.writeOpenTag(tag,
                            new Attribute[]{
                                    new Attribute.Default<String>(XMLParser.ID, BranchSpecificFixedEffectsParser.FIXED_EFFECTS_LIKELIHOOD + counter)
                            }
                    );
                    writer.writeOpenTag(DistributionLikelihoodParser.DATA);
                    writer.writeIDref(PARAMETER, BranchSpecificFixedEffectsParser.COEFFICIENT + counter);
                    writer.writeCloseTag(DistributionLikelihoodParser.DATA);
                    writer.writeOpenTag(DistributionLikelihoodParser.DISTRIBUTION);
                    writer.writeOpenTag(NormalDistributionModelParser.NORMAL_DISTRIBUTION_MODEL);
                    writer.writeOpenTag(LogNormalDistributionModelParser.MEAN);
                    writer.writeTag(PARAMETER,
                            new Attribute[]{
                                    new Attribute.Default<String>(ParameterParser.VALUE, "0.0"),
                                    new Attribute.Default<String>(ParameterParser.LOWER, "0.0")
                            }, true);
                    writer.writeCloseTag(LogNormalDistributionModelParser.MEAN);
                    writer.writeOpenTag(LogNormalDistributionModelParser.STDEV);
                    writer.writeTag(PARAMETER,
                            new Attribute[]{
                                    new Attribute.Default<String>(XMLParser.ID, "betweenGroupStd" + counter),
                                    new Attribute.Default<String>(ParameterParser.VALUE, "100.0"),
                                    new Attribute.Default<String>(ParameterParser.LOWER, "0.0")
                            }, true);
                    writer.writeCloseTag(LogNormalDistributionModelParser.STDEV);
                    writer.writeCloseTag(NormalDistributionModelParser.NORMAL_DISTRIBUTION_MODEL);
                    writer.writeCloseTag(DistributionLikelihoodParser.DISTRIBUTION);
                    writer.writeCloseTag(tag);
                    counter++;
                }

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
        writer.writeIDref(DefaultTreeModel.TREE_MODEL, treePrefix + DefaultTreeModel.TREE_MODEL);
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
        writer.writeIDref(DefaultTreeModel.TREE_MODEL, treePrefix + DefaultTreeModel.TREE_MODEL);
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
        writer.writeIDref(DefaultTreeModel.TREE_MODEL, treePrefix + DefaultTreeModel.TREE_MODEL);
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

            case HMC_CLOCK:
                tag = ArbitraryBranchRates.BRANCH_RATES;
                id = model.getPrefix() + BranchRateModel.BRANCH_RATES;
                break;

            case SHRINKAGE_LOCAL_CLOCK:
                writer.writeIDref(GAMMA_PRIOR, "globalScalePrior");
                writer.writeIDref(AutoCorrelatedBranchRatesDistributionParser.AUTO_CORRELATED_RATES, "substBranchRatesPrior");
                tag = ScaledByTreeTimeBranchRateModelParser.TREE_TIME_BRANCH_RATES;
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

            case MIXED_EFFECTS_CLOCK:
                //always write distribution likelihoods for rate, scale and intercept
                //writer.writeIDref(DistributionLikelihood.DISTRIBUTION_LIKELIHOOD, BranchSpecificFixedEffects.RATES_PRIOR);
                //writer.writeIDref(DistributionLikelihood.DISTRIBUTION_LIKELIHOOD, BranchSpecificFixedEffects.SCALE_PRIOR);
                //writer.writeIDref(DistributionLikelihood.DISTRIBUTION_LIKELIHOOD, BranchSpecificFixedEffects.INTERCEPT_PRIOR);
                //check for coefficients
                /*String coeff = BranchSpecificFixedEffectsParser.COEFFICIENT;
                int number = 1;
                String concat = coeff + number;
                while (model.hasParameter(concat)) {
                    writer.writeIDref(DistributionLikelihood.DISTRIBUTION_LIKELIHOOD, BranchSpecificFixedEffectsParser.FIXED_EFFECTS_LIKELIHOOD + number);
                    number++;
                    concat = coeff + number;
                }*/
                tag = ArbitraryBranchRatesParser.ARBITRARY_BRANCH_RATES;
                id = model.getPrefix() + ArbitraryBranchRates.BRANCH_RATES;
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
        String parameterName = options.useNuRelativeRates() ? "allNus" : "allMus";

        Parameter allMus = model.getParameter(parameterName);
        if (allMus.getSubParameters().size() > 1) {
            writer.writeComment("Collecting together relative rates for partitions");

            // allMus is global for each gene
            writer.writeOpenTag(CompoundParameterParser.COMPOUND_PARAMETER,
                    new Attribute[]{new Attribute.Default<String>(XMLParser.ID, model.getPrefix() + parameterName)});

            for (Parameter parameter : allMus.getSubParameters()) {
                writer.writeIDref(PARAMETER, parameter.getName());
            }

            writer.writeCloseTag(CompoundParameterParser.COMPOUND_PARAMETER);
            writer.writeText("");
        }
    }

    public void writeAllClockRateRefs(PartitionClockModel model, XMLWriter writer) {
        writer.writeIDref(PARAMETER, getClockRateString(model));
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
            case HMC_CLOCK:
                switch (model.getClockDistributionType()) {
                    case LOGNORMAL:
                        return prefix + ClockType.HMC_CLOCK_LOCATION;
                }
            case SHRINKAGE_LOCAL_CLOCK:
                return prefix + ClockType.SHRINKAGE_CLOCK_LOCATION;
            case MIXED_EFFECTS_CLOCK:
                return prefix + ClockType.ME_CLOCK_LOCATION;
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

        if (options.useNuRelativeRates()) {
            Parameter allNus = model.getParameter("allNus");
            if (allNus.getSubParameters().size() > 1) {
                // The mu's are the more relevant parameter and allow comparisons with the old parameterization
                // May be confusing to log the nus and mus, but necessary for use with generalized stepping-stone sampling
                writer.writeIDref(CompoundParameterParser.COMPOUND_PARAMETER, prefix + "allNus");

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
                writer.writeIDref(PARAMETER, prefix + "clock.rate");
                break;

            case FIXED_LOCAL_CLOCK:
                writer.writeIDref(PARAMETER, prefix + "clock.rate");
                for (Taxa taxonSet : options.taxonSets) {
                    if (options.taxonSetsMono.get(taxonSet)) {
                        String parameterName = taxonSet.getId() + ".rate";
                        writer.writeIDref(PARAMETER, model.getPrefix() + parameterName);
                    }
                }
                break;

            case UNCORRELATED:

                if (model.performModelAveraging()) {

                    writer.writeIDref(PARAMETER, prefix + ClockType.UCLD_MEAN);
                    writer.writeIDref(PARAMETER, prefix + ClockType.UCLD_STDEV);
                    writer.writeIDref(PARAMETER, prefix + ClockType.UCGD_MEAN);
                    writer.writeIDref(PARAMETER, prefix + ClockType.UCGD_SHAPE);
                    writer.writeIDref(PARAMETER, prefix + ClockType.UCED_MEAN);

                    writer.writeIDref(PARAMETER, "branchRates.distributionIndex");
                    writer.writeIDref(PARAMETER, "branchRates.quantiles");

                } else {

                    switch (model.getClockDistributionType()) {
                        case LOGNORMAL:
                            writer.writeIDref(PARAMETER, prefix + ClockType.UCLD_MEAN);
                            writer.writeIDref(PARAMETER, prefix + ClockType.UCLD_STDEV);
                            break;
                        case GAMMA:
                            writer.writeIDref(PARAMETER, prefix + ClockType.UCGD_MEAN);
                            writer.writeIDref(PARAMETER, prefix + ClockType.UCGD_SHAPE);
                            break;
                        case CAUCHY:
                            throw new UnsupportedOperationException("Uncorrelated Couchy model not supported yet");
//                        break;
                        case EXPONENTIAL:
                            writer.writeIDref(PARAMETER, prefix + ClockType.UCED_MEAN);
                            break;
                    }

                }
                break;

            case SHRINKAGE_LOCAL_CLOCK:
                writer.writeIDref(PARAMETER, prefix + ClockType.SHRINKAGE_CLOCK_LOCATION);
                break;

            case MIXED_EFFECTS_CLOCK:
                writer.writeIDref(PARAMETER, prefix + ClockType.ME_CLOCK_LOCATION);
                break;

            case HMC_CLOCK:
                switch (model.getClockDistributionType()) {
                    case LOGNORMAL:
                        writer.writeIDref(PARAMETER, prefix + ClockType.HMC_CLOCK_LOCATION);
                        writer.writeIDref(PARAMETER, prefix + ClockType.HMCLN_SCALE);
                        break;
                }
                break;

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
            case HMC_CLOCK:
            case FIXED_LOCAL_CLOCK:
                writer.writeIDref(RateStatisticParser.RATE_STATISTIC, prefix + "meanRate");
                writer.writeIDref(RateStatisticParser.RATE_STATISTIC, prefix + RateStatisticParser.COEFFICIENT_OF_VARIATION);
                writer.writeIDref(RateCovarianceStatisticParser.RATE_COVARIANCE_STATISTIC, prefix + "covariance");
                break;

            case MIXED_EFFECTS_CLOCK:
            case SHRINKAGE_LOCAL_CLOCK:
                break;

            case AUTOCORRELATED:
                writer.writeIDref(RateStatisticParser.RATE_STATISTIC, prefix + "meanRate");
                writer.writeIDref(RateStatisticParser.RATE_STATISTIC, prefix + RateStatisticParser.COEFFICIENT_OF_VARIATION);
                writer.writeIDref(RateCovarianceStatisticParser.RATE_COVARIANCE_STATISTIC, prefix + "covariance");
                writer.writeIDref(PARAMETER, prefix + "branchRates.var");
                writer.writeIDref(PARAMETER, prefix + "treeModel.rootRate");
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
