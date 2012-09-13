/*
 * Copyright (C) 2002-2009 Alexei Drummond and Andrew Rambaut
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
 * BEAST is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package dr.app.beauti.generator;

import dr.app.beauti.components.ComponentFactory;
import dr.app.beauti.options.*;
import dr.app.beauti.types.ClockType;
import dr.app.beauti.util.XMLWriter;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.clock.RateEvolutionLikelihood;
import dr.evomodel.tree.TreeModel;
import dr.evomodelxml.branchratemodel.DiscretizedBranchRatesParser;
import dr.evomodelxml.branchratemodel.RandomLocalClockModelParser;
import dr.evomodelxml.branchratemodel.StrictClockBranchRatesParser;
import dr.evomodelxml.clock.ACLikelihoodParser;
import dr.evomodelxml.tree.RateCovarianceStatisticParser;
import dr.evomodelxml.tree.RateStatisticParser;
import dr.evomodelxml.tree.TreeModelParser;
import dr.inference.distribution.ExponentialDistributionModel;
import dr.inference.model.ParameterParser;
import dr.inferencexml.distribution.LogNormalDistributionModelParser;
import dr.inferencexml.model.CompoundParameterParser;
import dr.inferencexml.model.SumStatisticParser;
import dr.util.Attribute;
import dr.xml.XMLParser;

import java.util.List;

/**
 * @author Alexei Drummond
 * @author Andrew Rambaut
 */
public class BranchRatesModelGenerator extends Generator {

    public BranchRatesModelGenerator(BeautiOptions options, ComponentFactory[] components) {
        super(options, components);
    }

    /**
     * Write the relaxed clock branch rates block.
     *
     * @param model  PartitionClockModel
     * @param writer the writer
     */
    public void writeBranchRatesModel(PartitionClockModel model, XMLWriter writer) {

        setModelPrefix(model.getPrefix());

        Attribute[] attributes;
        int categoryCount = 0;
        String treePrefix;
        List<PartitionTreeModel> activeTrees = options.getPartitionTreeModels(options.getDataPartitions(model));

        switch (model.getClockType()) {
            case STRICT_CLOCK:
                writer.writeComment("The strict clock (Uniform rates across branches)");

                writer.writeOpenTag(
                        StrictClockBranchRatesParser.STRICT_CLOCK_BRANCH_RATES,
                        new Attribute[]{new Attribute.Default<String>(XMLParser.ID, modelPrefix + BranchRateModel.BRANCH_RATES)}
                );
                writeParameter("rate", "clock.rate", model, writer);
                writer.writeCloseTag(StrictClockBranchRatesParser.STRICT_CLOCK_BRANCH_RATES);

                break;

            case UNCORRELATED:
                writer.writeComment("The uncorrelated relaxed clock (Drummond, Ho, Phillips & Rambaut (2006) PLoS Biology 4, e88 )");

                for (PartitionTreeModel tree : activeTrees) {
                    treePrefix = tree.getPrefix();

                    PartitionClockModelTreeModelLink clockTree = options.getPartitionClockTreeLink(model, tree);
                    if (clockTree == null) {
                        throw new IllegalArgumentException("Cannot find PartitionClockTreeLink, given clock model = " + model.getName()
                                + ", tree model = " + tree.getName());
                    }

                    //if (options.isFixedSubstitutionRate()) {
                    //    attributes = new Attribute[]{
                    //            new Attribute.Default<String>(XMLParser.ID, BranchRateModel.BRANCH_RATES),
                    //            new Attribute.Default<Double>(DiscretizedBranchRatesParser.NORMALIZED_MEAN, options.getMeanSubstitutionRate())
                    //    };
                    //} else {
                    attributes = new Attribute[]{new Attribute.Default<String>(XMLParser.ID, options.noDuplicatedPrefix(modelPrefix, treePrefix)
                            + BranchRateModel.BRANCH_RATES)};
                    //}
                    writer.writeOpenTag(DiscretizedBranchRatesParser.DISCRETIZED_BRANCH_RATES, attributes);
                    // tree
                    writer.writeIDref(TreeModel.TREE_MODEL, treePrefix + TreeModel.TREE_MODEL);

                    writer.writeOpenTag("distribution");

                    switch (model.getClockDistributionType()) {

                        case LOGNORMAL:
                            writer.writeOpenTag(LogNormalDistributionModelParser.LOGNORMAL_DISTRIBUTION_MODEL,
                                    new Attribute.Default<String>(LogNormalDistributionModelParser.MEAN_IN_REAL_SPACE, "true"));

                            if (activeTrees.indexOf(tree) < 1) {
                                writeParameter("mean", ClockType.UCLD_MEAN, model, writer);
                                writeParameter("stdev", ClockType.UCLD_STDEV, model, writer);
                            } else {
                                writeParameterRef("mean", modelPrefix + ClockType.UCLD_MEAN, writer);
                                writeParameterRef("stdev", modelPrefix + ClockType.UCLD_STDEV, writer);
                            }

                            writer.writeCloseTag(LogNormalDistributionModelParser.LOGNORMAL_DISTRIBUTION_MODEL);
                            break;
                        case GAMMA:
                            throw new UnsupportedOperationException("Uncorrelated gamma model not implemented yet");
//                            break;
                        case CAUCHY:
                            throw new UnsupportedOperationException("Uncorrelated Cauchy model not implemented yet");
//                            break;
                        case EXPONENTIAL:
                            writer.writeOpenTag(ExponentialDistributionModel.EXPONENTIAL_DISTRIBUTION_MODEL);

                            if (activeTrees.indexOf(tree) < 1) {
                                writeParameter("mean", ClockType.UCED_MEAN, model, writer);
                            } else {
                                writeParameterRef("mean", modelPrefix + ClockType.UCED_MEAN, writer);
                            }

                            writer.writeCloseTag(ExponentialDistributionModel.EXPONENTIAL_DISTRIBUTION_MODEL);
                            break;
                    }

                    writer.writeCloseTag("distribution");

                    writer.writeOpenTag(DiscretizedBranchRatesParser.RATE_CATEGORIES);
                    // AR - this parameter will now set its dimension automatically when BEAST is run
//                    if (!options.hasIdenticalTaxa()) {
//                        for (AbstractPartitionData dataPartition : options.dataPartitions) {
//                            if (dataPartition.getPartitionClockModel().equals(model)) {
//                                categoryCount = (dataPartition.getTaxonCount() - 1) * 2;
//                            }
//                        }
//                    } else {
//                        categoryCount = (options.taxonList.getTaxonCount() - 1) * 2;
//                    }
//                    writeParameter(clockTree.getParameter("branchRates.categories"), categoryCount, writer);
                    writeParameter(clockTree.getParameter("branchRates.categories"), -1, writer);
                    writer.writeCloseTag(DiscretizedBranchRatesParser.RATE_CATEGORIES);
                    writer.writeCloseTag(DiscretizedBranchRatesParser.DISCRETIZED_BRANCH_RATES);

                    writer.writeText("");
                    writer.writeOpenTag(
                            RateStatisticParser.RATE_STATISTIC,
                            new Attribute[]{
                                    new Attribute.Default<String>(XMLParser.ID, options.noDuplicatedPrefix(modelPrefix, treePrefix) + "meanRate"),
                                    new Attribute.Default<String>("name", options.noDuplicatedPrefix(modelPrefix, treePrefix) + "meanRate"),
                                    new Attribute.Default<String>("mode", "mean"),
                                    new Attribute.Default<String>("internal", "true"),
                                    new Attribute.Default<String>("external", "true")
                            }
                    );
                    writer.writeIDref(TreeModel.TREE_MODEL, treePrefix + TreeModel.TREE_MODEL);
                    writer.writeIDref(DiscretizedBranchRatesParser.DISCRETIZED_BRANCH_RATES, options.noDuplicatedPrefix(modelPrefix, treePrefix)
                            + BranchRateModel.BRANCH_RATES);
                    writer.writeCloseTag(RateStatisticParser.RATE_STATISTIC);

                    writer.writeText("");
                    writer.writeOpenTag(
                            RateStatisticParser.RATE_STATISTIC,
                            new Attribute[]{
                                    new Attribute.Default<String>(XMLParser.ID, options.noDuplicatedPrefix(modelPrefix, treePrefix) + RateStatisticParser.COEFFICIENT_OF_VARIATION),
                                    new Attribute.Default<String>("name", options.noDuplicatedPrefix(modelPrefix, treePrefix) + RateStatisticParser.COEFFICIENT_OF_VARIATION),
                                    new Attribute.Default<String>("mode", RateStatisticParser.COEFFICIENT_OF_VARIATION),
                                    new Attribute.Default<String>("internal", "true"),
                                    new Attribute.Default<String>("external", "true")
                            }
                    );
                    writer.writeIDref(TreeModel.TREE_MODEL, treePrefix + TreeModel.TREE_MODEL);
                    writer.writeIDref(DiscretizedBranchRatesParser.DISCRETIZED_BRANCH_RATES, options.noDuplicatedPrefix(modelPrefix, treePrefix)
                            + BranchRateModel.BRANCH_RATES);
                    writer.writeCloseTag(RateStatisticParser.RATE_STATISTIC);

                    writer.writeText("");
                    writer.writeOpenTag(
                            RateCovarianceStatisticParser.RATE_COVARIANCE_STATISTIC,
                            new Attribute[]{
                                    new Attribute.Default<String>(XMLParser.ID, options.noDuplicatedPrefix(modelPrefix, treePrefix) + "covariance"),
                                    new Attribute.Default<String>("name", options.noDuplicatedPrefix(modelPrefix, treePrefix) + "covariance")
                            }
                    );
                    writer.writeIDref(TreeModel.TREE_MODEL, treePrefix + TreeModel.TREE_MODEL);
                    writer.writeIDref(DiscretizedBranchRatesParser.DISCRETIZED_BRANCH_RATES, options.noDuplicatedPrefix(modelPrefix, treePrefix) + BranchRateModel.BRANCH_RATES);
                    writer.writeCloseTag(RateCovarianceStatisticParser.RATE_COVARIANCE_STATISTIC);
                }

                break;

            case AUTOCORRELATED:
                writer.writeComment("The autocorrelated relaxed clock (Rannala & Yang, 2007)");

                for (PartitionTreeModel tree : activeTrees) {
                    treePrefix = tree.getPrefix();

                    PartitionClockModelTreeModelLink clockTree = options.getPartitionClockTreeLink(model, tree);
                    if (clockTree == null) {
                        throw new IllegalArgumentException("Cannot find PartitionClockTreeLink, given clock model = " + model.getName()
                                + ", tree model = " + tree.getName());
                    }

                    attributes = new Attribute[]{
                            new Attribute.Default<String>(XMLParser.ID, options.noDuplicatedPrefix(modelPrefix, treePrefix) + BranchRateModel.BRANCH_RATES),
                            new Attribute.Default<String>("episodic", "false"),
                            new Attribute.Default<String>("logspace", "true"),
                    };

                    writer.writeOpenTag(ACLikelihoodParser.AC_LIKELIHOOD, attributes);
                    writer.writeIDref(TreeModel.TREE_MODEL, treePrefix + TreeModel.TREE_MODEL);

                    if (!model.isEstimatedRate()) { //TODO move to options or panel select method
                        Parameter para = tree.getParameter(TreeModel.TREE_MODEL + "." + RateEvolutionLikelihood.ROOTRATE);//"treeModel.rootRate"
                        para.isFixed = true;
                        para.initial = model.getRate();
                    }

                    writer.writeOpenTag(RateEvolutionLikelihood.RATES,
                            new Attribute[]{
                                    new Attribute.Default<String>(TreeModelParser.ROOT_NODE, "false"),
                                    new Attribute.Default<String>(TreeModelParser.INTERNAL_NODES, "true"),
                                    new Attribute.Default<String>(TreeModelParser.LEAF_NODES, "true")
                            });
                    writer.writeTag(ParameterParser.PARAMETER,
                            new Attribute.Default<String>(XMLParser.ID, options.noDuplicatedPrefix(modelPrefix, treePrefix) + TreeModel.TREE_MODEL + "."
                                    + TreeModelParser.NODE_RATES), true);
                    writer.writeCloseTag(RateEvolutionLikelihood.RATES);

                    writer.writeOpenTag(RateEvolutionLikelihood.ROOTRATE,
                            new Attribute[]{
                                    new Attribute.Default<String>(TreeModelParser.ROOT_NODE, "true"),
                                    new Attribute.Default<String>(TreeModelParser.INTERNAL_NODES, "false"),
                                    new Attribute.Default<String>(TreeModelParser.LEAF_NODES, "false")
                            });
                    writer.writeTag(ParameterParser.PARAMETER,
                            new Attribute.Default<String>(XMLParser.ID, options.noDuplicatedPrefix(modelPrefix, treePrefix) + TreeModel.TREE_MODEL + "."
                                    + RateEvolutionLikelihood.ROOTRATE), true);
                    writer.writeCloseTag(RateEvolutionLikelihood.ROOTRATE);
                    //                writeParameterRef("rates", treePrefix + "treeModel.nodeRates", writer);
                    //                writeParameterRef(RateEvolutionLikelihood.ROOTRATE, treePrefix + "treeModel.rootRate", writer);
                    writeParameter("variance", "branchRates.var", clockTree, writer);

                    writer.writeCloseTag(ACLikelihoodParser.AC_LIKELIHOOD);

                    if (model.isEstimatedRate()) {//TODO
                        writer.writeText("");
                        writer.writeOpenTag(CompoundParameterParser.COMPOUND_PARAMETER,
                                new Attribute[]{new Attribute.Default<String>(XMLParser.ID, options.noDuplicatedPrefix(modelPrefix, treePrefix) + TreeModel.TREE_MODEL
                                        + "." + "allRates")});
                        writer.writeIDref(ParameterParser.PARAMETER, options.noDuplicatedPrefix(modelPrefix, treePrefix) + TreeModel.TREE_MODEL + "."
                                + TreeModelParser.NODE_RATES);
                        writer.writeIDref(ParameterParser.PARAMETER, options.noDuplicatedPrefix(modelPrefix, treePrefix) + TreeModel.TREE_MODEL + "."
                                + RateEvolutionLikelihood.ROOTRATE);
                        writer.writeCloseTag(CompoundParameterParser.COMPOUND_PARAMETER);
                    }

                    writer.writeText("");
                    writer.writeOpenTag(
                            RateStatisticParser.RATE_STATISTIC,
                            new Attribute[]{
                                    new Attribute.Default<String>(XMLParser.ID, options.noDuplicatedPrefix(modelPrefix, treePrefix) + "meanRate"),
                                    new Attribute.Default<String>("name", options.noDuplicatedPrefix(modelPrefix, treePrefix) + "meanRate"),
                                    new Attribute.Default<String>("mode", "mean"),
                                    new Attribute.Default<String>("internal", "true"),
                                    new Attribute.Default<String>("external", "true")
                            }
                    );
                    writer.writeIDref(TreeModel.TREE_MODEL, treePrefix + TreeModel.TREE_MODEL);
                    writer.writeIDref(ACLikelihoodParser.AC_LIKELIHOOD, modelPrefix + BranchRateModel.BRANCH_RATES);
                    writer.writeCloseTag(RateStatisticParser.RATE_STATISTIC);

                    writer.writeText("");
                    writer.writeOpenTag(
                            RateStatisticParser.RATE_STATISTIC,
                            new Attribute[]{
                                    new Attribute.Default<String>(XMLParser.ID, modelPrefix + RateStatisticParser.COEFFICIENT_OF_VARIATION),
                                    new Attribute.Default<String>("name", options.noDuplicatedPrefix(modelPrefix, treePrefix) + RateStatisticParser.COEFFICIENT_OF_VARIATION),
                                    new Attribute.Default<String>("mode", RateStatisticParser.COEFFICIENT_OF_VARIATION),
                                    new Attribute.Default<String>("internal", "true"),
                                    new Attribute.Default<String>("external", "true")
                            }
                    );
                    writer.writeIDref(TreeModel.TREE_MODEL, treePrefix + TreeModel.TREE_MODEL);
                    writer.writeIDref(ACLikelihoodParser.AC_LIKELIHOOD, modelPrefix + BranchRateModel.BRANCH_RATES);
                    writer.writeCloseTag(RateStatisticParser.RATE_STATISTIC);

                    writer.writeText("");
                    writer.writeOpenTag(
                            RateCovarianceStatisticParser.RATE_COVARIANCE_STATISTIC,
                            new Attribute[]{
                                    new Attribute.Default<String>(XMLParser.ID, options.noDuplicatedPrefix(modelPrefix, treePrefix) + "covariance"),
                                    new Attribute.Default<String>("name", options.noDuplicatedPrefix(modelPrefix, treePrefix) + "covariance")
                            }
                    );
                    writer.writeIDref(TreeModel.TREE_MODEL, treePrefix + TreeModel.TREE_MODEL);
                    writer.writeIDref(ACLikelihoodParser.AC_LIKELIHOOD, modelPrefix + BranchRateModel.BRANCH_RATES);
                    writer.writeCloseTag(RateCovarianceStatisticParser.RATE_COVARIANCE_STATISTIC);
                }

                break;

            case RANDOM_LOCAL_CLOCK: // 1 random local clock CANNOT have different tree models
                writer.writeComment("The random local clock model (Drummond & Suchard, 2010)");

                if (activeTrees == null || activeTrees.size() != 1) {
                    throw new IllegalArgumentException("A single random local clock cannot be applied to multiple trees.");
                }
                treePrefix = activeTrees.get(0).getPrefix();

                writer.writeOpenTag(
                        RandomLocalClockModelParser.LOCAL_BRANCH_RATES,
                        new Attribute[]{
                                // 1 random local clock CANNOT have different tree models,
                                // so use modelPrefix not noDuplicatedPrefix(modelPrefix, treePrefix)
                                new Attribute.Default<String>(XMLParser.ID, modelPrefix + BranchRateModel.BRANCH_RATES),
                                new Attribute.Default<String>("ratesAreMultipliers", "false")
                        }
                );
                writer.writeIDref(TreeModel.TREE_MODEL, treePrefix + TreeModel.TREE_MODEL);

                writer.writeOpenTag("rates");
                writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default<String>
                        (XMLParser.ID, modelPrefix + ClockType.LOCAL_CLOCK + ".relativeRates")
                        , true);
                writer.writeCloseTag("rates");
                writer.writeOpenTag("rateIndicator");
                writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default<String>
                        (XMLParser.ID, modelPrefix + ClockType.LOCAL_CLOCK + ".changes")
                        , true);
                writer.writeCloseTag("rateIndicator");

                writeParameter("clockRate", "clock.rate", model, writer);

                writer.writeCloseTag(RandomLocalClockModelParser.LOCAL_BRANCH_RATES);

                writer.writeText("");
                writer.writeOpenTag(
                        SumStatisticParser.SUM_STATISTIC,
                        new Attribute[]{
                                new Attribute.Default<String>(XMLParser.ID, modelPrefix + "rateChanges"),
                                new Attribute.Default<String>("name", modelPrefix + "rateChangeCount"),
                                new Attribute.Default<String>("elementwise", "true"),
                        }
                );
                writer.writeIDref(ParameterParser.PARAMETER, modelPrefix + ClockType.LOCAL_CLOCK + ".changes");
                writer.writeCloseTag(SumStatisticParser.SUM_STATISTIC);

                writer.writeText("");

                writer.writeOpenTag(
                        RateStatisticParser.RATE_STATISTIC,
                        new Attribute[]{
                                new Attribute.Default<String>(XMLParser.ID, options.noDuplicatedPrefix(modelPrefix, treePrefix) + "meanRate"),
                                new Attribute.Default<String>("name", options.noDuplicatedPrefix(modelPrefix, treePrefix) + "meanRate"),
                                new Attribute.Default<String>("mode", "mean"),
                                new Attribute.Default<String>("internal", "true"),
                                new Attribute.Default<String>("external", "true")
                        }
                );
                writer.writeIDref(TreeModel.TREE_MODEL, treePrefix + TreeModel.TREE_MODEL);
                writer.writeIDref(RandomLocalClockModelParser.LOCAL_BRANCH_RATES, modelPrefix + BranchRateModel.BRANCH_RATES);
                writer.writeCloseTag(RateStatisticParser.RATE_STATISTIC);

                writer.writeText("");
                writer.writeOpenTag(
                        RateStatisticParser.RATE_STATISTIC,
                        new Attribute[]{
                                new Attribute.Default<String>(XMLParser.ID, options.noDuplicatedPrefix(modelPrefix, treePrefix) + RateStatisticParser.COEFFICIENT_OF_VARIATION),
                                new Attribute.Default<String>("name", options.noDuplicatedPrefix(modelPrefix, treePrefix) + RateStatisticParser.COEFFICIENT_OF_VARIATION),
                                new Attribute.Default<String>("mode", RateStatisticParser.COEFFICIENT_OF_VARIATION),
                                new Attribute.Default<String>("internal", "true"),
                                new Attribute.Default<String>("external", "true")
                        }
                );
                writer.writeIDref(TreeModel.TREE_MODEL, treePrefix + TreeModel.TREE_MODEL);
                writer.writeIDref(RandomLocalClockModelParser.LOCAL_BRANCH_RATES, modelPrefix + BranchRateModel.BRANCH_RATES);
                writer.writeCloseTag(RateStatisticParser.RATE_STATISTIC);

                writer.writeText("");
                writer.writeOpenTag(
                        RateCovarianceStatisticParser.RATE_COVARIANCE_STATISTIC,
                        new Attribute[]{
                                new Attribute.Default<String>(XMLParser.ID, options.noDuplicatedPrefix(modelPrefix, treePrefix) + "covariance"),
                                new Attribute.Default<String>("name", options.noDuplicatedPrefix(modelPrefix, treePrefix) + "covariance")
                        }
                );
                writer.writeIDref(TreeModel.TREE_MODEL, treePrefix + TreeModel.TREE_MODEL);
                writer.writeIDref(RandomLocalClockModelParser.LOCAL_BRANCH_RATES, modelPrefix + BranchRateModel.BRANCH_RATES);
                writer.writeCloseTag(RateCovarianceStatisticParser.RATE_COVARIANCE_STATISTIC);

                break;

            default:
                throw new IllegalArgumentException("Unknown clock model");
        }

    }

    public void writeAllClockRateRefs(PartitionClockModel model, XMLWriter writer) {
        writer.writeIDref(ParameterParser.PARAMETER, getClockRateString(model));
    }

    public String getClockRateString(PartitionClockModel model) {
        setModelPrefix(model.getPrefix());

        switch (model.getClockType()) {
            case STRICT_CLOCK:
            case RANDOM_LOCAL_CLOCK:
                return modelPrefix + "clock.rate";

            case UNCORRELATED:
                switch (model.getClockDistributionType()) {

                    case LOGNORMAL:
                        return modelPrefix + ClockType.UCLD_MEAN;
                    case GAMMA:
                        throw new UnsupportedOperationException("Uncorrelated gamma model not supported yet");
//                        return modelPrefix + ClockType.UCGD_SCALE;
                    case CAUCHY:
                        throw new UnsupportedOperationException("Uncorrelated Cauchy model not supported yet");
//                        return modelPrefix + ClockType.UCCD_MEAN;
                    case EXPONENTIAL:
                        return modelPrefix + ClockType.UCED_MEAN;
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
        setModelPrefix(model.getPrefix());

        switch (model.getClockType()) {
            case STRICT_CLOCK:
            case RANDOM_LOCAL_CLOCK:
                writer.writeIDref(ParameterParser.PARAMETER, model.getPrefix() + "clock.rate");
                break;

            case UNCORRELATED:
                switch (model.getClockDistributionType()) {
                    case LOGNORMAL:
                        writer.writeIDref(ParameterParser.PARAMETER, model.getPrefix() + ClockType.UCLD_MEAN);
                        writer.writeIDref(ParameterParser.PARAMETER, model.getPrefix() + ClockType.UCLD_STDEV);
                        break;
                    case GAMMA:
                        throw new UnsupportedOperationException("Uncorrelated gamma model not supported yet");
                    case CAUCHY:
                        throw new UnsupportedOperationException("Uncorrelated Couchy model not supported yet");
                    case EXPONENTIAL:
                        writer.writeIDref(ParameterParser.PARAMETER, model.getPrefix() + ClockType.UCED_MEAN);
                        break;
                }

            case AUTOCORRELATED:
// TODO
                break;

            default:
                throw new IllegalArgumentException("Unknown clock model");
        }

    }

    public void writeLogStatistic(PartitionClockModel model, XMLWriter writer) {
        setModelPrefix(model.getPrefix());

        switch (model.getClockType()) {
            case STRICT_CLOCK:
                break;

            case UNCORRELATED:
                for (PartitionTreeModel tree : options.getPartitionTreeModels(options.getDataPartitions(model))) {
                    writer.writeIDref(RateStatisticParser.RATE_STATISTIC, options.noDuplicatedPrefix(model.getPrefix(), tree.getPrefix()) + "meanRate");
                    writer.writeIDref(RateStatisticParser.RATE_STATISTIC, options.noDuplicatedPrefix(model.getPrefix(), tree.getPrefix()) + RateStatisticParser.COEFFICIENT_OF_VARIATION);
                    writer.writeIDref(RateCovarianceStatisticParser.RATE_COVARIANCE_STATISTIC, options.noDuplicatedPrefix(model.getPrefix(), tree.getPrefix()) + "covariance");
                }
                break;

            case AUTOCORRELATED:
// TODO
                for (PartitionTreeModel tree : options.getPartitionTreeModels(options.getDataPartitions(model))) {
                    writer.writeIDref(RateStatisticParser.RATE_STATISTIC, options.noDuplicatedPrefix(model.getPrefix(), tree.getPrefix()) + "meanRate");
                    writer.writeIDref(RateStatisticParser.RATE_STATISTIC, options.noDuplicatedPrefix(model.getPrefix(), tree.getPrefix()) + RateStatisticParser.COEFFICIENT_OF_VARIATION);
                    writer.writeIDref(RateCovarianceStatisticParser.RATE_COVARIANCE_STATISTIC, options.noDuplicatedPrefix(model.getPrefix(), tree.getPrefix()) + "covariance");
                    writer.writeIDref(ParameterParser.PARAMETER, options.noDuplicatedPrefix(model.getPrefix(), tree.getPrefix()) + "branchRates.var");
                    writer.writeIDref(ParameterParser.PARAMETER, options.noDuplicatedPrefix(model.getPrefix(), tree.getPrefix()) + "treeModel.rootRate");
                }
                break;

            case RANDOM_LOCAL_CLOCK:
                writer.writeIDref(SumStatisticParser.SUM_STATISTIC, model.getPrefix() + "rateChanges");
                for (PartitionTreeModel tree : options.getPartitionTreeModels(options.getDataPartitions(model))) {
                    writer.writeIDref(RateStatisticParser.RATE_STATISTIC, options.noDuplicatedPrefix(model.getPrefix(), tree.getPrefix()) + "meanRate");
                    writer.writeIDref(RateStatisticParser.RATE_STATISTIC, options.noDuplicatedPrefix(model.getPrefix(), tree.getPrefix()) + RateStatisticParser.COEFFICIENT_OF_VARIATION);
                    writer.writeIDref(RateCovarianceStatisticParser.RATE_COVARIANCE_STATISTIC, options.noDuplicatedPrefix(model.getPrefix(), tree.getPrefix()) + "covariance");
                }
                break;

            default:
                throw new IllegalArgumentException("Unknown clock model");
        }

    }

    public void writeClockLikelihoodReferences(XMLWriter writer) {
        for (AbstractPartitionData partition : options.dataPartitions) { // Each PD has one TreeLikelihood
            PartitionClockModel clockModel = partition.getPartitionClockModel();

            if (clockModel != null && clockModel.getClockType() == ClockType.AUTOCORRELATED) {
                throw new UnsupportedOperationException("Autocorrelated relaxed clock model not implemented yet");
//                writer.writeIDref(ACLikelihoodParser.AC_LIKELIHOOD, clockModel.getPrefix() + BranchRateModel.BRANCH_RATES);
            }
        }
    }

}
