package dr.app.beauti.generator;

import dr.app.beauti.XMLWriter;
import dr.app.beauti.options.BeautiOptions;
import dr.app.beauti.options.ClockType;
import dr.evomodel.branchratemodel.RandomLocalClockModel;
import dr.evomodel.branchratemodel.StrictClockBranchRates;
import dr.evomodel.tree.RateCovarianceStatistic;
import dr.evomodel.tree.RateStatistic;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.clock.ACLikelihood;
import dr.evomodelxml.DiscretizedBranchRatesParser;
import dr.inference.distribution.ExponentialDistributionModel;
import dr.inference.distribution.LogNormalDistributionModel;
import dr.inference.model.SumStatistic;
import dr.util.Attribute;

/**
 * @author Alexei Drummond
 * @author Andrew Rambaut
 */
public class BranchRatesModelGenerator extends Generator {

    public BranchRatesModelGenerator(BeautiOptions options) {
        super(options);
    }

    /**
     * Write the relaxed clock branch rates block.
     *
     * @param writer the writer
     */
    public void writeBranchRatesModel(XMLWriter writer) {

        Attribute[] attributes;
        int categoryCount;

        switch (options.clockType) {
            case STRICT_CLOCK:
                if (options.isFixedSubstitutionRate()) {

                    fixParameter("clock.rate", options.getMeanSubstitutionRate());
                }

                writer.writeComment("The strict clock (Uniform rates across branches)");
                writer.writeOpenTag(
                        StrictClockBranchRates.STRICT_CLOCK_BRANCH_RATES,
                        new Attribute[]{new Attribute.Default<String>("id", "branchRates")}
                );
                writeParameter("rate", "clock.rate", options, writer);
                writer.writeCloseTag(StrictClockBranchRates.STRICT_CLOCK_BRANCH_RATES);
                break;

            case UNCORRELATED_EXPONENTIAL:
            case UNCORRELATED_LOGNORMAL:
                writer.writeComment("The uncorrelated relaxed clock (Drummond, Ho, Phillips & Rambaut, 2006)");

                //if (options.isFixedSubstitutionRate()) {
                //    attributes = new Attribute[]{
                //            new Attribute.Default<String>("id", "branchRates"),
                //            new Attribute.Default<Double>(DiscretizedBranchRatesParser.NORMALIZED_MEAN, options.getMeanSubstitutionRate())
                //    };
                //} else {
                attributes = new Attribute[]{new Attribute.Default<String>("id", "branchRates")};
                //}
                writer.writeOpenTag(DiscretizedBranchRatesParser.DISCRETIZED_BRANCH_RATES, attributes);
                writer.writeTag(TreeModel.TREE_MODEL, new Attribute.Default<String>("idref", "treeModel"), true);
                writer.writeOpenTag("distribution");
                if (options.clockType == ClockType.UNCORRELATED_EXPONENTIAL) {
                    if (options.isFixedSubstitutionRate()) {

                        fixParameter("uced.mean", options.getMeanSubstitutionRate());
                    }

                    final String eModelName = ExponentialDistributionModel.EXPONENTIAL_DISTRIBUTION_MODEL;
                    writer.writeOpenTag(eModelName);
                    writeParameter("mean", "uced.mean", options, writer);
                    writer.writeCloseTag(eModelName);
                } else if (options.clockType == ClockType.UNCORRELATED_LOGNORMAL) {
                    if (options.isFixedSubstitutionRate()) {

                        fixParameter("ucld.mean", options.getMeanSubstitutionRate());
                    }

                    writer.writeOpenTag("logNormalDistributionModel", new Attribute.Default<String>(LogNormalDistributionModel.MEAN_IN_REAL_SPACE, "true"));
                    writeParameter("mean", "ucld.mean", options, writer);
                    writeParameter("stdev", "ucld.stdev", options, writer);
                    writer.writeCloseTag("logNormalDistributionModel");
                } else {
                    throw new RuntimeException("Unrecognised relaxed clock model");
                }
                writer.writeCloseTag("distribution");
                writer.writeOpenTag(DiscretizedBranchRatesParser.RATE_CATEGORIES);
                categoryCount = (options.taxonList.getTaxonCount() - 1) * 2;
                writeParameter("branchRates.categories", categoryCount, writer);
                writer.writeCloseTag(DiscretizedBranchRatesParser.RATE_CATEGORIES);
                writer.writeCloseTag(DiscretizedBranchRatesParser.DISCRETIZED_BRANCH_RATES);

                writer.writeText("");
                writer.writeOpenTag(
                        RateStatistic.RATE_STATISTIC,
                        new Attribute[]{
                                new Attribute.Default<String>("id", "meanRate"),
                                new Attribute.Default<String>("name", "meanRate"),
                                new Attribute.Default<String>("mode", "mean"),
                                new Attribute.Default<String>("internal", "true"),
                                new Attribute.Default<String>("external", "true")
                        }
                );
                writer.writeTag(TreeModel.TREE_MODEL, new Attribute.Default<String>("idref", "treeModel"), true);
                writer.writeTag(DiscretizedBranchRatesParser.DISCRETIZED_BRANCH_RATES, new Attribute.Default<String>("idref", "branchRates"), true);
                writer.writeCloseTag(RateStatistic.RATE_STATISTIC);

                writer.writeText("");
                writer.writeOpenTag(
                        RateStatistic.RATE_STATISTIC,
                        new Attribute[]{
                                new Attribute.Default<String>("id", "coefficientOfVariation"),
                                new Attribute.Default<String>("name", "coefficientOfVariation"),
                                new Attribute.Default<String>("mode", "coefficientOfVariation"),
                                new Attribute.Default<String>("internal", "true"),
                                new Attribute.Default<String>("external", "true")
                        }
                );
                writer.writeTag(TreeModel.TREE_MODEL, new Attribute.Default<String>("idref", "treeModel"), true);
                writer.writeTag(DiscretizedBranchRatesParser.DISCRETIZED_BRANCH_RATES, new Attribute.Default<String>("idref", "branchRates"), true);
                writer.writeCloseTag(RateStatistic.RATE_STATISTIC);

                writer.writeText("");
                writer.writeOpenTag(
                        RateCovarianceStatistic.RATE_COVARIANCE_STATISTIC,
                        new Attribute[]{
                                new Attribute.Default<String>("id", "covariance"),
                                new Attribute.Default<String>("name", "covariance")
                        }
                );
                writer.writeTag(TreeModel.TREE_MODEL, new Attribute.Default<String>("idref", "treeModel"), true);
                writer.writeTag(DiscretizedBranchRatesParser.DISCRETIZED_BRANCH_RATES, new Attribute.Default<String>("idref", "branchRates"), true);
                writer.writeCloseTag(RateCovarianceStatistic.RATE_COVARIANCE_STATISTIC);
                break;

            case AUTOCORRELATED_LOGNORMAL:
                writer.writeComment("The autocorrelated relaxed clock (Rannala & Yang, 2007)");

                attributes = new Attribute[]{
                        new Attribute.Default<String>("id", "branchRates"),
                        new Attribute.Default<String>("episodic", "false"),
                        new Attribute.Default<String>("logspace", "true"),
                };

                writer.writeOpenTag(ACLikelihood.AC_LIKELIHOOD, attributes);
                writer.writeTag(TreeModel.TREE_MODEL, new Attribute.Default<String>("idref", "treeModel"), true);

                if (options.isFixedSubstitutionRate()) {
                    fixParameter("treeModel.rootRate", options.getMeanSubstitutionRate());
                }

                writeParameterRef("rates", "treeModel.nodeRates", writer);
                writeParameterRef("rootRate", "treeModel.rootRate", writer);
                writeParameter("variance", "branchRates.var", options, writer);

                writer.writeCloseTag(ACLikelihood.AC_LIKELIHOOD);

                writer.writeText("");
                writer.writeOpenTag(
                        RateStatistic.RATE_STATISTIC,
                        new Attribute[]{
                                new Attribute.Default<String>("id", "meanRate"),
                                new Attribute.Default<String>("name", "meanRate"),
                                new Attribute.Default<String>("mode", "mean"),
                                new Attribute.Default<String>("internal", "true"),
                                new Attribute.Default<String>("external", "true")
                        }
                );
                writer.writeTag(TreeModel.TREE_MODEL, new Attribute.Default<String>("idref", "treeModel"), true);
                writer.writeTag(ACLikelihood.AC_LIKELIHOOD, new Attribute.Default<String>("idref", "branchRates"), true);
                writer.writeCloseTag(RateStatistic.RATE_STATISTIC);

                writer.writeText("");
                writer.writeOpenTag(
                        RateStatistic.RATE_STATISTIC,
                        new Attribute[]{
                                new Attribute.Default<String>("id", "coefficientOfVariation"),
                                new Attribute.Default<String>("name", "coefficientOfVariation"),
                                new Attribute.Default<String>("mode", "coefficientOfVariation"),
                                new Attribute.Default<String>("internal", "true"),
                                new Attribute.Default<String>("external", "true")
                        }
                );
                writer.writeTag(TreeModel.TREE_MODEL, new Attribute.Default<String>("idref", "treeModel"), true);
                writer.writeTag(ACLikelihood.AC_LIKELIHOOD, new Attribute.Default<String>("idref", "branchRates"), true);
                writer.writeCloseTag(RateStatistic.RATE_STATISTIC);

                writer.writeText("");
                writer.writeOpenTag(
                        RateCovarianceStatistic.RATE_COVARIANCE_STATISTIC,
                        new Attribute[]{
                                new Attribute.Default<String>("id", "covariance"),
                                new Attribute.Default<String>("name", "covariance")
                        }
                );
                writer.writeTag(TreeModel.TREE_MODEL, new Attribute.Default<String>("idref", "treeModel"), true);
                writer.writeTag(ACLikelihood.AC_LIKELIHOOD, new Attribute.Default<String>("idref", "branchRates"), true);
                writer.writeCloseTag(RateCovarianceStatistic.RATE_COVARIANCE_STATISTIC);

                break;

            case RANDOM_LOCAL_CLOCK:
                if (options.isFixedSubstitutionRate()) {

                    fixParameter("clock.rate", options.getMeanSubstitutionRate());
                }

                writer.writeComment("The random local clock model (Drummond & Suchard, 2007)");
                writer.writeOpenTag(
                        RandomLocalClockModel.LOCAL_BRANCH_RATES,
                        new Attribute[]{
                                new Attribute.Default<String>("id", "branchRates"),
                                new Attribute.Default<String>("ratesAreMultipliers", "false")
                        }
                );
                writer.writeTag(TreeModel.TREE_MODEL, new Attribute.Default<String>("idref", "treeModel"), true);

                writer.writeOpenTag("rates");
                writer.writeTag("parameter", new Attribute.Default<String>("idref", "localClock.rates"), true);
                writer.writeCloseTag("rates");

                writer.writeOpenTag("rateIndicator");
                writer.writeTag("parameter", new Attribute.Default<String>("idref", "localClock.changes"), true);
                writer.writeCloseTag("rateIndicator");

                writeParameter("clockRate", "clock.rate", options, writer);

                writer.writeCloseTag(RandomLocalClockModel.LOCAL_BRANCH_RATES);

                writer.writeText("");
                writer.writeOpenTag(
                        SumStatistic.SUM_STATISTIC,
                        new Attribute[]{
                                new Attribute.Default<String>("id", "rateChanges"),
                                new Attribute.Default<String>("name", "rateChangeCount"),
                                new Attribute.Default<String>("elementwise", "true"),
                        }
                );
                writer.writeTag("parameter", new Attribute.Default<String>("idref", "localClock.changes"), true);
                writer.writeCloseTag(SumStatistic.SUM_STATISTIC);

                writer.writeText("");

                writer.writeOpenTag(
                        RateStatistic.RATE_STATISTIC,
                        new Attribute[]{
                                new Attribute.Default<String>("id", "meanRate"),
                                new Attribute.Default<String>("name", "meanRate"),
                                new Attribute.Default<String>("mode", "mean"),
                                new Attribute.Default<String>("internal", "true"),
                                new Attribute.Default<String>("external", "true")
                        }
                );
                writer.writeTag(TreeModel.TREE_MODEL, new Attribute.Default<String>("idref", "treeModel"), true);
                writer.writeTag(RandomLocalClockModel.LOCAL_BRANCH_RATES, new Attribute.Default<String>("idref", "branchRates"), true);
                writer.writeCloseTag(RateStatistic.RATE_STATISTIC);

                writer.writeText("");
                writer.writeOpenTag(
                        RateStatistic.RATE_STATISTIC,
                        new Attribute[]{
                                new Attribute.Default<String>("id", "coefficientOfVariation"),
                                new Attribute.Default<String>("name", "coefficientOfVariation"),
                                new Attribute.Default<String>("mode", "coefficientOfVariation"),
                                new Attribute.Default<String>("internal", "true"),
                                new Attribute.Default<String>("external", "true")
                        }
                );
                writer.writeTag(TreeModel.TREE_MODEL, new Attribute.Default<String>("idref", "treeModel"), true);
                writer.writeTag(RandomLocalClockModel.LOCAL_BRANCH_RATES, new Attribute.Default<String>("idref", "branchRates"), true);
                writer.writeCloseTag(RateStatistic.RATE_STATISTIC);

                writer.writeText("");
                writer.writeOpenTag(
                        RateCovarianceStatistic.RATE_COVARIANCE_STATISTIC,
                        new Attribute[]{
                                new Attribute.Default<String>("id", "covariance"),
                                new Attribute.Default<String>("name", "covariance")
                        }
                );
                writer.writeTag(TreeModel.TREE_MODEL, new Attribute.Default<String>("idref", "treeModel"), true);
                writer.writeTag(RandomLocalClockModel.LOCAL_BRANCH_RATES, new Attribute.Default<String>("idref", "branchRates"), true);
                writer.writeCloseTag(RateCovarianceStatistic.RATE_COVARIANCE_STATISTIC);
                break;

            default:
                throw new IllegalArgumentException("Unknown clock model");
        }

    }

    public void writeBranchRatesModelOLD(XMLWriter writer) {
        if (options.clockType == ClockType.STRICT_CLOCK) {
            if (options.isFixedSubstitutionRate()) {

                fixParameter("clock.rate", options.getMeanSubstitutionRate());
            }

            writer.writeComment("The strict clock (Uniform rates across branches)");
            writer.writeOpenTag(
                    StrictClockBranchRates.STRICT_CLOCK_BRANCH_RATES,
                    new Attribute[]{new Attribute.Default<String>("id", "branchRates")}
            );
            writeParameter("rate", "clock.rate", options, writer);
            writer.writeCloseTag(StrictClockBranchRates.STRICT_CLOCK_BRANCH_RATES);
        } else if (options.clockType == ClockType.RANDOM_LOCAL_CLOCK) {
            if (options.isFixedSubstitutionRate()) {

                fixParameter("clock.rate", options.getMeanSubstitutionRate());
            }

            writer.writeComment("The random local clock model (Drummond & Suchard, 2007)");
            writer.writeOpenTag(
                    RandomLocalClockModel.LOCAL_BRANCH_RATES,
                    new Attribute[]{
                            new Attribute.Default<String>("id", "branchRates"),
                            new Attribute.Default<String>("ratesAreMultipliers", "false")
                    }
            );
            writer.writeTag(TreeModel.TREE_MODEL, new Attribute.Default<String>("idref", "treeModel"), true);

            writer.writeOpenTag("rates");
            writer.writeTag("parameter", new Attribute.Default<String>("idref", "localClock.rates"), true);
            writer.writeCloseTag("rates");

            writer.writeOpenTag("rateIndicator");
            writer.writeTag("parameter", new Attribute.Default<String>("idref", "localClock.changes"), true);
            writer.writeCloseTag("rateIndicator");

            writeParameter("clockRate", "clock.rate", options, writer);

            writer.writeCloseTag(RandomLocalClockModel.LOCAL_BRANCH_RATES);

            writer.writeText("");
            writer.writeOpenTag(
                    SumStatistic.SUM_STATISTIC,
                    new Attribute[]{
                            new Attribute.Default<String>("id", "rateChanges"),
                            new Attribute.Default<String>("name", "rateChangeCount"),
                            new Attribute.Default<String>("elementwise", "true"),
                    }
            );
            writer.writeTag("parameter", new Attribute.Default<String>("idref", "localClock.changes"), true);
            writer.writeCloseTag(SumStatistic.SUM_STATISTIC);

            writer.writeText("");

            writer.writeOpenTag(
                    RateStatistic.RATE_STATISTIC,
                    new Attribute[]{
                            new Attribute.Default<String>("id", "meanRate"),
                            new Attribute.Default<String>("name", "meanRate"),
                            new Attribute.Default<String>("mode", "mean"),
                            new Attribute.Default<String>("internal", "true"),
                            new Attribute.Default<String>("external", "true")
                    }
            );
            writer.writeTag(TreeModel.TREE_MODEL, new Attribute.Default<String>("idref", "treeModel"), true);
            writer.writeTag(RandomLocalClockModel.LOCAL_BRANCH_RATES, new Attribute.Default<String>("idref", "branchRates"), true);
            writer.writeCloseTag(RateStatistic.RATE_STATISTIC);

            writer.writeText("");
            writer.writeOpenTag(
                    RateStatistic.RATE_STATISTIC,
                    new Attribute[]{
                            new Attribute.Default<String>("id", "coefficientOfVariation"),
                            new Attribute.Default<String>("name", "coefficientOfVariation"),
                            new Attribute.Default<String>("mode", "coefficientOfVariation"),
                            new Attribute.Default<String>("internal", "true"),
                            new Attribute.Default<String>("external", "true")
                    }
            );
            writer.writeTag(TreeModel.TREE_MODEL, new Attribute.Default<String>("idref", "treeModel"), true);
            writer.writeTag(RandomLocalClockModel.LOCAL_BRANCH_RATES, new Attribute.Default<String>("idref", "branchRates"), true);
            writer.writeCloseTag(RateStatistic.RATE_STATISTIC);

            writer.writeText("");
            writer.writeOpenTag(
                    RateCovarianceStatistic.RATE_COVARIANCE_STATISTIC,
                    new Attribute[]{
                            new Attribute.Default<String>("id", "covariance"),
                            new Attribute.Default<String>("name", "covariance")
                    }
            );
            writer.writeTag(TreeModel.TREE_MODEL, new Attribute.Default<String>("idref", "treeModel"), true);
            writer.writeTag(RandomLocalClockModel.LOCAL_BRANCH_RATES, new Attribute.Default<String>("idref", "branchRates"), true);
            writer.writeCloseTag(RateCovarianceStatistic.RATE_COVARIANCE_STATISTIC);

        } else {
            writer.writeComment("The uncorrelated relaxed clock (Drummond, Ho, Phillips & Rambaut, 2006)");
            Attribute[] attributes;
            //if (options.isFixedSubstitutionRate()) {
            //    attributes = new Attribute[]{
            //            new Attribute.Default<String>("id", "branchRates"),
            //            new Attribute.Default<Double>(DiscretizedBranchRatesParser.NORMALIZED_MEAN, options.getMeanSubstitutionRate())
            //    };
            //} else {
            attributes = new Attribute[]{new Attribute.Default<String>("id", "branchRates")};
            //}
            writer.writeOpenTag(DiscretizedBranchRatesParser.DISCRETIZED_BRANCH_RATES, attributes);
            writer.writeTag(TreeModel.TREE_MODEL, new Attribute.Default<String>("idref", "treeModel"), true);
            writer.writeOpenTag("distribution");
            if (options.clockType == ClockType.UNCORRELATED_EXPONENTIAL) {
                if (options.isFixedSubstitutionRate()) {

                    fixParameter("uced.mean", options.getMeanSubstitutionRate());
                }

                final String eModelName = ExponentialDistributionModel.EXPONENTIAL_DISTRIBUTION_MODEL;
                writer.writeOpenTag(eModelName);
                writeParameter("mean", "uced.mean", options, writer);
                writer.writeCloseTag(eModelName);
            } else if (options.clockType == ClockType.UNCORRELATED_LOGNORMAL) {
                if (options.isFixedSubstitutionRate()) {

                    fixParameter("ucld.mean", options.getMeanSubstitutionRate());
                }

                writer.writeOpenTag("logNormalDistributionModel", new Attribute.Default<String>(LogNormalDistributionModel.MEAN_IN_REAL_SPACE, "true"));
                writeParameter("mean", "ucld.mean", options, writer);
                writeParameter("stdev", "ucld.stdev", options, writer);
                writer.writeCloseTag("logNormalDistributionModel");
            } else {
                throw new RuntimeException("Unrecognised relaxed clock model");
            }
            writer.writeCloseTag("distribution");
            writer.writeOpenTag(DiscretizedBranchRatesParser.RATE_CATEGORIES);
            int categoryCount = (options.taxonList.getTaxonCount() - 1) * 2;
            writeParameter("branchRates.categories", categoryCount, writer);
            writer.writeCloseTag(DiscretizedBranchRatesParser.RATE_CATEGORIES);
            writer.writeCloseTag(DiscretizedBranchRatesParser.DISCRETIZED_BRANCH_RATES);

            writer.writeText("");
            writer.writeOpenTag(
                    RateStatistic.RATE_STATISTIC,
                    new Attribute[]{
                            new Attribute.Default<String>("id", "meanRate"),
                            new Attribute.Default<String>("name", "meanRate"),
                            new Attribute.Default<String>("mode", "mean"),
                            new Attribute.Default<String>("internal", "true"),
                            new Attribute.Default<String>("external", "true")
                    }
            );
            writer.writeTag(TreeModel.TREE_MODEL, new Attribute.Default<String>("idref", "treeModel"), true);
            writer.writeTag(DiscretizedBranchRatesParser.DISCRETIZED_BRANCH_RATES, new Attribute.Default<String>("idref", "branchRates"), true);
            writer.writeCloseTag(RateStatistic.RATE_STATISTIC);

            writer.writeText("");
            writer.writeOpenTag(
                    RateStatistic.RATE_STATISTIC,
                    new Attribute[]{
                            new Attribute.Default<String>("id", "coefficientOfVariation"),
                            new Attribute.Default<String>("name", "coefficientOfVariation"),
                            new Attribute.Default<String>("mode", "coefficientOfVariation"),
                            new Attribute.Default<String>("internal", "true"),
                            new Attribute.Default<String>("external", "true")
                    }
            );
            writer.writeTag(TreeModel.TREE_MODEL, new Attribute.Default<String>("idref", "treeModel"), true);
            writer.writeTag(DiscretizedBranchRatesParser.DISCRETIZED_BRANCH_RATES, new Attribute.Default<String>("idref", "branchRates"), true);
            writer.writeCloseTag(RateStatistic.RATE_STATISTIC);

            writer.writeText("");
            writer.writeOpenTag(
                    RateCovarianceStatistic.RATE_COVARIANCE_STATISTIC,
                    new Attribute[]{
                            new Attribute.Default<String>("id", "covariance"),
                            new Attribute.Default<String>("name", "covariance")
                    }
            );
            writer.writeTag(TreeModel.TREE_MODEL, new Attribute.Default<String>("idref", "treeModel"), true);
            writer.writeTag(DiscretizedBranchRatesParser.DISCRETIZED_BRANCH_RATES, new Attribute.Default<String>("idref", "branchRates"), true);
            writer.writeCloseTag(RateCovarianceStatistic.RATE_COVARIANCE_STATISTIC);
        }
    }

}
