package dr.app.beauti.generator;

import dr.app.beauti.util.XMLWriter;
import dr.app.beauti.components.ComponentFactory;
import dr.app.beauti.options.BeautiOptions;
import dr.app.beauti.options.ClockType;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.branchratemodel.RandomLocalClockModel;
import dr.evomodel.branchratemodel.StrictClockBranchRates;
import dr.evomodel.clock.ACLikelihood;
import dr.evomodel.clock.RateEvolutionLikelihood;
import dr.evomodel.tree.RateCovarianceStatistic;
import dr.evomodel.tree.RateStatistic;
import dr.evomodel.tree.TreeModel;
import dr.evomodelxml.DiscretizedBranchRatesParser;
import dr.inference.distribution.ExponentialDistributionModel;
import dr.inference.distribution.LogNormalDistributionModel;
import dr.inference.model.SumStatistic;
import dr.util.Attribute;
import dr.xml.XMLParser;

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
                        new Attribute[]{new Attribute.Default<String>(XMLParser.ID, genePrefix + BranchRateModel.BRANCH_RATES)}
                );
                writeParameter("rate", "clock.rate", options, writer);
                writer.writeCloseTag(StrictClockBranchRates.STRICT_CLOCK_BRANCH_RATES);
                break;

            case UNCORRELATED_EXPONENTIAL:
            case UNCORRELATED_LOGNORMAL:
                writer.writeComment("The uncorrelated relaxed clock (Drummond, Ho, Phillips & Rambaut, 2006)");

                //if (options.isFixedSubstitutionRate()) {
                //    attributes = new Attribute[]{
                //            new Attribute.Default<String>(XMLParser.ID, BranchRateModel.BRANCH_RATES),
                //            new Attribute.Default<Double>(DiscretizedBranchRatesParser.NORMALIZED_MEAN, options.getMeanSubstitutionRate())
                //    };
                //} else {
                attributes = new Attribute[]{new Attribute.Default<String>(XMLParser.ID, genePrefix + BranchRateModel.BRANCH_RATES)};
                //}
                writer.writeOpenTag(DiscretizedBranchRatesParser.DISCRETIZED_BRANCH_RATES, attributes);
                writer.writeIDref(TreeModel.TREE_MODEL, genePrefix + TreeModel.TREE_MODEL);
                writer.writeOpenTag("distribution");
                if (options.clockType == ClockType.UNCORRELATED_EXPONENTIAL) {
                    if (options.isFixedSubstitutionRate()) {

                        fixParameter(ClockType.UCED_MEAN, options.getMeanSubstitutionRate());
                    }

                    final String eModelName = ExponentialDistributionModel.EXPONENTIAL_DISTRIBUTION_MODEL;
                    writer.writeOpenTag(eModelName);
                    writeParameter("mean", ClockType.UCED_MEAN, options, writer);
                    writer.writeCloseTag(eModelName);
                } else if (options.clockType == ClockType.UNCORRELATED_LOGNORMAL) {
                    if (options.isFixedSubstitutionRate()) {

                        fixParameter(ClockType.UCLD_MEAN, options.getMeanSubstitutionRate());
                    }

                    writer.writeOpenTag("logNormalDistributionModel",
                            new Attribute.Default<String>(genePrefix + LogNormalDistributionModel.MEAN_IN_REAL_SPACE, "true"));
                    writeParameter("mean", ClockType.UCLD_MEAN, options, writer);
                    writeParameter("stdev", ClockType.UCLD_STDEV, options, writer);
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
                                new Attribute.Default<String>(XMLParser.ID, genePrefix + "meanRate"),
                                new Attribute.Default<String>("name", "meanRate"),
                                new Attribute.Default<String>("mode", "mean"),
                                new Attribute.Default<String>("internal", "true"),
                                new Attribute.Default<String>("external", "true")
                        }
                );
                writer.writeIDref(TreeModel.TREE_MODEL, genePrefix + TreeModel.TREE_MODEL);
                writer.writeIDref(DiscretizedBranchRatesParser.DISCRETIZED_BRANCH_RATES, genePrefix + BranchRateModel.BRANCH_RATES);
                writer.writeCloseTag(RateStatistic.RATE_STATISTIC);

                writer.writeText("");
                writer.writeOpenTag(
                        RateStatistic.RATE_STATISTIC,
                        new Attribute[]{
                                new Attribute.Default<String>(XMLParser.ID, genePrefix + RateStatistic.COEFFICIENT_OF_VARIATION),
                                new Attribute.Default<String>("name", RateStatistic.COEFFICIENT_OF_VARIATION),
                                new Attribute.Default<String>("mode", RateStatistic.COEFFICIENT_OF_VARIATION),
                                new Attribute.Default<String>("internal", "true"),
                                new Attribute.Default<String>("external", "true")
                        }
                );
                writer.writeIDref(TreeModel.TREE_MODEL, genePrefix + TreeModel.TREE_MODEL);
                writer.writeIDref(DiscretizedBranchRatesParser.DISCRETIZED_BRANCH_RATES, genePrefix + BranchRateModel.BRANCH_RATES);
                writer.writeCloseTag(RateStatistic.RATE_STATISTIC);

                writer.writeText("");
                writer.writeOpenTag(
                        RateCovarianceStatistic.RATE_COVARIANCE_STATISTIC,
                        new Attribute[]{
                                new Attribute.Default<String>(XMLParser.ID, genePrefix + "covariance"),
                                new Attribute.Default<String>("name", "covariance")
                        }
                );
                writer.writeIDref(TreeModel.TREE_MODEL, genePrefix + TreeModel.TREE_MODEL);
                writer.writeIDref(DiscretizedBranchRatesParser.DISCRETIZED_BRANCH_RATES, genePrefix + BranchRateModel.BRANCH_RATES);
                writer.writeCloseTag(RateCovarianceStatistic.RATE_COVARIANCE_STATISTIC);
                break;

            case AUTOCORRELATED_LOGNORMAL:
                writer.writeComment("The autocorrelated relaxed clock (Rannala & Yang, 2007)");

                attributes = new Attribute[]{
                        new Attribute.Default<String>(XMLParser.ID, genePrefix + BranchRateModel.BRANCH_RATES),
                        new Attribute.Default<String>("episodic", "false"),
                        new Attribute.Default<String>("logspace", "true"),
                };

                writer.writeOpenTag(ACLikelihood.AC_LIKELIHOOD, attributes);
                writer.writeIDref(TreeModel.TREE_MODEL, genePrefix + TreeModel.TREE_MODEL);

                if (options.isFixedSubstitutionRate()) {
                    fixParameter("treeModel.rootRate", options.getMeanSubstitutionRate());
                }

                writeParameterRef("rates", "treeModel.nodeRates", writer);
                writeParameterRef(RateEvolutionLikelihood.ROOTRATE, "treeModel.rootRate", writer);
                writeParameter("variance", "branchRates.var", options, writer);

                writer.writeCloseTag(ACLikelihood.AC_LIKELIHOOD);

                writer.writeText("");
                writer.writeOpenTag(
                        RateStatistic.RATE_STATISTIC,
                        new Attribute[]{
                                new Attribute.Default<String>(XMLParser.ID, genePrefix + "meanRate"),
                                new Attribute.Default<String>("name", "meanRate"),
                                new Attribute.Default<String>("mode", "mean"),
                                new Attribute.Default<String>("internal", "true"),
                                new Attribute.Default<String>("external", "true")
                        }
                );
                writer.writeIDref(TreeModel.TREE_MODEL, genePrefix + TreeModel.TREE_MODEL);
                writer.writeIDref(ACLikelihood.AC_LIKELIHOOD, genePrefix + BranchRateModel.BRANCH_RATES);
                writer.writeCloseTag(RateStatistic.RATE_STATISTIC);

                writer.writeText("");
                writer.writeOpenTag(
                        RateStatistic.RATE_STATISTIC,
                        new Attribute[]{
                                new Attribute.Default<String>(XMLParser.ID, genePrefix + RateStatistic.COEFFICIENT_OF_VARIATION),
                                new Attribute.Default<String>("name", RateStatistic.COEFFICIENT_OF_VARIATION),
                                new Attribute.Default<String>("mode", RateStatistic.COEFFICIENT_OF_VARIATION),
                                new Attribute.Default<String>("internal", "true"),
                                new Attribute.Default<String>("external", "true")
                        }
                );
                writer.writeIDref(TreeModel.TREE_MODEL, genePrefix + TreeModel.TREE_MODEL);
                writer.writeIDref(ACLikelihood.AC_LIKELIHOOD, genePrefix + BranchRateModel.BRANCH_RATES);
                writer.writeCloseTag(RateStatistic.RATE_STATISTIC);

                writer.writeText("");
                writer.writeOpenTag(
                        RateCovarianceStatistic.RATE_COVARIANCE_STATISTIC,
                        new Attribute[]{
                                new Attribute.Default<String>(XMLParser.ID, genePrefix + "covariance"),
                                new Attribute.Default<String>("name", "covariance")
                        }
                );
                writer.writeIDref(TreeModel.TREE_MODEL, genePrefix + TreeModel.TREE_MODEL);
                writer.writeIDref(ACLikelihood.AC_LIKELIHOOD, genePrefix + BranchRateModel.BRANCH_RATES);
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
                                new Attribute.Default<String>(XMLParser.ID, genePrefix + BranchRateModel.BRANCH_RATES),
                                new Attribute.Default<String>("ratesAreMultipliers", "false")
                        }
                );
                writer.writeIDref(TreeModel.TREE_MODEL, genePrefix + TreeModel.TREE_MODEL);

                writer.writeOpenTag("rates");
                writer.writeIDref("parameter", genePrefix + ClockType.LOCAL_CLOCK + "rates");
                writer.writeCloseTag("rates");

                writer.writeOpenTag("rateIndicator");
                writer.writeIDref("parameter", genePrefix + ClockType.LOCAL_CLOCK + "changes");
                writer.writeCloseTag("rateIndicator");

                writeParameter("clockRate", "clock.rate", options, writer);

                writer.writeCloseTag(RandomLocalClockModel.LOCAL_BRANCH_RATES);

                writer.writeText("");
                writer.writeOpenTag(
                        SumStatistic.SUM_STATISTIC,
                        new Attribute[]{
                                new Attribute.Default<String>(XMLParser.ID, genePrefix + "rateChanges"),
                                new Attribute.Default<String>("name", "rateChangeCount"),
                                new Attribute.Default<String>("elementwise", "true"),
                        }
                );
                writer.writeIDref("parameter", genePrefix + ClockType.LOCAL_CLOCK + "changes");
                writer.writeCloseTag(SumStatistic.SUM_STATISTIC);

                writer.writeText("");

                writer.writeOpenTag(
                        RateStatistic.RATE_STATISTIC,
                        new Attribute[]{
                                new Attribute.Default<String>(XMLParser.ID, genePrefix + "meanRate"),
                                new Attribute.Default<String>("name", "meanRate"),
                                new Attribute.Default<String>("mode", "mean"),
                                new Attribute.Default<String>("internal", "true"),
                                new Attribute.Default<String>("external", "true")
                        }
                );
                writer.writeIDref(TreeModel.TREE_MODEL, genePrefix + TreeModel.TREE_MODEL);
                writer.writeIDref(RandomLocalClockModel.LOCAL_BRANCH_RATES, genePrefix + BranchRateModel.BRANCH_RATES);
                writer.writeCloseTag(RateStatistic.RATE_STATISTIC);

                writer.writeText("");
                writer.writeOpenTag(
                        RateStatistic.RATE_STATISTIC,
                        new Attribute[]{
                                new Attribute.Default<String>(XMLParser.ID, genePrefix + RateStatistic.COEFFICIENT_OF_VARIATION),
                                new Attribute.Default<String>("name", RateStatistic.COEFFICIENT_OF_VARIATION),
                                new Attribute.Default<String>("mode", RateStatistic.COEFFICIENT_OF_VARIATION),
                                new Attribute.Default<String>("internal", "true"),
                                new Attribute.Default<String>("external", "true")
                        }
                );
                writer.writeIDref(TreeModel.TREE_MODEL, genePrefix + TreeModel.TREE_MODEL);
                writer.writeIDref(RandomLocalClockModel.LOCAL_BRANCH_RATES, genePrefix + BranchRateModel.BRANCH_RATES);
                writer.writeCloseTag(RateStatistic.RATE_STATISTIC);

                writer.writeText("");
                writer.writeOpenTag(
                        RateCovarianceStatistic.RATE_COVARIANCE_STATISTIC,
                        new Attribute[]{
                                new Attribute.Default<String>(XMLParser.ID, genePrefix + "covariance"),
                                new Attribute.Default<String>("name", "covariance")
                        }
                );
                writer.writeIDref(TreeModel.TREE_MODEL, genePrefix + TreeModel.TREE_MODEL);
                writer.writeIDref(RandomLocalClockModel.LOCAL_BRANCH_RATES, genePrefix + BranchRateModel.BRANCH_RATES);
                writer.writeCloseTag(RateCovarianceStatistic.RATE_COVARIANCE_STATISTIC);
                break;

            default:
                throw new IllegalArgumentException("Unknown clock model");
        }

    }
/*
    public void writeBranchRatesModelOLD(XMLWriter writer) {
        if (options.clockType == ClockType.STRICT_CLOCK) {
            if (options.isFixedSubstitutionRate()) {

                fixParameter("clock.rate", options.getMeanSubstitutionRate());
            }

            writer.writeComment("The strict clock (Uniform rates across branches)");
            writer.writeOpenTag(
                    StrictClockBranchRates.STRICT_CLOCK_BRANCH_RATES,
                    new Attribute[]{new Attribute.Default<String>(XMLParser.ID, BranchRateModel.BRANCH_RATES)}
            );
            writeParameter("rate", "clock.rate", options, writer);
            writer.writeCloseTag(StrictClockBranchRates.STRICT_CLOCK_BRANCH_RATES);
        } else if( options.clockType == ClockType.RANDOM_LOCAL_CLOCK ) {
            if (options.isFixedSubstitutionRate()) {

                fixParameter("clock.rate", options.getMeanSubstitutionRate());
            }

            writer.writeComment("The random local clock model (Drummond & Suchard, 2007)");
            writer.writeOpenTag(
                    RandomLocalClockModel.LOCAL_BRANCH_RATES,
                    new Attribute[]{
                            new Attribute.Default<String>(XMLParser.ID, BranchRateModel.BRANCH_RATES),
                            new Attribute.Default<String>("ratesAreMultipliers", "false")
                    }
            );
            writer.writeIDref(TreeModel.TREE_MODEL, TreeModel.TREE_MODEL);

            writer.writeOpenTag("rates");
            writer.writeIDref("parameter", ClockType.LOCAL_CLOCK + "rates");
            writer.writeCloseTag("rates");

            writer.writeOpenTag("rateIndicator");
            writer.writeIDref("parameter", ClockType.LOCAL_CLOCK + "changes");
            writer.writeCloseTag("rateIndicator");

            writeParameter("clockRate", "clock.rate", options, writer);

            writer.writeCloseTag(RandomLocalClockModel.LOCAL_BRANCH_RATES);

            writer.writeText("");
            writer.writeOpenTag(
                    SumStatistic.SUM_STATISTIC,
                    new Attribute[]{
                            new Attribute.Default<String>(XMLParser.ID, "rateChanges"),
                            new Attribute.Default<String>("name", "rateChangeCount"),
                            new Attribute.Default<String>("elementwise", "true"),
                    }
            );
            writer.writeIDref("parameter", ClockType.LOCAL_CLOCK + "changes");
            writer.writeCloseTag(SumStatistic.SUM_STATISTIC);

            writer.writeText("");

            writer.writeOpenTag(
                    RateStatistic.RATE_STATISTIC,
                    new Attribute[]{
                            new Attribute.Default<String>(XMLParser.ID, "meanRate"),
                            new Attribute.Default<String>("name", "meanRate"),
                            new Attribute.Default<String>("mode", "mean"),
                            new Attribute.Default<String>("internal", "true"),
                            new Attribute.Default<String>("external", "true")
                    }
            );
            writer.writeIDref(TreeModel.TREE_MODEL, TreeModel.TREE_MODEL);
            writer.writeIDref(RandomLocalClockModel.LOCAL_BRANCH_RATES, BranchRateModel.BRANCH_RATES);
            writer.writeCloseTag(RateStatistic.RATE_STATISTIC);

            writer.writeText("");
            writer.writeOpenTag(
                    RateStatistic.RATE_STATISTIC,
                    new Attribute[]{
                            new Attribute.Default<String>(XMLParser.ID, RateStatistic.COEFFICIENT_OF_VARIATION),
                            new Attribute.Default<String>("name", RateStatistic.COEFFICIENT_OF_VARIATION),
                            new Attribute.Default<String>("mode", RateStatistic.COEFFICIENT_OF_VARIATION),
                            new Attribute.Default<String>("internal", "true"),
                            new Attribute.Default<String>("external", "true")
                    }
            );
            writer.writeIDref(TreeModel.TREE_MODEL, TreeModel.TREE_MODEL);
            writer.writeIDref(RandomLocalClockModel.LOCAL_BRANCH_RATES, BranchRateModel.BRANCH_RATES);
            writer.writeCloseTag(RateStatistic.RATE_STATISTIC);

            writer.writeText("");
            writer.writeOpenTag(
                    RateCovarianceStatistic.RATE_COVARIANCE_STATISTIC,
                    new Attribute[]{
                            new Attribute.Default<String>(XMLParser.ID, "covariance"),
                            new Attribute.Default<String>("name", "covariance")
                    }
            );
            writer.writeIDref(TreeModel.TREE_MODEL, TreeModel.TREE_MODEL);
            writer.writeIDref(RandomLocalClockModel.LOCAL_BRANCH_RATES, BranchRateModel.BRANCH_RATES);
            writer.writeCloseTag(RateCovarianceStatistic.RATE_COVARIANCE_STATISTIC);

        } else {
            writer.writeComment("The uncorrelated relaxed clock (Drummond, Ho, Phillips & Rambaut, 2006)");
            Attribute[] attributes;
            //if (options.isFixedSubstitutionRate()) {
            //    attributes = new Attribute[]{
            //            new Attribute.Default<String>(XMLParser.ID, BranchRateModel.BRANCH_RATES),
            //            new Attribute.Default<Double>(DiscretizedBranchRatesParser.NORMALIZED_MEAN, options.getMeanSubstitutionRate())
            //    };
            //} else {
            attributes = new Attribute[]{new Attribute.Default<String>(XMLParser.ID, BranchRateModel.BRANCH_RATES)};
            //}
            writer.writeOpenTag(DiscretizedBranchRatesParser.DISCRETIZED_BRANCH_RATES, attributes);
            writer.writeIDref(TreeModel.TREE_MODEL, TreeModel.TREE_MODEL);
            writer.writeOpenTag("distribution");
            if (options.clockType == ClockType.UNCORRELATED_EXPONENTIAL) {
                if (options.isFixedSubstitutionRate()) {

                    fixParameter(ClockType.UCED_MEAN, options.getMeanSubstitutionRate());
                }

                final String eModelName = ExponentialDistributionModel.EXPONENTIAL_DISTRIBUTION_MODEL;
                writer.writeOpenTag(eModelName);
                writeParameter("mean", ClockType.UCED_MEAN, options, writer);
                writer.writeCloseTag(eModelName);
            } else if (options.clockType == ClockType.UNCORRELATED_LOGNORMAL) {
                if (options.isFixedSubstitutionRate()) {

                    fixParameter(ClockType.UCLD_MEAN, options.getMeanSubstitutionRate());
                }

                writer.writeOpenTag("logNormalDistributionModel", new Attribute.Default<String>(LogNormalDistributionModel.MEAN_IN_REAL_SPACE, "true"));
                writeParameter("mean", ClockType.UCLD_MEAN, options, writer);
                writeParameter("stdev", ClockType.UCLD_STDEV, options, writer);
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
                            new Attribute.Default<String>(XMLParser.ID, "meanRate"),
                            new Attribute.Default<String>("name", "meanRate"),
                            new Attribute.Default<String>("mode", "mean"),
                            new Attribute.Default<String>("internal", "true"),
                            new Attribute.Default<String>("external", "true")
                    }
            );
            writer.writeIDref(TreeModel.TREE_MODEL, TreeModel.TREE_MODEL);
            writer.writeIDref(DiscretizedBranchRatesParser.DISCRETIZED_BRANCH_RATES, BranchRateModel.BRANCH_RATES);
            writer.writeCloseTag(RateStatistic.RATE_STATISTIC);

            writer.writeText("");
            writer.writeOpenTag(
                    RateStatistic.RATE_STATISTIC,
                    new Attribute[]{
                            new Attribute.Default<String>(XMLParser.ID, RateStatistic.COEFFICIENT_OF_VARIATION),
                            new Attribute.Default<String>("name", RateStatistic.COEFFICIENT_OF_VARIATION),
                            new Attribute.Default<String>("mode", RateStatistic.COEFFICIENT_OF_VARIATION),
                            new Attribute.Default<String>("internal", "true"),
                            new Attribute.Default<String>("external", "true")
                    }
            );
            writer.writeIDref(TreeModel.TREE_MODEL, TreeModel.TREE_MODEL);
            writer.writeIDref(DiscretizedBranchRatesParser.DISCRETIZED_BRANCH_RATES, BranchRateModel.BRANCH_RATES);
            writer.writeCloseTag(RateStatistic.RATE_STATISTIC);

            writer.writeText("");
            writer.writeOpenTag(
                    RateCovarianceStatistic.RATE_COVARIANCE_STATISTIC,
                    new Attribute[]{
                            new Attribute.Default<String>(XMLParser.ID, "covariance"),
                            new Attribute.Default<String>("name", "covariance")
                    }
            );
            writer.writeIDref(TreeModel.TREE_MODEL, TreeModel.TREE_MODEL);
            writer.writeIDref(DiscretizedBranchRatesParser.DISCRETIZED_BRANCH_RATES, BranchRateModel.BRANCH_RATES);
            writer.writeCloseTag(RateCovarianceStatistic.RATE_COVARIANCE_STATISTIC);
        }
    }
*/
}
