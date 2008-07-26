package dr.app.beauti.generator;

import dr.app.beauti.XMLWriter;
import dr.app.beauti.options.BeautiOptions;
import dr.app.beauti.options.ModelOptions;
import dr.evomodel.branchratemodel.DiscretizedBranchRates;
import dr.evomodel.branchratemodel.RandomLocalClockModel;
import dr.evomodel.branchratemodel.StrictClockBranchRates;
import dr.evomodel.tree.RateCovarianceStatistic;
import dr.evomodel.tree.RateStatistic;
import dr.evomodel.tree.TreeModel;
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
        if (options.clockModel == ModelOptions.STRICT_CLOCK) {
            if (options.fixedSubstitutionRate) {

                // TODO
                //fixParameter("clock.rate", options.meanSubstitutionRate);
            }

            writer.writeComment("The strict clock (Uniform rates across branches)");
            writer.writeOpenTag(
                    StrictClockBranchRates.STRICT_CLOCK_BRANCH_RATES,
                    new Attribute[]{new Attribute.Default<String>("id", "branchRates")}
            );
            writeParameter("rate", "clock.rate", writer, options);
            writer.writeCloseTag(StrictClockBranchRates.STRICT_CLOCK_BRANCH_RATES);
        } else if (options.clockModel == ModelOptions.RANDOM_LOCAL_CLOCK) {
            if (options.fixedSubstitutionRate) {

                // TODO
                //fixParameter("clock.rate", meanSubstitutionRate);
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

            writeParameter("clockRate", "clock.rate", writer, options);

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
            writer.writeOpenTag(
                    DiscretizedBranchRates.DISCRETIZED_BRANCH_RATES,
                    new Attribute[]{new Attribute.Default<String>("id", "branchRates")}
            );
            writer.writeTag(TreeModel.TREE_MODEL, new Attribute.Default<String>("idref", "treeModel"), true);
            writer.writeOpenTag("distribution");
            if (options.clockModel == ModelOptions.UNCORRELATED_EXPONENTIAL) {
                if (options.fixedSubstitutionRate) {

                    // TODO
                    //fixParameter("uced.mean", meanSubstitutionRate);
                }

                final String eModelName = ExponentialDistributionModel.EXPONENTIAL_DISTRIBUTION_MODEL;
                writer.writeOpenTag(eModelName);
                writeParameter("mean", "uced.mean", writer, options);
                writer.writeCloseTag(eModelName);
            } else if (options.clockModel == ModelOptions.UNCORRELATED_LOGNORMAL) {
                if (options.fixedSubstitutionRate) {

                    // TODO
                    //fixParameter("ucld.mean", meanSubstitutionRate);
                }

                writer.writeOpenTag("logNormalDistributionModel", new Attribute.Default<String>(LogNormalDistributionModel.MEAN_IN_REAL_SPACE, "true"));
                writeParameter("mean", "ucld.mean", writer, options);
                writeParameter("stdev", "ucld.stdev", writer, options);
                writer.writeCloseTag("logNormalDistributionModel");
            } else {
                throw new RuntimeException("Unrecognised relaxed clock model");
            }
            writer.writeCloseTag("distribution");
            writer.writeOpenTag("rateCategories");
            int categoryCount = (options.taxonList.getTaxonCount() - 1) * 2;
            writeParameter("branchRates.categories", categoryCount, writer);
            writer.writeCloseTag("rateCategories");
            writer.writeCloseTag(DiscretizedBranchRates.DISCRETIZED_BRANCH_RATES);

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
            writer.writeTag(DiscretizedBranchRates.DISCRETIZED_BRANCH_RATES, new Attribute.Default<String>("idref", "branchRates"), true);
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
            writer.writeTag(DiscretizedBranchRates.DISCRETIZED_BRANCH_RATES, new Attribute.Default<String>("idref", "branchRates"), true);
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
            writer.writeTag(DiscretizedBranchRates.DISCRETIZED_BRANCH_RATES, new Attribute.Default<String>("idref", "branchRates"), true);
            writer.writeCloseTag(RateCovarianceStatistic.RATE_COVARIANCE_STATISTIC);
        }
    }

}
