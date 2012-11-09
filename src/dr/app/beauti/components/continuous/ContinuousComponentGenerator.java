package dr.app.beauti.components.continuous;

import dr.app.beauti.generator.BaseComponentGenerator;
import dr.app.beauti.options.*;
import dr.app.beauti.util.XMLWriter;
import dr.evolution.datatype.ContinuousDataType;
import dr.evolution.util.Taxon;
import dr.evomodelxml.tree.TreeLoggerParser;
import dr.util.Attribute;
import dr.xml.AttributeParser;

/**
 * @author Andrew Rambaut
 * @version $Id$
 */

public class ContinuousComponentGenerator extends BaseComponentGenerator {

    protected ContinuousComponentGenerator(BeautiOptions options) {
        super(options);
    }

    public boolean usesInsertionPoint(InsertionPoint point) {

        if (options.getDataPartitions(ContinuousDataType.INSTANCE).size() == 0) {
            // Empty, so do nothing
            return false;
        }

        switch (point) {
//            case IN_TAXON:
            case AFTER_SITE_MODEL:
            case AFTER_TREE_LIKELIHOOD:
            case IN_OPERATORS:
            case IN_MCMC_PRIOR:
            case IN_MCMC_LIKELIHOOD:
            case IN_FILE_LOG_PARAMETERS:
            case IN_FILE_LOG_LIKELIHOODS:
            case IN_TREES_LOG:
                return true;
            default:
                return false;
        }

    }// END: usesInsertionPoint

    protected void generate(InsertionPoint point, Object item, XMLWriter writer) {

        ContinuousComponentOptions component = (ContinuousComponentOptions) options
                .getComponentOptions(ContinuousComponentOptions.class);

        switch (point) {
//          Don't need this because all traits are written for all taxa:
//            case IN_TAXON:
//                Taxon taxon = (Taxon)item;
//                writeTaxonTraits(taxon, writer);
//                break;
            case AFTER_SITE_MODEL:
                writeMultivariateDiffusionModels(writer, component);
                break;
            case AFTER_TREE_LIKELIHOOD:
                writeMultivariateTreeLikelihoods(writer, component);
                break;
            case IN_OPERATORS:
                writeRRWOperators(writer, component);
                writePrecisionGibbsOperators(writer, component);
                break;
            case IN_MCMC_PRIOR:
                writeMultivariatePriors(writer, component);
                break;
            case IN_MCMC_LIKELIHOOD:
                writeMultivariateTreeLikelihoodIdRefs(writer, component);
                break;
            case IN_FILE_LOG_PARAMETERS:
                writePrecisionMatrixIdRefs(writer, component);
                break;
            case IN_FILE_LOG_LIKELIHOODS:
                writeMultivariateTreeLikelihoodIdRefs(writer, component);
                break;
            case IN_TREES_LOG:
                writeTreeLogEntries((PartitionTreeModel)item, writer);
                break;
            default:
                throw new IllegalArgumentException(
                        "This insertion point is not implemented for "
                                + this.getClass().getName());
        }

    }// END: generate

    protected String getCommentLabel() {
        return "Multivariate diffusion model";
    }

    private void writeTaxonTraits(Taxon taxon, XMLWriter writer) {
        for (AbstractPartitionData partition : options.getDataPartitions(ContinuousDataType.INSTANCE)) {
            writer.writeOpenTag(AttributeParser.ATTRIBUTE, new Attribute[]{
                    new Attribute.Default<String>(Attribute.NAME, partition.getName())});

            StringBuilder sb = new StringBuilder();
            boolean first = true;
            for (TraitData trait : partition.getTraits()) {
                if (!first) {
                    sb.append(" ");
                }

                if (taxon.containsAttribute(trait.getName())) {
                    sb.append(taxon.getAttribute(trait.getName()).toString());
                } else {
                    sb.append("?");
                }
                first = false;
            }
            writer.writeText(sb.toString());
            writer.writeCloseTag(AttributeParser.ATTRIBUTE);
        }
    }

    private void writeMultivariateDiffusionModels(XMLWriter writer,
                                                  ContinuousComponentOptions component) {

        boolean first = false;

        for (PartitionSubstitutionModel model : component.getOptions().getPartitionSubstitutionModels(ContinuousDataType.INSTANCE)) {
            String precisionMatrixId = model.getName() + ".precision";

            if (!first) { writer.writeBlankLine(); } else {  first = false;  }

            writeMultivariateDiffusionModel(writer, model, precisionMatrixId);

            writer.writeBlankLine();

            writeMultivariateWishartPrior(writer, model, precisionMatrixId);
        }
    }

    private void writeMultivariateDiffusionModel(XMLWriter writer,
                                                 PartitionSubstitutionModel model,
                                                 String precisionMatrixId) {

        writer.writeOpenTag("multivariateDiffusionModel",
                new Attribute[] {
                        new Attribute.Default<String>("id", model.getName() + ".diffusionModel")
                });

        writer.writeOpenTag("precisionMatrix");
        writer.writeOpenTag("matrixParameter",
                new Attribute[] {
                        new Attribute.Default<String>("id", precisionMatrixId)
                });

        for (int i = 0; i < model.getContinuousTraitCount(); i++) {
            StringBuilder sb = new StringBuilder();
            for (int j = 0; j < model.getContinuousTraitCount(); j++) {
                if (j > 0) {
                    sb.append(" ");
                }
                if (i == j) {
                    sb.append(0.05);
                } else {
                    sb.append(0.002);
                }
            }
            writer.writeTag("parameter",
                    new Attribute[]{
                            new Attribute.Default<String>("id", "col" + (i + 1)),
                            new Attribute.Default<String>("value", sb.toString())
                    }, true);
        }

        writer.writeCloseTag("matrixParameter");
        writer.writeCloseTag("precisionMatrix");
        writer.writeCloseTag("multivariateDiffusionModel");
    }

    private void writeMultivariateWishartPrior(XMLWriter writer,
                                               PartitionSubstitutionModel model,
                                               String precisionMatrixId) {

        int n = model.getContinuousTraitCount();

        writer.writeOpenTag("multivariateWishartPrior",
                new Attribute[] {
                        new Attribute.Default<String>("id", model.getName() + ".precisionPrior"),
                        new Attribute.Default<String>("df", "" + n),
                });

        writer.writeOpenTag("scaleMatrix");
        writer.writeOpenTag("matrixParameter");

        for (int i = 0; i < n; i++) {
            StringBuilder sb = new StringBuilder();
            for (int j = 0; j < n; j++) {
                if (j > 0) {
                    sb.append(" ");
                }
                if (i == j) {
                    sb.append(1.0);
                } else {
                    sb.append(0.0);
                }
            }
            writer.writeTag("parameter",
                    new Attribute[]{
                            new Attribute.Default<String>("value", sb.toString())
                    }, true);

        }

        writer.writeCloseTag("matrixParameter");
        writer.writeCloseTag("scaleMatrix");

        writer.writeOpenTag("data");
        writer.writeIDref("parameter", precisionMatrixId);
        writer.writeCloseTag("data");

        writer.writeCloseTag("multivariateWishartPrior");
    }

    private void writeMultivariateTreeLikelihoods(XMLWriter writer,
                                                  ContinuousComponentOptions component) {

        boolean first = true;

        for (AbstractPartitionData partitionData : component.getOptions().getDataPartitions(ContinuousDataType.INSTANCE)) {
            PartitionSubstitutionModel model = partitionData.getPartitionSubstitutionModel();
            String diffusionModelId = model.getName() + ".diffusionModel";
            String treeModelId = partitionData.getPartitionTreeModel().getPrefix() + "treeModel";

            if (!first) { writer.writeBlankLine(); } else {  first = false;  }

            if (model.getContinuousSubstModelType() != ContinuousSubstModelType.HOMOGENOUS) {
                writeRelaxedBranchRateModel(writer, partitionData, treeModelId);
            }

            writer.writeBlankLine();

            writeMultivariateTreeLikelihood(writer, partitionData, diffusionModelId, treeModelId);

            if (partitionData.getTraits().size() == 2) {
                // if we are analysing bivariate traits we can add these special statistics...
                String precisionMatrixId = model.getName() + ".precision";
                write2DStatistics(writer, partitionData, precisionMatrixId, treeModelId);
            }
        }
    }


    private void writeRelaxedBranchRateModel(XMLWriter writer,
                                             AbstractPartitionData partitionData,
                                             String treeModelId) {

        String prefix = partitionData.getName() + ".";

        writer.writeOpenTag("discretizedBranchRates",
                new Attribute[] {
                        new Attribute.Default<String>("id",
                                prefix + "diffusionRates"),
                });

        writer.writeIDref("treeModel", treeModelId);

        writer.writeOpenTag("distribution");
        writer.writeOpenTag("onePGammaDistributionModel");
        // don't think this needs an id
//                new Attribute[] {
//                        new Attribute.Default<String>("id", prefix + "gamma"),
//                });

        writer.writeOpenTag("shape");
        switch (partitionData.getPartitionSubstitutionModel().getContinuousSubstModelType()) {

            case HOMOGENOUS:
                throw new IllegalArgumentException("Shouldn't be here");
            case CAUCHY_RRW:
                writer.writeComment("half DF (i.e., df = 1)");
                writer.writeTag("parameter",
                        new Attribute[]{
                                // don't think this needs an id
//                        new Attribute.Default<String>("id", "halfDF"),
                                new Attribute.Default<String>("value", "0.5")
                        }, true);
                break;
            case GAMMA_RRW:
                writer.writeComment("half DF");
                writer.writeTag("parameter",
                        new Attribute[]{
                                new Attribute.Default<String>("id", prefix + "halfDF"),
                                new Attribute.Default<String>("value", "0.5")
                        }, true);
                break;
            default:
                throw new IllegalArgumentException("Unknown continuous substitution type");
        }
        writer.writeCloseTag("shape");
        writer.writeCloseTag("onePGammaDistributionModel");
        writer.writeCloseTag("distribution");


        writer.writeOpenTag("rateCategories");
        writer.writeTag("parameter", new Attribute.Default<String>("id", prefix + "rrwCategories"), true);
        writer.writeCloseTag("rateCategories");

        writer.writeCloseTag("discretizedBranchRates");
    }

    private void writeMultivariateTreeLikelihood(XMLWriter writer,
                                                 AbstractPartitionData partitionData,
                                                 String diffusionModelId,
                                                 String treeModelId) {

        writer.writeOpenTag("multivariateTraitLikelihood",
                new Attribute[] {
                        new Attribute.Default<String>("id", partitionData.getName() + ".traitLikelihood"),
                        new Attribute.Default<String>("traitName", partitionData.getName()),
                        new Attribute.Default<String>("useTreeLength", "true"),
                        new Attribute.Default<String>("scaleByTime", "true"),
                        new Attribute.Default<String>("reportAsMultivariate", "true"),
                        new Attribute.Default<String>("reciprocalRates", "true"),
                        new Attribute.Default<String>("integrateInternalTraits", "true")
                });

        writer.writeIDref("multivariateDiffusionModel", diffusionModelId);
        writer.writeIDref("treeModel", treeModelId);

        writer.writeOpenTag("traitParameter");
        writer.writeTag("parameter", new Attribute.Default<String>("id", "leaf." + partitionData.getName()), true);
        writer.writeCloseTag("traitParameter");

        writer.writeOpenTag("conjugateRootPrior");

        writer.writeOpenTag("meanParameter");
        StringBuilder sb = new StringBuilder();
        for (int j = 0; j < partitionData.getTraits().size(); j++) {
            if (j > 0) {
                sb.append(" ");
            }
            sb.append(0.0);
        }

        writer.writeTag("parameter", new Attribute.Default<String>("value", sb.toString()), true);
        writer.writeCloseTag("meanParameter");

        writer.writeOpenTag("priorSampleSize");
        writer.writeTag("parameter", new Attribute.Default<String>("value", "0.001"), true);
        writer.writeCloseTag("priorSampleSize");

        writer.writeCloseTag("conjugateRootPrior");

        if (partitionData.getPartitionSubstitutionModel().getContinuousSubstModelType() != ContinuousSubstModelType.HOMOGENOUS) {
            writer.writeIDref("discretizedBranchRates", partitionData.getName() + "." + "diffusionRates");
        }

        writer.writeCloseTag("multivariateTraitLikelihood");
    }

    private void write2DStatistics(XMLWriter writer, AbstractPartitionData partitionData, String precisionMatrixId, String treeModelId) {
        String prefix = partitionData.getName() + ".";

        writer.writeOpenTag("correlation",
                new Attribute[] {
                        new Attribute.Default<String>("id", prefix + "correlation")
                });
        writer.writeIDref("matrixParameter", precisionMatrixId);
        writer.writeCloseTag("correlation");

        writer.writeOpenTag("treeLengthStatistic",
                new Attribute[] {
                        new Attribute.Default<String>("id", prefix + "treeLength")
                });
        writer.writeIDref("treeModel", treeModelId);
        writer.writeCloseTag("treeLengthStatistic");


        writer.writeOpenTag("productStatistic",
                new Attribute[] {
                        new Attribute.Default<String>("id", prefix + "treeLengthPrecision1")
                });
        writer.writeIDref("treeLengthStatistic", prefix + "treeLength");
        writer.writeOpenTag("subStatistic",
                new Attribute[] {
                        new Attribute.Default<String>("id", prefix + "precision1"),
                        new Attribute.Default<String>("dimension", "0")
                });
        writer.writeIDref("parameter", prefix + "col1");
        writer.writeCloseTag("subStatistic");
        writer.writeCloseTag("productStatistic");

        writer.writeOpenTag("productStatistic",
                new Attribute[] {
                        new Attribute.Default<String>("id", prefix + "treeLengthPrecision2")
                });
        writer.writeIDref("treeLengthStatistic", prefix + "treeLength");
        writer.writeOpenTag("subStatistic",
                new Attribute[] {
                        new Attribute.Default<String>("id", prefix + "precision2"),
                        new Attribute.Default<String>("dimension", "1")
                });
        writer.writeIDref("parameter", prefix + "col2");
        writer.writeCloseTag("subStatistic");
        writer.writeCloseTag("productStatistic");

    }

    private void write2DStatisticsIDrefs(XMLWriter writer, AbstractPartitionData partitionData) {
        String prefix = partitionData.getName() + ".";
        writer.writeIDref("correlation", prefix + "correlation");

        writer.writeIDref("treeLengthStatistic", prefix + "treeLength");

        writer.writeIDref("productStatistic", prefix + "treeLengthPrecision1");

        writer.writeIDref("productStatistic", prefix + "treeLengthPrecision2");
    }

    /*
 <correlation id="locationCorrelation" dimension1="1" dimension2="2">
 <matrixParameter idref="precisionMatrix"/>
 </correlation>

 <treeLengthStatistic id="treeLength">
 <treeModel idref="treeModel"/>
 </treeLengthStatistic>

 <productStatistic id="treeLengthPrecision1">
 <treeLengthStatistic idref="treeLength"/>
 <subStatistic id="precision1" dimension="0">  <!-- I do not  know why Joseph programmed these to start counting from 0 -->
 <parameter idref="col1"/>
 </subStatistic>
 </productStatistic>

 <productStatistic id="treeLengthPrecision2">
 <treeLengthStatistic idref="treeLength"/>
 <subStatistic id="precision2" dimension="1">
 <parameter idref="col2"/>
 </subStatistic>
 </productStatistic>

 <treeDispersionStatistic id="dispersionRate" greatCircleDistance="true">
 <treeModel idref="treeModel"/>
 <multivariateTraitLikelihood idref="traitLikelihood"/>
 </treeDispersionStatistic>
    */

    private void writeRRWOperators(XMLWriter writer,
                                   ContinuousComponentOptions component) {

        for (AbstractPartitionData partitionData : component.getOptions().getDataPartitions(ContinuousDataType.INSTANCE)) {
            ContinuousSubstModelType type = partitionData.getPartitionSubstitutionModel().getContinuousSubstModelType();

            if (type != ContinuousSubstModelType.HOMOGENOUS) {
                // this is now in the parameter table...
//                if (type == ContinuousSubstModelType.GAMMA_RRW) {
//                    writer.writeOpenTag("scaleOperator",
//                            new Attribute[] {
//                                    new Attribute.Default<String>("scaleFactor", "0.75"),
//                                    new Attribute.Default<String>("weight", "1")
//                            });
//                    writer.writeIDref("parameter", partitionData.getName() + ".halfDF");
//                    writer.writeCloseTag("scaleOperator");
//                }

                writer.writeOpenTag("swapOperator",
                        new Attribute[] {
                                new Attribute.Default<String>("size", "1"),
                                new Attribute.Default<String>("weight", "30"),
                                new Attribute.Default<String>("autoOptimize", "false")
                        });
                writer.writeIDref("parameter", partitionData.getName() + ".rrwCategories");
                writer.writeCloseTag("swapOperator");

                // See Issue 500:
                // http://code.google.com/p/beast-mcmc/issues/detail?id=500&can=1&start=400
//                writer.writeOpenTag("randomWalkIntegerOperator",
//                        new Attribute[] {
//                                new Attribute.Default<String>("windowSize", "2"),
//                                new Attribute.Default<String>("weight", "10")
//                        });
//                writer.writeIDref("parameter", partitionData.getName() + ".rrwCategories");
//                writer.writeCloseTag("randomWalkIntegerOperator");

                writer.writeOpenTag("uniformIntegerOperator",
                        new Attribute[] {
                                new Attribute.Default<String>("weight", "10")
                        });
                writer.writeIDref("parameter", partitionData.getName() + ".rrwCategories");
                writer.writeCloseTag("uniformIntegerOperator");
            }
        }
    }

    private void writePrecisionGibbsOperators(XMLWriter writer,
                                              ContinuousComponentOptions component) {

        for (AbstractPartitionData partitionData : component.getOptions().getDataPartitions(ContinuousDataType.INSTANCE)) {
            writePrecisionGibbsOperator(writer, component, partitionData);
        }
    }

    private void writePrecisionGibbsOperator(final XMLWriter writer,
                                             final ContinuousComponentOptions component,
                                             AbstractPartitionData partitionData
    ) {
        writer.writeOpenTag("precisionGibbsOperator",
                new Attribute[] {
                        new Attribute.Default<String>("weight", "" + partitionData.getTraits().size())
                });
        writer.writeIDref("multivariateTraitLikelihood", partitionData.getName() + ".traitLikelihood");
        writer.writeIDref("multivariateWishartPrior", partitionData.getPartitionSubstitutionModel().getName() + ".precisionPrior");
        writer.writeCloseTag("precisionGibbsOperator");

    }

    private void writePrecisionMatrixIdRefs(final XMLWriter writer, final ContinuousComponentOptions component) {
        for (PartitionSubstitutionModel model : component.getOptions().getPartitionSubstitutionModels(ContinuousDataType.INSTANCE)) {
            writer.writeIDref("matrixParameter", model.getName() + ".precision");

            if (model.getContinuousTraitCount() == 2) {
                // if we are analysing bivariate traits we can add these special statistics...
                write2DStatisticsIDrefs(writer, null);
            }
        }
    }

    private void writeMultivariatePriors(XMLWriter writer,
                                                       ContinuousComponentOptions component) {

        for (AbstractPartitionData partitionData : component.getOptions().getDataPartitions(ContinuousDataType.INSTANCE)) {
            writer.writeIDref("multivariateWishartPrior", partitionData.getName() + ".precisionPrior");
        }
    }

    private void writeMultivariateTreeLikelihoodIdRefs(XMLWriter writer,
                                                       ContinuousComponentOptions component) {

        for (AbstractPartitionData partitionData : component.getOptions().getDataPartitions(ContinuousDataType.INSTANCE)) {
            writer.writeIDref("multivariateTraitLikelihood", partitionData.getName() + ".traitLikelihood");
        }
    }

    private void writeTreeLogEntries(PartitionTreeModel treeModel, XMLWriter writer) {
        for (AbstractPartitionData partitionData : options.getDataPartitions(ContinuousDataType.INSTANCE)) {
            if (partitionData.getPartitionTreeModel() == treeModel) {
                PartitionSubstitutionModel model = partitionData.getPartitionSubstitutionModel();
                writer.writeIDref("multivariateDiffusionModel", model.getName() + ".diffusionModel");
                writer.writeIDref("multivariateTraitLikelihood", partitionData.getName() + ".traitLikelihood");
                if (model.getContinuousSubstModelType() != ContinuousSubstModelType.HOMOGENOUS) {
                    writer.writeOpenTag(TreeLoggerParser.TREE_TRAIT,
                            new Attribute[] {
                                    new Attribute.Default<String>(TreeLoggerParser.NAME, "rate"),
                                    new Attribute.Default<String>(TreeLoggerParser.TAG, partitionData.getName() + ".rate"),
                            });
                    writer.writeIDref("discretizedBranchRates",  partitionData.getName() + "." + "diffusionRates");
                    writer.writeCloseTag(TreeLoggerParser.TREE_TRAIT);
                }
            }
        }
    }
}
