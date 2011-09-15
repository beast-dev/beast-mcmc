package dr.app.beauti.components.continuous;

import dr.app.beauti.generator.BaseComponentGenerator;
import dr.app.beauti.options.*;
import dr.app.beauti.util.XMLWriter;
import dr.app.treespace.InputFile;
import dr.evolution.continuous.Continuous;
import dr.evolution.datatype.ContinuousDataType;
import dr.evolution.datatype.GeneralDataType;
import dr.evolution.util.Taxon;
import dr.evomodelxml.treelikelihood.AncestralStateTreeLikelihoodParser;
import dr.evomodelxml.treelikelihood.TreeLikelihoodParser;
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

        if (options.getAllPartitionData(ContinuousDataType.INSTANCE).size() == 0) {
            // Empty, so do nothing
            return false;
        }

        switch (point) {
            case IN_TAXON:
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
            case IN_TAXON:
                Taxon taxon = (Taxon)item;
                writeTaxonTraits(taxon, writer);
                break;
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
        for (AbstractPartitionData partition : options.getAllPartitionData(ContinuousDataType.INSTANCE)) {
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

        for (PartitionSubstitutionModel model : component.getOptions().getPartitionSubstitutionModels(ContinuousDataType.INSTANCE)) {
            String precisionMatrixId = model.getPrefix(ContinuousDataType.INSTANCE) + "precision";
            writeMultivariateDiffusionModel(writer, model, precisionMatrixId);
            writeMultivariateWishartPrior(writer, model, precisionMatrixId);
        }
    }

    private void writeMultivariateDiffusionModel(XMLWriter writer,
                                                 PartitionSubstitutionModel model,
                                                 String precisionMatrixId) {

        writer.writeOpenTag("multivariateDiffusionModel",
                new Attribute[] {
                        new Attribute.Default<String>("id", model.getPrefix(ContinuousDataType.INSTANCE) + "diffusionModel")
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

        writer.writeBlankLine();
    }

    private void writeMultivariateWishartPrior(XMLWriter writer,
                                               PartitionSubstitutionModel model,
                                               String precisionMatrixId) {

        int n = model.getContinuousTraitCount();

        writer.writeOpenTag("multivariateWishartPrior",
                new Attribute[] {
                        new Attribute.Default<String>("id", model.getPrefix(ContinuousDataType.INSTANCE) + "precisionPrior"),
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

        for (AbstractPartitionData partitionData : component.getOptions().getAllPartitionData(ContinuousDataType.INSTANCE)) {
            PartitionSubstitutionModel model = partitionData.getPartitionSubstitutionModel();
            String diffusionModelId = model.getPrefix(ContinuousDataType.INSTANCE) + "diffusionModel";
            String treeModelId = partitionData.getPartitionTreeModel().getPrefix() + "treeModel";

            if (model.getContinuousSubstModelType() != ContinuousSubstModelType.HOMOGENOUS) {
                writeRelaxedBranchRateModel(writer, partitionData, treeModelId);
            }

            writeMultivariateTreeLikelihood(writer, partitionData, diffusionModelId, treeModelId);
        }
    }

    private void writeRelaxedBranchRateModel(XMLWriter writer,
                                             AbstractPartitionData partitionData,
                                             String treeModelId) {

        String prefix = partitionData.getPrefix(ContinuousDataType.INSTANCE);

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

//        <discretizedBranchRates id="branchRatesDiffusion">
//            <treeModel idref="treeModel"/>
//            <distribution>
//                <onePGammaDistributionModel id="gammaModel">
//                     <shape>
//                        <parameter id="halfDF" value="0.5"/>
//                    </shape>
//                </onePGammaDistributionModel>
//            </distribution>
//            <rateCategories>
//                <parameter id="branchRatesDiffusion.categories" dimension="92"/>
//            </rateCategories>
//        </discretizedBranchRates>
//
//        <rateStatistic id="meanRateDiffusion" name="meanRate" mode="mean" internal="true" external="true">
//            <treeModel idref="treeModel"/>
//            <discretizedBranchRates idref="branchRatesDiffusion"/>
//        </rateStatistic>
//
//        <rateStatistic id="coefficientOfVariationDiffusion" name="coefficientOfVariationDiffusion" mode="coefficientOfVariation" internal="true" external="true">
//            <treeModel idref="treeModel"/>
//            <discretizedBranchRates idref="branchRatesDiffusion"/>
//        </rateStatistic>
//
//        <rateCovarianceStatistic id="covarianceDiffusion" name="covarianceDiffusion">
//            <treeModel idref="treeModel"/>
//            <discretizedBranchRates idref="branchRatesDiffusion"/>
//        </rateCovarianceStatistic>
    }

    private void writeMultivariateTreeLikelihood(XMLWriter writer,
                                                 AbstractPartitionData partitionData,
                                                 String diffusionModelId,
                                                 String treeModelId) {

        writer.writeOpenTag("multivariateTraitLikelihood",
                new Attribute[] {
                        new Attribute.Default<String>("id", partitionData.getPrefix(ContinuousDataType.INSTANCE) + "traitLikelihood"),
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

        writer.writeCloseTag("multivariateTraitLikelihood");
    }

    private void writeRRWOperators(XMLWriter writer,
                                   ContinuousComponentOptions component) {

        for (AbstractPartitionData partitionData : component.getOptions().getAllPartitionData(ContinuousDataType.INSTANCE)) {
            ContinuousSubstModelType type = partitionData.getPartitionSubstitutionModel().getContinuousSubstModelType();

            if (type != ContinuousSubstModelType.HOMOGENOUS) {
                // this is now in the parameter table...
//                if (type == ContinuousSubstModelType.GAMMA_RRW) {
//                    writer.writeOpenTag("scaleOperator",
//                            new Attribute[] {
//                                    new Attribute.Default<String>("scaleFactor", "0.75"),
//                                    new Attribute.Default<String>("weight", "1")
//                            });
//                    writer.writeIDref("parameter", partitionData.getPrefix(ContinuousDataType.INSTANCE) + "halfDF");
//                    writer.writeCloseTag("scaleOperator");
//                }

                writer.writeOpenTag("swapOperator",
                        new Attribute[] {
                                new Attribute.Default<String>("size", "1"),
                                new Attribute.Default<String>("weight", "30"),
                                new Attribute.Default<String>("autoOptimize", "false")
                        });
                writer.writeIDref("parameter", partitionData.getPrefix(ContinuousDataType.INSTANCE) + "rrwCategories");
                writer.writeCloseTag("swapOperator");

                writer.writeOpenTag("randomWalkIntegerOperator",
                        new Attribute[] {
                                new Attribute.Default<String>("size", "2"),
                                new Attribute.Default<String>("weight", "10")
                        });
                writer.writeIDref("parameter", partitionData.getPrefix(ContinuousDataType.INSTANCE) + "rrwCategories");
                writer.writeCloseTag("randomWalkIntegerOperator");

                writer.writeOpenTag("uniformIntegerOperator",
                        new Attribute[] {
                                new Attribute.Default<String>("weight", "10")
                        });
                writer.writeIDref("parameter", partitionData.getPrefix(ContinuousDataType.INSTANCE) + "rrwCategories");
                writer.writeCloseTag("uniformIntegerOperator");
            }
        }
    }

    private void writePrecisionGibbsOperators(XMLWriter writer,
                                              ContinuousComponentOptions component) {

        for (AbstractPartitionData partitionData : component.getOptions().getAllPartitionData(ContinuousDataType.INSTANCE)) {
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
        writer.writeIDref("multivariateTraitLikelihood", partitionData.getPrefix(ContinuousDataType.INSTANCE) + "traitLikelihood");
        writer.writeIDref("multivariateWishartPrior", partitionData.getPartitionSubstitutionModel().getPrefix(ContinuousDataType.INSTANCE) + "precisionPrior");
        writer.writeCloseTag("precisionGibbsOperator");

    }

    private void writePrecisionMatrixIdRefs(final XMLWriter writer, final ContinuousComponentOptions component) {
        for (PartitionSubstitutionModel model : component.getOptions().getPartitionSubstitutionModels(ContinuousDataType.INSTANCE)) {
            writer.writeIDref("matrixParameter", model.getPrefix(ContinuousDataType.INSTANCE) + "precision");
        }
    }

    private void writeMultivariatePriors(XMLWriter writer,
                                                       ContinuousComponentOptions component) {

        for (AbstractPartitionData partitionData : component.getOptions().getAllPartitionData(ContinuousDataType.INSTANCE)) {
            writer.writeIDref("multivariateWishartPrior", partitionData.getPrefix(ContinuousDataType.INSTANCE) + "precisionPrior");
        }
    }

    private void writeMultivariateTreeLikelihoodIdRefs(XMLWriter writer,
                                                       ContinuousComponentOptions component) {

        for (AbstractPartitionData partitionData : component.getOptions().getAllPartitionData(ContinuousDataType.INSTANCE)) {
            writer.writeIDref("multivariateTraitLikelihood", partitionData.getPrefix(ContinuousDataType.INSTANCE) + "traitLikelihood");
        }
    }

    private void writeTreeLogEntries(PartitionTreeModel treeModel, XMLWriter writer) {
        for (AbstractPartitionData partitionData : options.getAllPartitionData(ContinuousDataType.INSTANCE)) {
            if (partitionData.getPartitionTreeModel() == treeModel) {
                PartitionSubstitutionModel model = partitionData.getPartitionSubstitutionModel();
                writer.writeIDref("multivariateDiffusionModel", model.getPrefix(ContinuousDataType.INSTANCE) + "diffusionModel");
                writer.writeIDref("multivariateTraitLikelihood", partitionData.getPrefix(ContinuousDataType.INSTANCE) + "traitLikelihood");
            }
        }
    }
}
