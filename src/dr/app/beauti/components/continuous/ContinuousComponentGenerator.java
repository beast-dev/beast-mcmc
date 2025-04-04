/*
 * ContinuousComponentGenerator.java
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

package dr.app.beauti.components.continuous;

import dr.app.beauti.components.GeneratorHelper;
import dr.app.beauti.components.XMLWriterObject;
import dr.app.beauti.generator.BaseComponentGenerator;
import dr.app.beauti.options.*;
import dr.app.beauti.util.XMLWriter;
import dr.evolution.datatype.ContinuousDataType;
import dr.evolution.util.Taxon;
import dr.evomodel.continuous.TreeDataContinuousDiffusionStatistic;
import dr.evomodel.treedatalikelihood.continuous.IntegratedFactorAnalysisLikelihoodParser;
import dr.evomodelxml.treedatalikelihood.continuous.RepeatedMeasuresTraitDataModelParser;
import dr.evomodel.treedatalikelihood.continuous.RepeatedMeasuresWishartStatistics;
import dr.evomodel.treedatalikelihood.continuous.WishartStatisticsWrapper;
import dr.evomodelxml.tree.TreeLoggerParser;
import dr.evomodelxml.treedatalikelihood.ContinuousDataLikelihoodParser;
import dr.inference.model.ParameterParser;
import dr.util.Attribute;
import dr.xml.AttributeParser;

/**
 * @author Andrew Rambaut
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

    protected void generate(final InsertionPoint point, final Object item, final String prefix, final XMLWriter writer) {

        ContinuousComponentOptions component = (ContinuousComponentOptions) options
                .getComponentOptions(ContinuousComponentOptions.class);

        switch (point) {
            case IN_TAXON:
                Taxon taxon = (Taxon) item;
                writeTaxonTraits(taxon, writer);
                break;
            case AFTER_SITE_MODEL:
                writeMultivariateDiffusionModels(writer, component);
                writeContinuousExtensionModels(writer, component);
                break;
            case AFTER_TREE_LIKELIHOOD:
                writeMultivariateTreeLikelihoods(writer, component);
                break;
            case IN_OPERATORS:
                writePrecisionGibbsOperators(writer, component);
                break;
            case IN_MCMC_PRIOR:
                writeMultivariatePriors(writer, component);
                break;
            case IN_MCMC_LIKELIHOOD:
                writeMultivariateTreeLikelihoodIdRefs(writer, component);
                break;
            case IN_FILE_LOG_PARAMETERS:
                writeParameterIdRefs(writer, component);
                break;
            case IN_FILE_LOG_LIKELIHOODS:
                writeMultivariateTreeLikelihoodIdRefs(writer, component);
                break;
            case IN_TREES_LOG:
                writeTreeLogEntries((PartitionTreeModel) item, writer);
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

            if (!first) {
                writer.writeBlankLine();
            } else {
                first = false;
            }

            writeMultivariateDiffusionModel(writer, model, precisionMatrixId);

            writer.writeBlankLine();

            String wishartId = model.getName() + ".precisionPrior";

            writeMultivariateWishartPrior(writer, wishartId, precisionMatrixId, model.getContinuousTraitDimension());
        }
    }

    private void writeMultivariateDiffusionModel(XMLWriter writer,
                                                 PartitionSubstitutionModel model,
                                                 String precisionMatrixId) {

        writer.writeOpenTag("multivariateDiffusionModel",
                new Attribute[]{
                        new Attribute.Default<String>("id", model.getName() + ".diffusionModel")
                });

        writer.writeOpenTag("precisionMatrix");
        writer.writeOpenTag("matrixParameter",
                new Attribute[]{
                        new Attribute.Default<String>("id", precisionMatrixId)
                });

        double diagValue = (model.getContinuousExtensionType() == ContinuousModelExtensionType.LATENT_FACTORS) ? 1.0 : 0.05;
        double offDiagValue = (model.getContinuousExtensionType() == ContinuousModelExtensionType.LATENT_FACTORS) ? 0.0 : 0.002;

        for (int i = 0; i < model.getContinuousTraitDimension(); i++) {
            StringBuilder sb = new StringBuilder();
            for (int j = 0; j < model.getContinuousTraitDimension(); j++) {
                if (j > 0) {
                    sb.append(" ");
                }
                if (i == j) {
                    sb.append(diagValue);
                } else {
                    sb.append(offDiagValue);
                }
            }
            writer.writeTag("parameter",
                    new Attribute[]{
                            new Attribute.Default<String>("id", precisionMatrixId + ".col" + (i + 1)),
                            new Attribute.Default<String>("value", sb.toString())
                    }, true);
        }

        writer.writeCloseTag("matrixParameter");
        writer.writeCloseTag("precisionMatrix");
        writer.writeCloseTag("multivariateDiffusionModel");
    }


    private void writeNormalDistributionPrior(XMLWriter writer, String id, String parameterId) {
//        String distributionLikelihood = "distributionLikelihood";
//        writer.writeOpenTag(distributionLikelihood);
//
//        String data = "data";
//        writer.writeOpenTag(data);
//        writer.writeIDref("parameter", parameterId);
//        writer.writeCloseTag(data);
//        writer.writeCloseTag(distributionLikelihood);
        XMLWriterObject param = new XMLWriterObject("parameter", parameterId);
        param.setAlreadyWritten(true);
        XMLWriterObject data = new XMLWriterObject("data", param);
        XMLWriterObject distribution = getNormalDistribution(0, 1);

        XMLWriterObject prior = new XMLWriterObject(
                "distributionLikelihood",
                id,
                new XMLWriterObject[]{
                        data,
                        new XMLWriterObject(
                                "distribution",
                                distribution
                        )
                });
        prior.writeOrReference(writer);
    }

    private XMLWriterObject getNormalDistribution(double mean, double stdev) {
//        writer.writeOpenTag("distribution");
//        writer.writeOpenTag("mean");
//        writer.writeCloseTag("mean");
//        writer.writeOpenTag("stdev");
//        writer.writeCloseTag("stdev");
//        writer.writeCloseTag("");
        XMLWriterObject distribution = new XMLWriterObject(
                "normalDistributionModel",
                null,
                new XMLWriterObject[]{
                        new XMLWriterObject("mean",
                                new XMLWriterObject("parameter",
                                        null,
                                        null,
                                        new Attribute[]{
                                                new Attribute.Default("value", mean)
                                        })),
                        new XMLWriterObject("stdev",
                                new XMLWriterObject("parameter",
                                        null,
                                        null,
                                        new Attribute[]{
                                                new Attribute.Default("value", stdev),
                                                new Attribute.Default("lower", 0)
                                        }))
                }
        );

        return distribution;
    }

    private void writeMultivariateWishartPrior(XMLWriter writer,
                                               String id,
                                               String precisionMatrixId,
                                               int n) {


        writer.writeOpenTag("multivariateWishartPrior",
                new Attribute[]{
                        new Attribute.Default<String>("id", id),
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

    private void writeContinuousExtensionModels(XMLWriter writer, ContinuousComponentOptions component) {

        boolean first = true;

        for (AbstractPartitionData partitionData :
                component.getOptions().getDataPartitions(ContinuousDataType.INSTANCE)) {
            PartitionSubstitutionModel model = partitionData.getPartitionSubstitutionModel();

            ContinuousModelExtensionType extensionType = model.getContinuousExtensionType();

            if (extensionType != ContinuousModelExtensionType.NONE) {


                String treeModelId = partitionData.getPartitionTreeModel().getPrefix() + "treeModel";

                if (first) {
                    writer.writeBlankLine();
                } else {
                    first = false;
                }

                String precisionId;

                switch (extensionType) {
                    case RESIDUAL:
                        precisionId = repeatedMeasuresTraitDataModelParser.getBeautiParameterIDProvider(
                                "extensionPrecision").getId(model.getName());
                        writeResidualExtensionModel(writer, partitionData, treeModelId, precisionId);

                        writer.writeBlankLine();

                        String wishartId = repeatedMeasuresTraitDataModelParser.getBeautiParameterIDProvider("extensionPrecision").getPriorId(model.getName());
                        writeMultivariateWishartPrior(writer, wishartId, precisionId, model.getExtendedTraitCount());


                        break;
                    case LATENT_FACTORS:
                        precisionId = integratedFactorsParser.getBeautiParameterIDProvider(
                                integratedFactorsParser.PRECISION).getId(model.getName());
                        String precisionPriorId = integratedFactorsParser.getBeautiParameterIDProvider(
                                integratedFactorsParser.PRECISION).getPriorId(model.getName());

                        String loadingsId = integratedFactorsParser.getBeautiParameterIDProvider(
                                integratedFactorsParser.LOADINGS).getId(model.getName());
                        String loadingsPriorId = integratedFactorsParser.getBeautiParameterIDProvider(
                                integratedFactorsParser.LOADINGS).getPriorId(model.getName());
                        writeLatentFactorModel(writer, partitionData, treeModelId, precisionId, loadingsId);
                        writer.writeBlankLine();

                        writeNormalDistributionPrior(writer, loadingsPriorId, loadingsId);

                        XMLWriterObject precisionPrior = new XMLWriterObject(
                                "gammaPrior",
                                precisionPriorId,
                                new XMLWriterObject[]{
                                        new XMLWriterObject("parameter",
                                                null,
                                                null,
                                                new Attribute[]{
                                                        new Attribute.Default("idref", precisionId)
                                                })
                                },
                                new Attribute[]{
                                        new Attribute.Default("scale", 1),
                                        new Attribute.Default("shape", 1)
                                }
                        );

                        precisionPrior.writeOrReference(writer);
                        break;
                    case NONE:
                        throw new IllegalArgumentException("Shouldn't be here");
                    default:
                        throw new IllegalArgumentException("Unknown continuous model extension type");
                }
            }


        }

    }

    private void writeResidualExtensionModel(XMLWriter writer, AbstractPartitionData partitionData, String treeModelId,
                                             String precisionId) {
        PartitionSubstitutionModel model = partitionData.getPartitionSubstitutionModel();
        int p = model.getExtendedTraitCount();


        writer.writeOpenTag("repeatedMeasuresModel",
                new Attribute[]{
                        new Attribute.Default<String>("id", repeatedMeasuresTraitDataModelParser.getId(
                                model.getName())),
                        new Attribute.Default<String>("traitName", partitionData.getName())
                });

        writer.writeIDref("treeModel", treeModelId);

        writeTraitParameter(writer, partitionData);

        writer.writeOpenTag("samplingPrecision");

        GeneratorHelper.writeIdentityMatrixParameter(writer, precisionId, p);


        writer.writeCloseTag("samplingPrecision");

        writer.writeCloseTag("repeatedMeasuresModel");

    }

    private void writeLatentFactorModel(XMLWriter writer, AbstractPartitionData partitionData, String treeModelId,
                                        String precisionId, String loadingsId) {


        PartitionSubstitutionModel model = partitionData.getPartitionSubstitutionModel();
        int p = model.getExtendedTraitCount();


        writer.writeOpenTag(integratedFactorsParser.getParserTag(),
                new Attribute[]{
                        new Attribute.Default<String>("id", integratedFactorsParser.getId(model.getName())),
                        new Attribute.Default<String>("traitName", partitionData.getName())
                });

        writer.writeIDref("treeModel", treeModelId);

        writeTraitParameter(writer, partitionData);

        writer.writeOpenTag(integratedFactorsParser.PRECISION);
        GeneratorHelper.writeParameter(writer, precisionId, p, 1);
        writer.writeCloseTag(integratedFactorsParser.PRECISION);

        writer.writeOpenTag(integratedFactorsParser.LOADINGS);
        GeneratorHelper.writeMatrixParameter(writer, loadingsId, model.getContinuousTraitDimension(), p);
        writer.writeCloseTag(integratedFactorsParser.LOADINGS);

        writer.writeCloseTag(integratedFactorsParser.getParserTag());
    }


    private void writeMultivariateTreeLikelihoods(XMLWriter writer,
                                                  ContinuousComponentOptions component) {

        boolean first = true;

        for (AbstractPartitionData partitionData : component.getOptions().getDataPartitions(ContinuousDataType.INSTANCE)) {
            PartitionSubstitutionModel model = partitionData.getPartitionSubstitutionModel();
            String diffusionModelId = model.getName() + ".diffusionModel";
            String treeModelId = partitionData.getPartitionTreeModel().getPrefix() + "treeModel";

            if (!first) {
                writer.writeBlankLine();
            } else {
                first = false;
            }

            if (model.getContinuousSubstModelType() != ContinuousSubstModelType.HOMOGENOUS &&
                    model.getContinuousSubstModelType() != ContinuousSubstModelType.DRIFT) {

                if (ContinuousComponentOptions.USE_ARBITRARY_BRANCH_RATE_MODEL) {
                    writeArbitraryRateBranchRateModel(writer, partitionData, treeModelId);
                } else {
                    writeRelaxedBranchRateModel(writer, partitionData, treeModelId);
                }
            }

            writer.writeBlankLine();

            writeMultivariateTreeLikelihood(writer, partitionData, diffusionModelId, treeModelId);

            String precisionMatrixId = model.getName() + ".precision";

            writeDiffusionStatistics(writer, partitionData, treeModelId, precisionMatrixId,
                    continuousDataLikelihoodParser.getId(partitionData.getName()));
        }
    }


    private void writeRelaxedBranchRateModel(XMLWriter writer,
                                             AbstractPartitionData partitionData,
                                             String treeModelId) {

        String prefix = partitionData.getName() + ".";

        writer.writeOpenTag("discretizedBranchRates",
                new Attribute[]{
                        new Attribute.Default<String>("id",
                                prefix + "diffusion.branchRates"),
                });

        writer.writeIDref("treeModel", treeModelId);

        writer.writeOpenTag("distribution");

        if (partitionData.getPartitionSubstitutionModel().getContinuousSubstModelType() == ContinuousSubstModelType.LOGNORMAL_RRW) {
            writer.writeOpenTag("logNormalDistributionModel",
                    new Attribute[]{new Attribute.Default<String>("meanInRealSpace", "true")});

            writer.writeOpenTag("mean");
            writer.writeTag("parameter",
                    new Attribute[]{
                            new Attribute.Default<String>("value", "1.0")
                    }, true);

            writer.writeCloseTag("mean");

            writer.writeOpenTag("stdev");
            writer.writeTag("parameter",
                    new Attribute[]{
                            new Attribute.Default<String>("id", prefix + ContinuousComponentOptions.STDEV),
                            new Attribute.Default<String>("value", "1.0"),
                            new Attribute.Default<String>("lower", "0.0")
                    }, true);

            writer.writeCloseTag("stdev");
            writer.writeCloseTag("logNormalDistributionModel");
        } else {
            writer.writeOpenTag("onePGammaDistributionModel");
            writer.writeOpenTag("shape");
            switch (partitionData.getPartitionSubstitutionModel().getContinuousSubstModelType()) {
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
                                    new Attribute.Default<String>("id", prefix + ContinuousComponentOptions.HALF_DF),
                                    new Attribute.Default<String>("value", "0.5")
                            }, true);
                    break;
                case LOGNORMAL_RRW:
                case DRIFT:
                case HOMOGENOUS:
                    throw new IllegalArgumentException("Shouldn't be here");
                default:
                    throw new IllegalArgumentException("Unknown continuous substitution type");
            }
            writer.writeCloseTag("shape");
            writer.writeCloseTag("onePGammaDistributionModel");
        }
        writer.writeCloseTag("distribution");


        writer.writeOpenTag("rateCategories");
        writer.writeTag("parameter", new Attribute.Default<String>("id", prefix + "rrwCategories"), true);
        writer.writeCloseTag("rateCategories");

        writer.writeCloseTag("discretizedBranchRates");
    }

    private void writeArbitraryRateBranchRateModel(XMLWriter writer,
                                                   AbstractPartitionData partitionData,
                                                   String treeModelId) {

        String prefix = partitionData.getName() + ".";

        writer.writeOpenTag("arbitraryBranchRates",
                new Attribute[]{
                        new Attribute.Default<String>("id",
                                prefix + "diffusion.branchRates"),
                });

        writer.writeIDref("treeModel", treeModelId);
        writer.writeOpenTag("rates");
        writer.writeTag("parameter",
                new Attribute[]{
                        new Attribute.Default<String>("id", prefix + "diffusion.rates"),
                        new Attribute.Default<String>("lower", "0.0")
                }, true);
        writer.writeCloseTag("rates");

        writer.writeCloseTag("arbitraryBranchRates");

        writer.writeOpenTag("distributionLikelihood",
                new Attribute[]{
                        new Attribute.Default<String>("id", prefix + "diffusion.prior")
                });

        writer.writeOpenTag("data");
        writer.writeIDref("parameter", prefix + "diffusion.rates");
        writer.writeCloseTag("data");

        writer.writeOpenTag("distribution");

        if (partitionData.getPartitionSubstitutionModel().getContinuousSubstModelType() == ContinuousSubstModelType.LOGNORMAL_RRW) {
            writer.writeOpenTag("logNormalDistributionModel",
                    new Attribute[]{new Attribute.Default<String>("meanInRealSpace", "true")});

            writer.writeOpenTag("mean");
            writer.writeTag("parameter",
                    new Attribute[]{
                            new Attribute.Default<String>("value", "1.0")
                    }, true);

            writer.writeCloseTag("mean");

            writer.writeOpenTag("stdev");
            writer.writeTag("parameter",
                    new Attribute[]{
                            new Attribute.Default<String>("id", prefix + ContinuousComponentOptions.STDEV),
                            new Attribute.Default<String>("value", "1.0"),
                            new Attribute.Default<String>("lower", "0.0")
                    }, true);

            writer.writeCloseTag("stdev");
            writer.writeCloseTag("logNormalDistributionModel");
        } else {
            writer.writeOpenTag("onePGammaDistributionModel");
            writer.writeOpenTag("shape");
            switch (partitionData.getPartitionSubstitutionModel().getContinuousSubstModelType()) {
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
                                    new Attribute.Default<String>("id", prefix + ContinuousComponentOptions.HALF_DF),
                                    new Attribute.Default<String>("value", "0.5")
                            }, true);
                    break;
                case LOGNORMAL_RRW:
                case DRIFT:
                case HOMOGENOUS:
                    throw new IllegalArgumentException("Shouldn't be here");
                default:
                    throw new IllegalArgumentException("Unknown continuous substitution type");
            }
            writer.writeCloseTag("shape");
            writer.writeCloseTag("onePGammaDistributionModel");
        }
        writer.writeCloseTag("distribution");

        writer.writeCloseTag("distributionLikelihood");
    }

    private void writeTraitParameter(XMLWriter writer, AbstractPartitionData partitionData) {
        writer.writeOpenTag("traitParameter");
        writer.writeTag("parameter", new Attribute.Default<String>("id", "leaf." + partitionData.getName()), true);
        writer.writeCloseTag("traitParameter");
    }

    private void writeMultivariateTreeLikelihood(XMLWriter writer,
                                                 AbstractPartitionData partitionData,
                                                 String diffusionModelId,
                                                 String treeModelId) {

        int traitDimension = partitionData.getPartitionSubstitutionModel().getContinuousTraitDimension(); //TODO: update


        writer.writeOpenTag(continuousDataLikelihoodParser.getParserTag(),
                new Attribute[]{
                        new Attribute.Default<String>("id", continuousDataLikelihoodParser.getId(partitionData.getName())),
                        new Attribute.Default<String>("traitName", partitionData.getName()),
                        new Attribute.Default<String>("useTreeLength", "true"),
                        new Attribute.Default<String>("scaleByTime", "true"),
                        new Attribute.Default<String>("reportAsMultivariate", "true"),
                        new Attribute.Default<String>("reciprocalRates", "false"),
                        new Attribute.Default<String>("integrateInternalTraits", "true")
                });

        writer.writeIDref("multivariateDiffusionModel", diffusionModelId);


        ContinuousComponentOptions component = (ContinuousComponentOptions) options
                .getComponentOptions(ContinuousComponentOptions.class);

        PartitionSubstitutionModel model = partitionData.getPartitionSubstitutionModel();

        if (component.useLambda(model)) {
            writer.writeOpenTag("transformedTreeModel");
            writer.writeIDref("treeModel", treeModelId);
            writer.writeTag("parameter", new Attribute[]{
                    new Attribute.Default<String>("id", partitionData.getName() + "." + ContinuousComponentOptions.LAMBDA),
                    new Attribute.Default<String>("value", "0.5"),
                    new Attribute.Default<String>("lower", "0.0"),
                    new Attribute.Default<String>("upper", "1.0")
            }, true);
            writer.writeCloseTag("transformedTreeModel");
        } else {
            writer.writeIDref("treeModel", treeModelId);
        }

        switch (model.getContinuousExtensionType()) {
            case NONE:
                writeTraitParameter(writer, partitionData);
                break;
            case RESIDUAL:
                writer.writeIDref("repeatedMeasuresModel", repeatedMeasuresTraitDataModelParser.getId(model.getName()));
                break;
            case LATENT_FACTORS:
                writer.writeIDref("integratedFactorModel", integratedFactorsParser.getId(model.getName()));
                break;
            default:
                throw new IllegalArgumentException("Unknown model extension type");
        }


        if (partitionData.getPartitionSubstitutionModel().getContinuousSubstModelType() == ContinuousSubstModelType.DRIFT) {
            writer.writeOpenTag("driftModels");
            for (int i = 0; i < traitDimension; i++) {
                writer.writeOpenTag("strictClockBranchRates");
                writer.writeOpenTag("rate");
                writer.writeTag("parameter", new Attribute[]{
                        new Attribute.Default<String>("id", partitionData.getName() + "." + ContinuousComponentOptions.DRIFT_RATE +
                                (traitDimension > 1 ? "." + i : "")),
                        new Attribute.Default<String>("value", "0.0"),
                }, true);
                writer.writeCloseTag("rate");
                writer.writeCloseTag("strictClockBranchRates");
            }
            writer.writeCloseTag("driftModels");
        }

        if (model.getJitterWindow() > 0.0) {
            StringBuilder sb = new StringBuilder(Double.toString(model.getJitterWindow()));
            for (int i = 1; i < model.getContinuousTraitCount(); i++) {
                sb.append(" ").append(Double.toString(model.getJitterWindow()));
            }
            writer.writeOpenTag("jitter", new Attribute[]{
                    new Attribute.Default<String>("window", sb.toString()),
                    new Attribute.Default<String>("duplicatesOnly", "true")
            });
            writer.writeTag("parameter", new Attribute.Default<String>("idref", "leaf." + partitionData.getName()), true);
            writer.writeCloseTag("jitter");
        }

        writer.writeOpenTag("conjugateRootPrior");

        writer.writeOpenTag("meanParameter");
        StringBuilder sb = new StringBuilder();
        for (int j = 0; j < partitionData.getPartitionSubstitutionModel().getContinuousTraitDimension(); j++) {
            if (j > 0) {
                sb.append(" ");
            }
            sb.append(0.0);
        }

        writer.writeTag("parameter", new Attribute.Default<String>("value", sb.toString()), true);
        writer.writeCloseTag("meanParameter");

        writer.writeOpenTag("priorSampleSize");
        writer.writeTag("parameter", new Attribute.Default<String>("value", "0.000001"), true);
        writer.writeCloseTag("priorSampleSize");

        writer.writeCloseTag("conjugateRootPrior");

        if (partitionData.getPartitionSubstitutionModel().getContinuousSubstModelType() != ContinuousSubstModelType.HOMOGENOUS &&
                partitionData.getPartitionSubstitutionModel().getContinuousSubstModelType() != ContinuousSubstModelType.DRIFT) {
            if (ContinuousComponentOptions.USE_ARBITRARY_BRANCH_RATE_MODEL) {
                writer.writeIDref("arbitraryBranchRates", partitionData.getName() + "." + "diffusion.branchRates");
            } else {
                writer.writeIDref("discretizedBranchRates", partitionData.getName() + "." + "diffusion.branchRates");
            }
        }

        writer.writeCloseTag(continuousDataLikelihoodParser.getParserTag());

        if (traitDimension > 1 && partitionData.getPartitionSubstitutionModel().getContinuousSubstModelType() == ContinuousSubstModelType.DRIFT) {
            writer.writeOpenTag("compoundParameter",
                    new Attribute.Default<String>("id", partitionData.getName() + "." + ContinuousComponentOptions.DRIFT_RATE));
            for (int i = 0; i < traitDimension; i++) { // todo iterate over dimension of trait
                writer.writeTag("parameter", new Attribute.Default<String>("idref", partitionData.getName() + "." + ContinuousComponentOptions.DRIFT_RATE + "." + i), true);
            }
            writer.writeCloseTag("compoundParameter");
        }

    }

    private void writeDiffusionStatistics(XMLWriter writer, AbstractPartitionData partitionData,
                                          String treeModelId, String precisionMatrixId, String traitLikelihoodId) {
        String prefix = partitionData.getName() + ".";

        if (partitionData.getTraits().size() == 2) {
            writer.writeOpenTag("correlation",
                    new Attribute[]{
                            new Attribute.Default<String>("id", prefix + "correlation"),
                            new Attribute.Default<Integer>("dimension1", 1),
                            new Attribute.Default<Integer>("dimension2", 2)
                    });
            writer.writeIDref("matrixParameter", precisionMatrixId);
            writer.writeCloseTag("correlation");
        }

        /*
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
        writer.writeIDref("parameter", prefix + "precision.col1");
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
        writer.writeIDref("parameter", prefix + "precision.col2");
        writer.writeCloseTag("subStatistic");
        writer.writeCloseTag("productStatistic");
        */

        GeneratorHelper.writeMatrixInverse(writer, prefix + "varCovar", precisionMatrixId);
        PartitionSubstitutionModel model = partitionData.getPartitionSubstitutionModel();
        String modelName = model.getName();
        switch (partitionData.getPartitionSubstitutionModel().getContinuousExtensionType()) {
            case RESIDUAL:
                GeneratorHelper.writeMatrixInverse(writer,
                        repeatedMeasuresTraitDataModelParser.getBeautiParameterIDProvider(
                                RepeatedMeasuresTraitDataModelParser.EXTENSION_VARIANCE).getId(modelName),
                        repeatedMeasuresTraitDataModelParser.getBeautiParameterIDProvider(
                                RepeatedMeasuresTraitDataModelParser.EXTENSION_PRECISION).getId(modelName)
                );
                break;
            case LATENT_FACTORS:
            case NONE:
                // do nothing
                break;
            default:
                throw new IllegalArgumentException("Unknown extension type");
        }

        writer.writeOpenTag(TreeDataContinuousDiffusionStatistic.CONTINUOUS_DIFFUSION_STATISTIC, (partitionData.getPartitionSubstitutionModel().isLatitudeLongitude() ?
                new Attribute[]{
                        new Attribute.Default<String>("id", prefix + "diffusionRate"),
                        new Attribute.Default<String>("traitName", partitionData.getName()),
                        new Attribute.Default<String>("displacementScheme", "greatCircleDistance"),
                        new Attribute.Default<String>("scalingScheme", "dependent"),
                        new Attribute.Default<String>("weightingScheme", "weighted")
                } :
                new Attribute[]{
                        new Attribute.Default<String>("id", prefix + "diffusionRate"),
                        new Attribute.Default<String>("traitName", partitionData.getName()),
                        new Attribute.Default<String>("displacementScheme", "linear"),
                        new Attribute.Default<String>("scalingScheme", "dependent"),
                        new Attribute.Default<String>("weightingScheme", "weighted")
                }));
        continuousDataLikelihoodParser.writeIDrefFromID(writer, traitLikelihoodId);
        writer.writeCloseTag(TreeDataContinuousDiffusionStatistic.CONTINUOUS_DIFFUSION_STATISTIC);
    }

    private void writePrecisionGibbsOperators(XMLWriter writer,
                                              ContinuousComponentOptions component) {

        for (AbstractPartitionData partitionData : component.getOptions().getDataPartitions(ContinuousDataType.INSTANCE)) {

            switch (partitionData.getPartitionSubstitutionModel().getContinuousExtensionType()) {
                case NONE:
                    writePrecisionGibbsOperator(writer, component, partitionData, ContinuousModelExtensionType.NONE);
                    break;
                case RESIDUAL:
                    writePrecisionGibbsOperator(writer, component, partitionData, ContinuousModelExtensionType.NONE);
                    writePrecisionGibbsOperator(writer, component, partitionData,
                            ContinuousModelExtensionType.RESIDUAL);
                    break;
                case LATENT_FACTORS:
                    writeLatentFactorOperators(writer, component, partitionData);
                    break;
                default:
                    throw new IllegalArgumentException("Unknown extension type");
            }
        }
    }

    private void writeLatentFactorOperators(final XMLWriter writer,
                                            final ContinuousComponentOptions component,
                                            final AbstractPartitionData partitionData
    ) {
        PartitionSubstitutionModel model = partitionData.getPartitionSubstitutionModel();
        XMLWriterObject loadingsGibbsOp = new XMLWriterObject(
                "loadingsGibbsOperator",
                null,
                new XMLWriterObject[]{
                        XMLWriterObject.refObj("integratedFactorModel",
                                integratedFactorsParser.getId(model.getName())),
                        XMLWriterObject.refObj(continuousDataLikelihoodParser.getParserTag(),
                                continuousDataLikelihoodParser.getId(partitionData.getName())),
                        XMLWriterObject.refObj("distributionLikelihood",
                                integratedFactorsParser.getBeautiParameterIDProvider(
                                        "loadings").getPriorId(model.getName()))
                },
                new Attribute[]{
                        new Attribute.Default("weight", 3),
                        new Attribute.Default("newMode", true),
                        new Attribute.Default("sparsity", "upperTriangular")
                }
        );
        loadingsGibbsOp.writeOrReference(writer);

        XMLWriterObject precisionGibbsOp = new XMLWriterObject(
                "normalGammaPrecisionGibbsOperator",
                null,
                new XMLWriterObject[]{
                        new XMLWriterObject("prior", XMLWriterObject.refObj("gammaPrior",
                                integratedFactorsParser.getBeautiParameterIDProvider(
                                        "precision").getPriorId(model.getName()))),
                        new XMLWriterObject("normalExtension",
                                null,
                                new XMLWriterObject[]{
                                        XMLWriterObject.refObj("integratedFactorModel",
                                                integratedFactorsParser.getId(model.getName())),
                                        XMLWriterObject.refObj(continuousDataLikelihoodParser.getParserTag(),
                                                continuousDataLikelihoodParser.getId(partitionData.getName()))
                                },
                                new Attribute[]{
                                        new Attribute.Default("treeTraitName", partitionData.getName())
                                })
                },
                new Attribute[]{
                        new Attribute.Default("weight", 1.0)
                }
        );
        precisionGibbsOp.writeOrReference(writer);
    }

    private void writePrecisionGibbsOperator(final XMLWriter writer,
                                             final ContinuousComponentOptions component,
                                             AbstractPartitionData partitionData,
                                             ContinuousModelExtensionType extensionType
    ) {
        writer.writeOpenTag(ContinuousComponentOptions.PRECISION_GIBBS_OPERATOR,
                new Attribute[]{
                        new Attribute.Default<String>("weight", "" + partitionData.getTraits().size())
                });

        String wrapperName;

        String priorId;
        switch (extensionType) {
            case NONE:
                wrapperName = WishartStatisticsWrapper.PARSER_NAME;
                priorId = partitionData.getPartitionSubstitutionModel().getName() + ".precisionPrior";
                break;
            case RESIDUAL:
                wrapperName = RepeatedMeasuresWishartStatistics.RM_WISHART_STATISTICS;
                PartitionSubstitutionModel model = partitionData.getPartitionSubstitutionModel();

                priorId = repeatedMeasuresTraitDataModelParser.getBeautiParameterIDProvider("extensionPrecision").getPriorId(model.getName());
                break;
            default:
                throw new IllegalArgumentException("Unknown or unsupported extension type");
        }

        writer.writeOpenTag(wrapperName,
                new Attribute[]{
                        new Attribute.Default<String>("traitName", "" + partitionData.getName())
                });
        continuousDataLikelihoodParser.writeIDrefFromName(writer, partitionData.getName());

        switch (extensionType) {
            case NONE:
                break;
            case RESIDUAL:
                writer.writeIDref(
                        RepeatedMeasuresTraitDataModelParser.REPEATED_MEASURES_MODEL,
                        repeatedMeasuresTraitDataModelParser.getId(partitionData.getPartitionSubstitutionModel().getName())
                );
                break;
            default:
                throw new IllegalArgumentException("Unknown or unsupported extension type");
        }

        writer.writeCloseTag(wrapperName);
        writer.writeIDref("multivariateWishartPrior", priorId);
        writer.writeCloseTag(ContinuousComponentOptions.PRECISION_GIBBS_OPERATOR);

    }

    private void writeParameterIdRefs(final XMLWriter writer, final ContinuousComponentOptions component) {
        for (AbstractPartitionData partitionData : component.getOptions().getDataPartitions(ContinuousDataType.INSTANCE)) {
            PartitionSubstitutionModel model = partitionData.getPartitionSubstitutionModel();
            String prefix = partitionData.getName() + ".";

            if (model.getContinuousExtensionType() != ContinuousModelExtensionType.LATENT_FACTORS) {
                writer.writeIDref("matrixParameter", model.getName() + ".precision");

                if (partitionData.getTraits().size() == 2) {
                    writer.writeIDref("correlation", prefix + "correlation");
//            writer.writeIDref("treeLengthStatistic", prefix + "treeLength");
//            writer.writeIDref("productStatistic", prefix + "treeLengthPrecision1");
//            writer.writeIDref("productStatistic", prefix + "treeLengthPrecision2");
                }
                writer.writeIDref("matrixInverse", prefix + "varCovar");
            }

            writer.writeIDref(TreeDataContinuousDiffusionStatistic.CONTINUOUS_DIFFUSION_STATISTIC, prefix + "diffusionRate");

            if (partitionData.getPartitionSubstitutionModel().getContinuousSubstModelType() == ContinuousSubstModelType.GAMMA_RRW) {
                writer.writeIDref(ParameterParser.PARAMETER, prefix + ContinuousComponentOptions.HALF_DF);
            } else if (partitionData.getPartitionSubstitutionModel().getContinuousSubstModelType() == ContinuousSubstModelType.LOGNORMAL_RRW) {
                writer.writeIDref(ParameterParser.PARAMETER, prefix + ContinuousComponentOptions.STDEV);
            } else if (partitionData.getPartitionSubstitutionModel().getContinuousSubstModelType() == ContinuousSubstModelType.DRIFT) {
                writer.writeIDref(ParameterParser.PARAMETER, prefix + ContinuousComponentOptions.DRIFT_RATE);
            }

            if (component.useLambda(model)) {
                writer.writeIDref("parameter", model.getName() + "." + ContinuousComponentOptions.LAMBDA);
            }

            switch (partitionData.getPartitionSubstitutionModel().getContinuousExtensionType()) {
                case RESIDUAL:
                    writer.writeIDref("matrixParameter",
                            repeatedMeasuresTraitDataModelParser.getBeautiParameterIDProvider(
                                    "extensionPrecision").getId(model.getName()));
                    writer.writeIDref("matrixInverse",
                            repeatedMeasuresTraitDataModelParser.getBeautiParameterIDProvider(
                                    RepeatedMeasuresTraitDataModelParser.EXTENSION_VARIANCE).getId(model.getName())
                    );
                    break;
                case LATENT_FACTORS:
                    writer.writeIDref("matrixParameter",
                            integratedFactorsParser.getBeautiParameterIDProvider(
                                    "loadings").getId(model.getName()));
                    writer.writeIDref("parameter",
                            integratedFactorsParser.getBeautiParameterIDProvider(
                                    "precision").getId(model.getName()));
                    break;
                case NONE:
                    break;
                default:
                    throw new IllegalArgumentException("Unknown extension type");
            }
        }
    }

    private void writeMultivariatePriors(XMLWriter writer,
                                         ContinuousComponentOptions component) {

        if (ContinuousComponentOptions.USE_ARBITRARY_BRANCH_RATE_MODEL) {
            for (AbstractPartitionData partitionData : component.getOptions().getDataPartitions(ContinuousDataType.INSTANCE)) {
                PartitionSubstitutionModel model = partitionData.getPartitionSubstitutionModel();
                if (model.getContinuousSubstModelType() != ContinuousSubstModelType.HOMOGENOUS &&
                        model.getContinuousSubstModelType() != ContinuousSubstModelType.DRIFT) {
                    writer.writeIDref("distributionLikelihood", partitionData.getName() + "." + "diffusion.prior");
                }
            }
        }
        for (AbstractPartitionData partitionData : component.getOptions().getDataPartitions(ContinuousDataType.INSTANCE)) {
            writer.writeIDref("multivariateWishartPrior", partitionData.getName() + ".precisionPrior");

            PartitionSubstitutionModel model = partitionData.getPartitionSubstitutionModel();

            switch (partitionData.getPartitionSubstitutionModel().getContinuousExtensionType()) {
                case RESIDUAL:
                    writer.writeIDref("multivariateWishartPrior",
                            repeatedMeasuresTraitDataModelParser.getBeautiParameterIDProvider("extensionPrecision").getPriorId(model.getName()));
                    break;
                case LATENT_FACTORS:

                    writer.writeIDref("distributionLikelihood",
                            integratedFactorsParser.getBeautiParameterIDProvider(
                                    "loadings").getPriorId(model.getName()));
                    writer.writeIDref("gammaPrior",
                            integratedFactorsParser.getBeautiParameterIDProvider(
                                    "precision").getPriorId(model.getName()));
                    break;
                case NONE:
                    break;
                default:
                    throw new IllegalArgumentException("Unknown extension type");
            }
        }

    }

    private void writeMultivariateTreeLikelihoodIdRefs(XMLWriter writer,
                                                       ContinuousComponentOptions component) {

        for (AbstractPartitionData partitionData : component.getOptions().getDataPartitions(ContinuousDataType.INSTANCE)) {
            continuousDataLikelihoodParser.writeIDrefFromName(writer, partitionData.getName());
        }
    }

    private void writeTreeLogEntries(PartitionTreeModel treeModel, XMLWriter writer) {
        for (AbstractPartitionData partitionData : options.getDataPartitions(ContinuousDataType.INSTANCE)) {
            if (partitionData.getPartitionTreeModel() == treeModel) {
                PartitionSubstitutionModel model = partitionData.getPartitionSubstitutionModel();
                writer.writeIDref("multivariateDiffusionModel", model.getName() + ".diffusionModel");
                continuousDataLikelihoodParser.writeIDrefFromName(writer, partitionData.getName());
                if (model.getContinuousSubstModelType() != ContinuousSubstModelType.HOMOGENOUS &&
                        model.getContinuousSubstModelType() != ContinuousSubstModelType.DRIFT) {
                    writer.writeOpenTag(TreeLoggerParser.TREE_TRAIT,
                            new Attribute[]{
                                    new Attribute.Default<String>(TreeLoggerParser.NAME, "rate"),
                                    new Attribute.Default<String>(TreeLoggerParser.TAG, partitionData.getName() + ".rate"),
                            });
                    if (ContinuousComponentOptions.USE_ARBITRARY_BRANCH_RATE_MODEL) {
                        writer.writeIDref("arbitraryBranchRates", partitionData.getName() + "." + "diffusion.branchRates");
                    } else {
                        writer.writeIDref("discretizedBranchRates", partitionData.getName() + "." + "diffusion.branchRates");
                    }

                    writer.writeCloseTag(TreeLoggerParser.TREE_TRAIT);
                }
            }
        }
    }

    private static final ContinuousDataLikelihoodParser continuousDataLikelihoodParser = new ContinuousDataLikelihoodParser();
    private static final RepeatedMeasuresTraitDataModelParser repeatedMeasuresTraitDataModelParser = new RepeatedMeasuresTraitDataModelParser();
    private static final IntegratedFactorAnalysisLikelihoodParser integratedFactorsParser = new IntegratedFactorAnalysisLikelihoodParser();
}
