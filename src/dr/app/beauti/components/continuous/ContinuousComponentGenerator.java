/*
 * ContinuousComponentGenerator.java
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

package dr.app.beauti.components.continuous;

import dr.app.beauti.generator.BaseComponentGenerator;
import dr.app.beauti.options.*;
import dr.app.beauti.types.OperatorType;
import dr.app.beauti.util.XMLWriter;
import dr.evolution.datatype.ContinuousDataType;
import dr.evolution.util.Taxon;
import dr.evomodel.continuous.ContinuousDiffusionStatistic;
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
                // the RRW operators are added to the operator list
//                writeRRWOperators(writer, component);
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
                            new Attribute.Default<String>("id", precisionMatrixId + ".col" + (i + 1)),
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

            if (model.getContinuousSubstModelType() != ContinuousSubstModelType.HOMOGENOUS &&
                    model.getContinuousSubstModelType() != ContinuousSubstModelType.DRIFT) {
                writeRelaxedBranchRateModel(writer, partitionData, treeModelId);
            }

            writer.writeBlankLine();

            writeMultivariateTreeLikelihood(writer, partitionData, diffusionModelId, treeModelId);

            String precisionMatrixId = model.getName() + ".precision";

            writeDiffusionStatistics(writer, partitionData, treeModelId, precisionMatrixId,
                    partitionData.getName() + ".traitLikelihood");
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

        if (partitionData.getPartitionSubstitutionModel().getContinuousSubstModelType() == ContinuousSubstModelType.LOGNORMAL_RRW) {
            writer.writeOpenTag("logNormalDistributionModel",
                    new Attribute[]{ new Attribute.Default<String>("meanInRealSpace", "true") });

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

    private void writeMultivariateTreeLikelihood(XMLWriter writer,
                                                 AbstractPartitionData partitionData,
                                                 String diffusionModelId,
                                                 String treeModelId) {

        int traitDimension = 1; // todo - set this to trait dimension
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


        ContinuousComponentOptions component = (ContinuousComponentOptions) options
                .getComponentOptions(ContinuousComponentOptions.class);

        PartitionSubstitutionModel model = partitionData.getPartitionSubstitutionModel();

        if (component.useLambda(model)) {
            writer.writeOpenTag("transformedTreeModel");
            writer.writeIDref("treeModel", treeModelId);
            writer.writeTag("parameter", new Attribute[] {
                    new Attribute.Default<String>("id", partitionData.getName() + "." + ContinuousComponentOptions.LAMBDA),
                    new Attribute.Default<String>("value", "0.5"),
                    new Attribute.Default<String>("lower", "0.0"),
                    new Attribute.Default<String>("upper", "1.0")
            }, true);
            writer.writeCloseTag("transformedTreeModel");
        } else {
            writer.writeIDref("treeModel", treeModelId);
        }

        writer.writeOpenTag("traitParameter");
        writer.writeTag("parameter", new Attribute.Default<String>("id", "leaf." + partitionData.getName()), true);
        writer.writeCloseTag("traitParameter");

        if (partitionData.getPartitionSubstitutionModel().getContinuousSubstModelType() == ContinuousSubstModelType.DRIFT) {
            writer.writeOpenTag("driftModels");
            for (int i = 0; i < traitDimension; i++) { 
                writer.writeOpenTag("strictClockBranchRates");
                writer.writeTag("parameter", new Attribute[]{
                        new Attribute.Default<String>("id", partitionData.getName() + "." + ContinuousComponentOptions.DRIFT_RATE +
                                (traitDimension > 1 ? "." + i : "")),
                        new Attribute.Default<String>("value", "0.0"),
                }, true);
                writer.writeCloseTag("strictClockBranchRates");
            }
            writer.writeCloseTag("driftModels");
        }

        if (model.getJitterWindow() > 0.0) {
            StringBuilder sb = new StringBuilder(Double.toString(model.getJitterWindow()));
            for (int i = 1; i < model.getContinuousTraitCount(); i++) {
                sb.append(" ").append(Double.toString(model.getJitterWindow()));
            }
            writer.writeOpenTag("jitter", new Attribute[] {
                    new Attribute.Default<String>("window", sb.toString()),
                    new Attribute.Default<String>("duplicatesOnly", "true")
            });
            writer.writeTag("parameter", new Attribute.Default<String>("idref", "leaf." + partitionData.getName()), true);
            writer.writeCloseTag("jitter");
        }

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
        writer.writeTag("parameter", new Attribute.Default<String>("value", "0.000001"), true);
        writer.writeCloseTag("priorSampleSize");

        writer.writeCloseTag("conjugateRootPrior");

        if (partitionData.getPartitionSubstitutionModel().getContinuousSubstModelType() != ContinuousSubstModelType.HOMOGENOUS) {
            writer.writeIDref("discretizedBranchRates", partitionData.getName() + "." + "diffusionRates");
        }

        writer.writeCloseTag("multivariateTraitLikelihood");

        if (traitDimension > 1) {
            writer.writeOpenTag("compoundParameter",
                    new Attribute.Default<String>("id", partitionData.getName() + "." + ContinuousComponentOptions.DRIFT_RATE));
            for (int i = 0; i < traitDimension; i++) { // todo iterate over dimension of trait
                writer.writeTag("parameter", new Attribute.Default<String>("idref", partitionData.getName() + "." + ContinuousComponentOptions.DRIFT_RATE + "." + i), true);
            }
            writer.writeCloseTag("priorSampleSize");
        }

    }

    private void writeDiffusionStatistics(XMLWriter writer, AbstractPartitionData partitionData,
                                          String treeModelId, String precisionMatrixId, String traitLikelihoodId) {
        String prefix = partitionData.getName() + ".";

        if (partitionData.getTraits().size() == 2) {
            writer.writeOpenTag("correlation",
                    new Attribute[] {
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

        writer.writeOpenTag("matrixInverse",
                new Attribute[] {
                        new Attribute.Default<String>("id", prefix + "varCovar")
                });
        writer.writeIDref("matrixParameter", precisionMatrixId);
        writer.writeCloseTag("matrixInverse");

        writer.writeOpenTag(ContinuousDiffusionStatistic.CONTINUOUS_DIFFUSION_STATISTIC, (partitionData.getPartitionSubstitutionModel().isLatitudeLongitude() ?
                new Attribute[] {
                        new Attribute.Default<String>("id", prefix + "diffusionRate"),
                        new Attribute.Default<String>("greatCircleDistance", "true")
                } :
                new Attribute[] {
                        new Attribute.Default<String>("id", prefix + "diffusionRate"),
                }));
        writer.writeIDref("multivariateTraitLikelihood", traitLikelihoodId);
        writer.writeCloseTag(ContinuousDiffusionStatistic.CONTINUOUS_DIFFUSION_STATISTIC);
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
        writer.writeOpenTag(ContinuousComponentOptions.PRECISION_GIBBS_OPERATOR,
                new Attribute[] {
                        new Attribute.Default<String>("weight", "" + partitionData.getTraits().size())
                });
        writer.writeIDref("multivariateTraitLikelihood", partitionData.getName() + ".traitLikelihood");
        writer.writeIDref("multivariateWishartPrior", partitionData.getPartitionSubstitutionModel().getName() + ".precisionPrior");
        writer.writeCloseTag(ContinuousComponentOptions.PRECISION_GIBBS_OPERATOR);

    }

    private void writeParameterIdRefs(final XMLWriter writer, final ContinuousComponentOptions component) {
        for (AbstractPartitionData partitionData : component.getOptions().getDataPartitions(ContinuousDataType.INSTANCE)) {
            PartitionSubstitutionModel model = partitionData.getPartitionSubstitutionModel();
            writer.writeIDref("matrixParameter", model.getName() + ".precision");

            String prefix = partitionData.getName() + ".";
            if (partitionData.getTraits().size() == 2) {
                writer.writeIDref("correlation", prefix + "correlation");
//            writer.writeIDref("treeLengthStatistic", prefix + "treeLength");
//            writer.writeIDref("productStatistic", prefix + "treeLengthPrecision1");
//            writer.writeIDref("productStatistic", prefix + "treeLengthPrecision2");
            }
            writer.writeIDref("matrixInverse", prefix + "varCovar");
            writer.writeIDref(ContinuousDiffusionStatistic.CONTINUOUS_DIFFUSION_STATISTIC, prefix + "diffusionRate");

            if (component.useLambda(model)) {
                writer.writeIDref("parameter", model.getName() + "." + ContinuousComponentOptions.LAMBDA);
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
