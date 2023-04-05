/*
 * ParameterPriorGenerator.java
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
import dr.app.beauti.types.ClockType;
import dr.app.beauti.types.PriorType;
import dr.app.beauti.util.XMLWriter;
import dr.evolution.util.Taxa;
import dr.evomodel.tree.DefaultTreeModel;
import dr.evomodelxml.tree.CTMCScalePriorParser;
import dr.evomodelxml.tree.MonophylyStatisticParser;
import dr.inference.model.ParameterParser;
import dr.inferencexml.distribution.*;
import dr.inferencexml.model.BooleanLikelihoodParser;
import dr.inferencexml.model.OneOnXPriorParser;
import dr.util.Attribute;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @author Walter Xie
 */
public class ParameterPriorGenerator extends Generator {

    //map parameters to prior IDs, for use with HMC
    private HashMap<String, String> mapParameterToPrior;

    public ParameterPriorGenerator(BeautiOptions options, ComponentFactory[] components) {
        super(options, components);
        //TODO don't like this being here, but will see how things pan out as more HMC approaches are added
        mapParameterToPrior = new HashMap<String, String>();
        mapParameterToPrior.put("skygrid.precision", "skygrid.precision.prior");
    }

    /**
     * Write the priors for each parameter
     *
     * @param writer       the writer
     */
    public void writeParameterPriors(XMLWriter writer) {
        boolean first = true;

        for (Map.Entry<Taxa, Boolean> taxaBooleanEntry : options.taxonSetsMono.entrySet()) {
            if (taxaBooleanEntry.getValue()) {
                if (first) {
                    writer.writeOpenTag(BooleanLikelihoodParser.BOOLEAN_LIKELIHOOD);
                    first = false;
                }
                final String taxaRef = "monophyly(" + taxaBooleanEntry.getKey().getId() + ")";
                writer.writeIDref(MonophylyStatisticParser.MONOPHYLY_STATISTIC, taxaRef);
            }
        }
        if (!first) {
            writer.writeCloseTag(BooleanLikelihoodParser.BOOLEAN_LIKELIHOOD);
        }

        List<Parameter> parameters = options.selectParameters();

        for (Parameter parameter : parameters) {
            if (!(parameter.priorType == PriorType.NONE_TREE_PRIOR ||
                    parameter.priorType == PriorType.NONE_FIXED ||
                    parameter.priorType == PriorType.NONE_STATISTIC)) {
                if (parameter.isCached) {
                    writeCachedParameterPrior(parameter, writer);
                    //if (parameter.priorType != PriorType.UNIFORM_PRIOR || parameter.isNodeHeight) {
                } else if (!(options.treeModelOptions.isNodeCalibrated(parameter) && parameter.isCalibratedYule)) {
                    writeParameterPrior(parameter, writer);
                }
            }
        }

    }

    private void writeCachedParameterPrior(Parameter parameter, XMLWriter writer) {
        writer.writeOpenTag(CachedDistributionLikelihoodParser.CACHED_PRIOR);

        writeParameterPrior(parameter, writer);
        writeParameterIdref(writer, parameter);

        writer.writeCloseTag(CachedDistributionLikelihoodParser.CACHED_PRIOR);
    }

    /**
     * Write the priors for each parameter
     *
     * @param parameter the parameter
     * @param writer    the writer
     */
    public void writeParameterPrior(Parameter parameter, XMLWriter writer) {

        //if models need to have a prior defined before the priors block
        /*if (reservedParameters.contains(parameter.getName())) {

            return;
        }*/

        if (mapParameterToPrior.keySet().contains(parameter.getName())) {
            writePriorIdref(writer, parameter, mapParameterToPrior.get(parameter.getName()));
            return;
        }

        if (parameter.priorType == PriorType.NONE_FIXED) {
            return;
        }

        boolean isTruncated = parameter.isTruncated;
        if (isTruncated) {
            Attribute[] attributes = null;

            if (!Double.isInfinite(parameter.truncationLower)) {
                if (!Double.isInfinite(parameter.truncationUpper)) {
                    attributes = new Attribute[]{
                            new Attribute.Default<String>(PriorParsers.LOWER, "" + parameter.truncationLower),
                            new Attribute.Default<String>(PriorParsers.UPPER, "" + parameter.truncationUpper)
                    };
                } else {
                    attributes = new Attribute[]{
                            new Attribute.Default<String>(PriorParsers.LOWER, "" + parameter.truncationLower),
                    };
                }
            } else {
                if (!Double.isInfinite(parameter.truncationUpper)) {
                    attributes = new Attribute[]{
                            new Attribute.Default<String>(PriorParsers.UPPER, "" + parameter.truncationUpper)
                    };
                } else {
                    // both are infinite so there is no trunction
                    isTruncated = false;
                }
            }

            if (isTruncated) {
                writer.writeOpenTag(PriorParsers.TRUNCATED, attributes);
            }
        }

        switch (parameter.priorType) {
            case NONE_FIXED:
                break;
            case NONE_IMPROPER:
                writer.writeComment("Improper uniform prior: " + parameter.getName());
                break;
            case DISCRETE_UNIFORM_PRIOR:
                writer.writeOpenTag(PriorParsers.DISCRETE_UNIFORM_PRIOR,
                        new Attribute[]{
                                new Attribute.Default<String>(PriorParsers.LOWER, "" + parameter.getLowerBound()),
                                new Attribute.Default<String>(PriorParsers.UPPER, "" + parameter.getUpperBound())
                        });
                writeParameterIdref(writer, parameter);
                writer.writeCloseTag(PriorParsers.DISCRETE_UNIFORM_PRIOR);
                break;
            case UNIFORM_PRIOR:
                if (parameter.isPriorImproper()) {
                    throw new IllegalArgumentException("Uniform priors cannot have infinite bounds (use 'NONE_IMPROPER')");
//                    writer.writeComment("Improper uniform prior: " + parameter.getName());
                } else {
                    writer.writeOpenTag(PriorParsers.UNIFORM_PRIOR,
                            new Attribute[]{
                                    new Attribute.Default<String>(PriorParsers.LOWER, "" + parameter.getLowerBound()),
                                    new Attribute.Default<String>(PriorParsers.UPPER, "" + parameter.getUpperBound())
                            });
                    writeParameterIdref(writer, parameter);
                    writer.writeCloseTag(PriorParsers.UNIFORM_PRIOR);
                }
                break;
            case EXPONENTIAL_PRIOR:
                writer.writeOpenTag(PriorParsers.EXPONENTIAL_PRIOR,
                        new Attribute[]{
                                new Attribute.Default<String>(PriorParsers.MEAN, "" + parameter.mean),
                                new Attribute.Default<String>(PriorParsers.OFFSET, "" + parameter.offset)
                        });
                writeParameterIdref(writer, parameter);
                writer.writeCloseTag(PriorParsers.EXPONENTIAL_PRIOR);
                break;
            case LAPLACE_PRIOR:
                writer.writeOpenTag(PriorParsers.LAPLACE_PRIOR,
                        new Attribute[]{
                                new Attribute.Default<String>(PriorParsers.MEAN, "" + parameter.mean),
                                new Attribute.Default<String>(PriorParsers.SCALE, "" + parameter.scale)
                        });
                writeParameterIdref(writer, parameter);
                writer.writeCloseTag(PriorParsers.LAPLACE_PRIOR);
                break;
            case NORMAL_PRIOR:
                writer.writeOpenTag(PriorParsers.NORMAL_PRIOR,
                        new Attribute[]{
                                new Attribute.Default<String>(PriorParsers.MEAN, "" + parameter.mean),
                                new Attribute.Default<String>(PriorParsers.STDEV, "" + parameter.stdev)
                        });
                writeParameterIdref(writer, parameter);
                writer.writeCloseTag(PriorParsers.NORMAL_PRIOR);
                break;
            case LOGNORMAL_PRIOR:
                if (parameter.isInRealSpace()) {
                    writer.writeOpenTag(PriorParsers.LOG_NORMAL_PRIOR,
                            new Attribute[]{
                                    new Attribute.Default<String>(PriorParsers.MEAN, "" + parameter.mean),
                                    new Attribute.Default<String>(PriorParsers.STDEV, "" + parameter.stdev),
                                    new Attribute.Default<String>(PriorParsers.OFFSET, "" + parameter.offset),
                            });
                } else {
                    // if the log normal parameters are not set in real space then the parameter use the mu
                    // and sigma parameters (the mean and stdev of the underlying normal).
                    writer.writeOpenTag(PriorParsers.LOG_NORMAL_PRIOR,
                            new Attribute[]{
                                    new Attribute.Default<String>(PriorParsers.MU, "" + parameter.mean),
                                    new Attribute.Default<String>(PriorParsers.SIGMA, "" + parameter.stdev),
                                    new Attribute.Default<String>(PriorParsers.OFFSET, "" + parameter.offset),
                            });
                }
                writeParameterIdref(writer, parameter);
                writer.writeCloseTag(PriorParsers.LOG_NORMAL_PRIOR);
                break;
            case GAMMA_PRIOR:
                writer.writeOpenTag(PriorParsers.GAMMA_PRIOR,
                        new Attribute[]{
                                new Attribute.Default<String>(PriorParsers.SHAPE, "" + parameter.shape),
                                new Attribute.Default<String>(PriorParsers.SCALE, "" + parameter.scale),
                                new Attribute.Default<String>(PriorParsers.OFFSET, "" + parameter.offset)
                        });
                writeParameterIdref(writer, parameter);
                writer.writeCloseTag(PriorParsers.GAMMA_PRIOR);
                break;
            case INVERSE_GAMMA_PRIOR:
                writer.writeOpenTag(PriorParsers.INVGAMMA_PRIOR,
                        new Attribute[]{
                                new Attribute.Default<String>(PriorParsers.SHAPE, "" + parameter.shape),
                                new Attribute.Default<String>(PriorParsers.SCALE, "" + parameter.scale),
                                new Attribute.Default<String>(PriorParsers.OFFSET, "" + parameter.offset)
                        });
                writeParameterIdref(writer, parameter);
                writer.writeCloseTag(PriorParsers.INVGAMMA_PRIOR);
                break;
            case ONE_OVER_X_PRIOR:
                writer.writeOpenTag(OneOnXPriorParser.ONE_ONE_X_PRIOR);
                writeParameterIdref(writer, parameter);
                writer.writeCloseTag(OneOnXPriorParser.ONE_ONE_X_PRIOR);
                break;
            case POISSON_PRIOR:
                writer.writeOpenTag(PriorParsers.POISSON_PRIOR,
                        new Attribute[]{
                                new Attribute.Default<String>(PriorParsers.MEAN, "" + parameter.mean),
                                new Attribute.Default<String>(PriorParsers.OFFSET, "" + parameter.offset)
                        });
                writeParameterIdref(writer, parameter);
                writer.writeCloseTag(PriorParsers.POISSON_PRIOR);
                break;
            case BETA_PRIOR:
                writer.writeOpenTag(PriorParsers.BETA_PRIOR,
                        new Attribute[]{
                                new Attribute.Default<String>(PriorParsers.SHAPE, "" + parameter.shape),
                                new Attribute.Default<String>(PriorParsers.SHAPEB, "" + parameter.shapeB),
                                new Attribute.Default<String>(PriorParsers.OFFSET, "" + parameter.offset)
                        });
                writeParameterIdref(writer, parameter);
                writer.writeCloseTag(PriorParsers.BETA_PRIOR);
                break;
            case CTMC_RATE_REFERENCE_PRIOR:
                writer.writeOpenTag(CTMCScalePriorParser.MODEL_NAME);
                writer.writeOpenTag(CTMCScalePriorParser.SCALEPARAMETER);
                writeParameterIdref(writer, parameter);
                writer.writeCloseTag(CTMCScalePriorParser.SCALEPARAMETER);
                // Find correct tree for this rate parameter

                PartitionTreeModel treeModel = null;

                if (parameter.getTaxonSet() != null) {
                    treeModel = options.taxonSetsTreeModel.get(parameter.getTaxonSet());
                } else {
                    for (PartitionClockModel pcm : options.getPartitionClockModels()) {
                        if (pcm.performModelAveraging()) {
                            treeModel = pcm.getPartitionTreeModel();
                        } else if (pcm.getClockRateParameter() == parameter) {
                            for (AbstractPartitionData pd : options.getDataPartitions(pcm)) {
                                treeModel = pd.getPartitionTreeModel();
                                break; // todo - This breaks after the first iteration. Why a loop?
                            }
                        }
                    }
                }
                if (treeModel == null) {
                    throw new IllegalArgumentException("No tree model found for clock model");
                }

                PartitionClockModel pcm = (PartitionClockModel)parameter.getOptions();
                if (pcm.getClockType() == ClockType.FIXED_LOCAL_CLOCK) {
                    if (parameter.getTaxonSet() != null) {
                        writer.writeIDref("taxa", parameter.getTaxonSet().getId());
                    } else {
                        writer.writeOpenTag("taxa");
                        writer.writeIDref("taxa", "taxa");
                        writer.writeOpenTag("exclude");
                        for (Taxa taxonSet : options.taxonSets) {
                            if (options.taxonSetsMono.get(taxonSet)) {
                                writer.writeIDref("taxa", taxonSet.getId());
                            }
                        }
                        writer.writeCloseTag("exclude");
                        writer.writeCloseTag("taxa");

                    }
                }

                writer.writeIDref(DefaultTreeModel.TREE_MODEL, treeModel.getPrefix() + DefaultTreeModel.TREE_MODEL);
                writer.writeCloseTag(CTMCScalePriorParser.MODEL_NAME);
                break;
            case NORMAL_HPM_PRIOR:
            case LOGNORMAL_HPM_PRIOR:
                // Do nothing, densities are already in a distributionLikelihood
                break;
            case DIRICHLET_PRIOR:
                writer.writeOpenTag(PriorParsers.DIRICHLET_PRIOR,
                        new Attribute[]{
                                new Attribute.Default<String>(PriorParsers.ALPHA, "1.0"),
                                new Attribute.Default<Double>(PriorParsers.SUMS_TO, parameter.maintainedSum)
                        }
                );
                writeParameterIdref(writer, parameter);
                writer.writeCloseTag(PriorParsers.DIRICHLET_PRIOR);
                break;
            default:
                throw new IllegalArgumentException("Unknown priorType");
        }

        if (isTruncated) {
            writer.writeCloseTag(PriorParsers.TRUNCATED);
        }

    }

    private void writeParameterIdref(XMLWriter writer, Parameter parameter) {
        if (parameter.isStatistic) {
            writer.writeIDref("statistic", parameter.getName());
        } else {
            writer.writeIDref(ParameterParser.PARAMETER, parameter.getName());
        }
    }

    private void writePriorIdref(XMLWriter writer, Parameter parameter, String priorID) {
        switch (parameter.priorType) {
            case GAMMA_PRIOR:
                writer.writeIDref(PriorParsers.GAMMA_PRIOR, priorID);
                break;
            case LOGNORMAL_PRIOR:
                writer.writeIDref(PriorParsers.LOG_NORMAL_PRIOR, priorID);
                break;
            case EXPONENTIAL_PRIOR:
                writer.writeIDref(PriorParsers.EXPONENTIAL_PRIOR, priorID);
                break;
            case INVERSE_GAMMA_PRIOR:
                writer.writeIDref(PriorParsers.INVGAMMA_PRIOR_CORRECT, priorID);
                break;
            default:
                throw new IllegalArgumentException("Unknown or invalid prior defined on " + parameter.getName());
        }
    }

}
