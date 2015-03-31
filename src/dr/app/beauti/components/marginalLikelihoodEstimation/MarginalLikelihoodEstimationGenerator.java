/*
 * MarginalLikelihoodEstimationGenerator.java
 *
 * Copyright (C) 2002-2013 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.app.beauti.components.marginalLikelihoodEstimation;

import dr.app.beauti.generator.BaseComponentGenerator;
import dr.app.beauti.options.BeautiOptions;
import dr.app.beauti.options.Parameter;
import dr.app.beauti.types.PriorType;
import dr.app.beauti.types.TreePriorType;
import dr.app.beauti.util.XMLWriter;
import dr.evolution.util.Units;
import dr.evomodel.coalescent.ExponentialGrowthModel;
import dr.evomodel.tree.TreeModel;
import dr.evomodelxml.TreeWorkingPriorParsers;
import dr.evomodelxml.coalescent.*;
import dr.inference.mcmc.MarginalLikelihoodEstimator;
import dr.inference.model.ParameterParser;
import dr.inference.model.PathLikelihood;
import dr.inference.trace.GeneralizedSteppingStoneSamplingAnalysis;
import dr.inference.trace.PathSamplingAnalysis;
import dr.inference.trace.SteppingStoneSamplingAnalysis;
import dr.inferencexml.distribution.WorkingPriorParsers;
import dr.inferencexml.model.CompoundLikelihoodParser;
import dr.util.Attribute;
import dr.xml.XMLParser;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Andrew Rambaut
 * @version $Id$
 */
public class MarginalLikelihoodEstimationGenerator extends BaseComponentGenerator {

    public static final boolean DEBUG = true;

    private BeautiOptions beautiOptions = null;

    MarginalLikelihoodEstimationGenerator(final BeautiOptions options) {
        super(options);
        this.beautiOptions = options;
    }

    public boolean usesInsertionPoint(final InsertionPoint point) {
        MarginalLikelihoodEstimationOptions component = (MarginalLikelihoodEstimationOptions) options.getComponentOptions(MarginalLikelihoodEstimationOptions.class);

        if (!component.performMLE && !component.performMLEGSS) {
            return false;
        }

        switch (point) {
            case AFTER_MCMC:
                return true;
        }
        return false;
    }

    protected void generate(final InsertionPoint point, final Object item, final String prefix, final XMLWriter writer) {
        MarginalLikelihoodEstimationOptions component = (MarginalLikelihoodEstimationOptions) options.getComponentOptions(MarginalLikelihoodEstimationOptions.class);

        switch (point) {
            case AFTER_MCMC:
                writeMLE(writer, component);
                break;
            default:
                throw new IllegalArgumentException("This insertion point is not implemented for " + this.getClass().getName());
        }

    }

    protected String getCommentLabel() {
        return "Marginal Likelihood Estimator";
    }

    /**
     * Write the marginalLikelihoodEstimator, pathSamplingAnalysis and steppingStoneSamplingAnalysis blocks.
     *
     * @param writer XMLWriter
     */
    public void writeMLE(XMLWriter writer, MarginalLikelihoodEstimationOptions options) {

        if (options.performMLE) {

            writer.writeComment("Define marginal likelihood estimator (PS/SS) settings");

            List<Attribute> attributes = new ArrayList<Attribute>();
            //attributes.add(new Attribute.Default<String>(XMLParser.ID, "mcmc"));
            attributes.add(new Attribute.Default<Integer>("chainLength", options.mleChainLength));
            attributes.add(new Attribute.Default<Integer>("pathSteps", options.pathSteps));
            attributes.add(new Attribute.Default<String>("pathScheme", options.pathScheme));
            if (!options.pathScheme.equals("linear")) {
                attributes.add(new Attribute.Default<Double>("alpha", options.schemeParameter));
            }

            writer.writeOpenTag(MarginalLikelihoodEstimator.MARGINAL_LIKELIHOOD_ESTIMATOR, attributes);

            writer.writeOpenTag("samplers");
            writer.writeIDref("mcmc", "mcmc");
            writer.writeCloseTag("samplers");

            attributes = new ArrayList<Attribute>();
            attributes.add(new Attribute.Default<String>(XMLParser.ID, "pathLikelihood"));
            writer.writeOpenTag(PathLikelihood.PATH_LIKELIHOOD, attributes);
            writer.writeOpenTag(PathLikelihood.SOURCE);
            writer.writeIDref(CompoundLikelihoodParser.POSTERIOR, CompoundLikelihoodParser.POSTERIOR);
            writer.writeCloseTag(PathLikelihood.SOURCE);
            writer.writeOpenTag(PathLikelihood.DESTINATION);
            writer.writeIDref(CompoundLikelihoodParser.PRIOR, CompoundLikelihoodParser.PRIOR);
            writer.writeCloseTag(PathLikelihood.DESTINATION);
            writer.writeCloseTag(PathLikelihood.PATH_LIKELIHOOD);

            attributes = new ArrayList<Attribute>();
            attributes.add(new Attribute.Default<String>(XMLParser.ID, "MLELog"));
            attributes.add(new Attribute.Default<Integer>("logEvery", options.mleLogEvery));
            attributes.add(new Attribute.Default<String>("fileName", options.mleFileName));
            writer.writeOpenTag("log", attributes);
            writer.writeIDref("pathLikelihood", "pathLikelihood");
            writer.writeCloseTag("log");

            writer.writeCloseTag(MarginalLikelihoodEstimator.MARGINAL_LIKELIHOOD_ESTIMATOR);

            writer.writeComment("Path sampling estimator from collected samples");
            attributes = new ArrayList<Attribute>();
            attributes.add(new Attribute.Default<String>("fileName", options.mleFileName));
            writer.writeOpenTag(PathSamplingAnalysis.PATH_SAMPLING_ANALYSIS, attributes);
            writer.writeTag("likelihoodColumn", new Attribute.Default<String>("name", "pathLikelihood.delta"), true);
            writer.writeTag("thetaColumn", new Attribute.Default<String>("name", "pathLikelihood.theta"), true);
            writer.writeCloseTag(PathSamplingAnalysis.PATH_SAMPLING_ANALYSIS);

            writer.writeComment("Stepping-stone sampling estimator from collected samples");
            attributes = new ArrayList<Attribute>();
            attributes.add(new Attribute.Default<String>("fileName", options.mleFileName));
            writer.writeOpenTag(SteppingStoneSamplingAnalysis.STEPPING_STONE_SAMPLING_ANALYSIS, attributes);
            writer.writeTag("likelihoodColumn", new Attribute.Default<String>("name", "pathLikelihood.delta"), true);
            writer.writeTag("thetaColumn", new Attribute.Default<String>("name", "pathLikelihood.theta"), true);
            writer.writeCloseTag(SteppingStoneSamplingAnalysis.STEPPING_STONE_SAMPLING_ANALYSIS);

        } else if (options.performMLEGSS) {

            //First define necessary components for the tree working prior
            if (options.choiceTreeWorkingPrior.equals("Product of exponential distributions")) {
                //more general product of exponentials needs to be constructed

                if (DEBUG) {
                    System.err.println("productOfExponentials selected: " + options.choiceTreeWorkingPrior);
                }

                List<Attribute> attributes = new ArrayList<Attribute>();
                attributes.add(new Attribute.Default<String>(XMLParser.ID, "exponentials"));
                attributes.add(new Attribute.Default<String>("fileName", beautiOptions.logFileName));
                attributes.add(new Attribute.Default<String>("burnin", "" + beautiOptions.chainLength*0.10));
                attributes.add(new Attribute.Default<String>("parameterColumn", "coalescentEventsStatistic"));
                attributes.add(new Attribute.Default<String>("dimension", "" + (beautiOptions.taxonList.getTaxonCount()-1)));

                writer.writeOpenTag(TreeWorkingPriorParsers.PRODUCT_OF_EXPONENTIALS_POSTERIOR_MEANS_LOESS, attributes);
                writer.writeTag(TreeModel.TREE_MODEL, new Attribute.Default<String>(XMLParser.ID, TreeModel.TREE_MODEL), true);
                writer.writeCloseTag(TreeWorkingPriorParsers.PRODUCT_OF_EXPONENTIALS_POSTERIOR_MEANS_LOESS);

            } else {
                //matching coalescent model has to be constructed
                //getting the coalescent model
                if (DEBUG) {
                    System.err.println("matching coalescent model selected: " + options.choiceTreeWorkingPrior);
                    System.err.println(beautiOptions.getPartitionTreePriors().get(0).getNodeHeightPrior());
                }
                /*for (PartitionTreePrior prior : options.getPartitionTreePriors()) {
                    treePriorGenerator.writeTreePriorModel(prior, writer);
                    writer.writeText("");
                }*/
                //TODO: extend for more than 1 coalescent model?
                TreePriorType nodeHeightPrior = beautiOptions.getPartitionTreePriors().get(0).getNodeHeightPrior();

                switch (nodeHeightPrior) {
                    case CONSTANT:

                        writer.writeComment("A working prior for the constant population size model.");
                        writer.writeOpenTag(
                                ConstantPopulationModelParser.CONSTANT_POPULATION_MODEL,
                                new Attribute[]{
                                        new Attribute.Default<String>(XMLParser.ID, modelPrefix + "constantReference"),
                                        new Attribute.Default<String>("units", Units.Utils.getDefaultUnitName(beautiOptions.units))
                                }
                        );

                        writer.writeOpenTag(ConstantPopulationModelParser.POPULATION_SIZE);
                        writeParameter("constantReference.popSize", "constant.popSize", beautiOptions.logFileName, (int) (options.mleChainLength * 0.10), writer);
                        writer.writeCloseTag(ConstantPopulationModelParser.POPULATION_SIZE);
                        writer.writeCloseTag(ConstantPopulationModelParser.CONSTANT_POPULATION_MODEL);

                        writer.writeComment("A working prior for the coalescent.");
                        writer.writeOpenTag(
                                CoalescentLikelihoodParser.COALESCENT_LIKELIHOOD,
                                new Attribute[]{
                                        new Attribute.Default<String>(XMLParser.ID, modelPrefix + "coalescentReference")
                                }
                        );
                        writer.writeOpenTag(CoalescentLikelihoodParser.MODEL);
                        writer.writeIDref(ConstantPopulationModelParser.CONSTANT_POPULATION_MODEL, beautiOptions.getPartitionTreePriors().get(0).getPrefix() + "constantReference");
                        writer.writeCloseTag(CoalescentLikelihoodParser.MODEL);
                        writer.writeOpenTag(CoalescentLikelihoodParser.POPULATION_TREE);
                        writer.writeIDref(TreeModel.TREE_MODEL, modelPrefix + TreeModel.TREE_MODEL);
                        writer.writeCloseTag(CoalescentLikelihoodParser.POPULATION_TREE);
                        writer.writeCloseTag(CoalescentLikelihoodParser.COALESCENT_LIKELIHOOD);

                        break;

                    case EXPONENTIAL:



                        writer.writeComment("A working prior for the coalescent.");
                        writer.writeOpenTag(
                                CoalescentLikelihoodParser.COALESCENT_LIKELIHOOD,
                                new Attribute[]{
                                        new Attribute.Default<String>(XMLParser.ID, modelPrefix + "coalescentReference")
                                }
                        );
                        writer.writeOpenTag(CoalescentLikelihoodParser.MODEL);
                        writer.writeIDref(ExponentialGrowthModelParser.EXPONENTIAL_GROWTH_MODEL, beautiOptions.getPartitionTreePriors().get(0).getPrefix() + "constantReference");
                        writer.writeCloseTag(CoalescentLikelihoodParser.MODEL);
                        writer.writeOpenTag(CoalescentLikelihoodParser.POPULATION_TREE);
                        writer.writeIDref(TreeModel.TREE_MODEL, modelPrefix + TreeModel.TREE_MODEL);
                        writer.writeCloseTag(CoalescentLikelihoodParser.POPULATION_TREE);
                        writer.writeCloseTag(CoalescentLikelihoodParser.COALESCENT_LIKELIHOOD);

                        break;

                    case LOGISTIC:



                        writer.writeComment("A working prior for the coalescent.");
                        writer.writeOpenTag(
                                CoalescentLikelihoodParser.COALESCENT_LIKELIHOOD,
                                new Attribute[]{
                                        new Attribute.Default<String>(XMLParser.ID, modelPrefix + "coalescentReference")
                                }
                        );
                        writer.writeOpenTag(CoalescentLikelihoodParser.MODEL);
                        writer.writeIDref(LogisticGrowthModelParser.LOGISTIC_GROWTH_MODEL, beautiOptions.getPartitionTreePriors().get(0).getPrefix() + "constantReference");
                        writer.writeCloseTag(CoalescentLikelihoodParser.MODEL);
                        writer.writeOpenTag(CoalescentLikelihoodParser.POPULATION_TREE);
                        writer.writeIDref(TreeModel.TREE_MODEL, modelPrefix + TreeModel.TREE_MODEL);
                        writer.writeCloseTag(CoalescentLikelihoodParser.POPULATION_TREE);
                        writer.writeCloseTag(CoalescentLikelihoodParser.COALESCENT_LIKELIHOOD);

                        break;

                    case EXPANSION:



                        writer.writeComment("A working prior for the coalescent.");
                        writer.writeOpenTag(
                                CoalescentLikelihoodParser.COALESCENT_LIKELIHOOD,
                                new Attribute[]{
                                        new Attribute.Default<String>(XMLParser.ID, modelPrefix + "coalescentReference")
                                }
                        );
                        writer.writeOpenTag(CoalescentLikelihoodParser.MODEL);
                        writer.writeIDref(ExpansionModelParser.EXPANSION_MODEL, beautiOptions.getPartitionTreePriors().get(0).getPrefix() + "constantReference");
                        writer.writeCloseTag(CoalescentLikelihoodParser.MODEL);
                        writer.writeOpenTag(CoalescentLikelihoodParser.POPULATION_TREE);
                        writer.writeIDref(TreeModel.TREE_MODEL, modelPrefix + TreeModel.TREE_MODEL);
                        writer.writeCloseTag(CoalescentLikelihoodParser.POPULATION_TREE);
                        writer.writeCloseTag(CoalescentLikelihoodParser.COALESCENT_LIKELIHOOD);

                        break;

                    default:

                        //TODO: show menu that explains mismatch between prior and working prior

                }

                //TODO: if not a simple coalescent model, switch to product of exponentials

            }

            writer.writeComment("Define marginal likelihood estimator (GSS) settings");

            List<Attribute> attributes = new ArrayList<Attribute>();
            attributes.add(new Attribute.Default<Integer>("chainLength", options.mleChainLength));
            attributes.add(new Attribute.Default<Integer>("pathSteps", options.pathSteps));
            attributes.add(new Attribute.Default<String>("pathScheme", options.pathScheme));
            if (!options.pathScheme.equals("linear")) {
                attributes.add(new Attribute.Default<Double>("alpha", options.schemeParameter));
            }

            writer.writeOpenTag(MarginalLikelihoodEstimator.MARGINAL_LIKELIHOOD_ESTIMATOR, attributes);

            writer.writeOpenTag("samplers");
            writer.writeIDref("mcmc", "mcmc");
            writer.writeCloseTag("samplers");

            attributes = new ArrayList<Attribute>();
            attributes.add(new Attribute.Default<String>(XMLParser.ID, "pathLikelihood"));
            writer.writeOpenTag(PathLikelihood.PATH_LIKELIHOOD, attributes);
            writer.writeOpenTag(PathLikelihood.SOURCE);
            writer.writeIDref(CompoundLikelihoodParser.POSTERIOR, CompoundLikelihoodParser.POSTERIOR);
            writer.writeCloseTag(PathLikelihood.SOURCE);
            writer.writeOpenTag(PathLikelihood.DESTINATION);
            writer.writeOpenTag(CompoundLikelihoodParser.WORKING_PRIOR);

            ArrayList<Parameter> parameters = beautiOptions.selectParameters();

            for (Parameter param : parameters) {
                if (DEBUG) {
                    System.err.println(param.toString() + "   " + param.priorType.toString());
                }
                //should leave out those parameters set by the coalescent
                if (param.priorType != PriorType.NONE_TREE_PRIOR) {
                    //TODO: frequencies is multidimensional, is that automatically dealt with?
                    writer.writeOpenTag(WorkingPriorParsers.NORMAL_REFERENCE_PRIOR,
                            new Attribute[]{
                                    new Attribute.Default<String>("fileName", beautiOptions.logFileName),
                                    new Attribute.Default<String>("parameterColumn", param.getName()),
                                    new Attribute.Default<String>("burnin", "" + beautiOptions.chainLength*0.10)
                            });
                    writeParameterIdref(writer, param);
                    writer.writeCloseTag(WorkingPriorParsers.NORMAL_REFERENCE_PRIOR);
                }
            }

            if (options.choiceTreeWorkingPrior.equals("Product of exponential distributions")) {
                writer.writeIDref("productOfExponentialsPosteriorMeansLoess", "exponentials");
            } else {
                writer.writeIDref(CoalescentLikelihoodParser.COALESCENT_LIKELIHOOD, "coalescentReference");
            }

            writer.writeCloseTag(CompoundLikelihoodParser.WORKING_PRIOR);
            writer.writeCloseTag(PathLikelihood.DESTINATION);
            writer.writeCloseTag(PathLikelihood.PATH_LIKELIHOOD);

            attributes = new ArrayList<Attribute>();
            attributes.add(new Attribute.Default<String>(XMLParser.ID, "MLELog"));
            attributes.add(new Attribute.Default<Integer>("logEvery", options.mleLogEvery));
            attributes.add(new Attribute.Default<String>("fileName", options.mleFileName));
            writer.writeOpenTag("log", attributes);
            writer.writeIDref("pathLikelihood", "pathLikelihood");
            writer.writeCloseTag("log");

            writer.writeCloseTag(MarginalLikelihoodEstimator.MARGINAL_LIKELIHOOD_ESTIMATOR);

            writer.writeComment("Generalized stepping-stone sampling estimator from collected samples");
            attributes = new ArrayList<Attribute>();
            attributes.add(new Attribute.Default<String>("fileName", options.mleFileName));
            writer.writeOpenTag(GeneralizedSteppingStoneSamplingAnalysis.GENERALIZED_STEPPING_STONE_SAMPLING_ANALYSIS, attributes);
            writer.writeTag("sourceColumn", new Attribute.Default<String>("name", "pathLikelihood.source"), true);
            writer.writeTag("destinationColumn", new Attribute.Default<String>("name", "pathLikelihood.destination"), true);
            writer.writeTag("thetaColumn", new Attribute.Default<String>("name", "pathLikelihood.theta"), true);
            writer.writeCloseTag(GeneralizedSteppingStoneSamplingAnalysis.GENERALIZED_STEPPING_STONE_SAMPLING_ANALYSIS);

        }

    }

    private void writeParameterIdref(XMLWriter writer, Parameter parameter) {
        if (parameter.isStatistic) {
            writer.writeIDref("statistic", parameter.getName());
        } else {
            writer.writeIDref(ParameterParser.PARAMETER, parameter.getName());
        }
    }

}
