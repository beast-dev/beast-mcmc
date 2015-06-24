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
import dr.app.beauti.options.*;
import dr.app.beauti.types.*;
import dr.app.beauti.util.XMLWriter;
import dr.evolution.datatype.DataType;
import dr.evolution.util.Units;
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
 * @author Guy Baele
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
            attributes.add(new Attribute.Default<Integer>(MarginalLikelihoodEstimator.CHAIN_LENGTH, options.mleChainLength));
            attributes.add(new Attribute.Default<Integer>(MarginalLikelihoodEstimator.PATH_STEPS, options.pathSteps));
            attributes.add(new Attribute.Default<String>(MarginalLikelihoodEstimator.PATH_SCHEME, options.pathScheme));
            if (!options.pathScheme.equals(MarginalLikelihoodEstimator.LINEAR)) {
                attributes.add(new Attribute.Default<Double>(MarginalLikelihoodEstimator.ALPHA, options.schemeParameter));
            }
            if (options.printOperatorAnalysis) {
                attributes.add(new Attribute.Default<Boolean>(MarginalLikelihoodEstimator.PRINT_OPERATOR_ANALYSIS, true));
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

                        writer.writeComment("A working prior for the exponential growth model.");
                        writer.writeOpenTag(
                                ExponentialGrowthModelParser.EXPONENTIAL_GROWTH_MODEL,
                                new Attribute[]{
                                        new Attribute.Default<String>(XMLParser.ID, modelPrefix + "exponentialReference"),
                                        new Attribute.Default<String>("units", Units.Utils.getDefaultUnitName(beautiOptions.units))
                                }
                        );

                        writer.writeOpenTag(ExponentialGrowthModelParser.POPULATION_SIZE);
                        writeParameter("exponentialReference.popSize", "exponential.popSize", beautiOptions.logFileName, (int) (options.mleChainLength * 0.10), writer);
                        writer.writeCloseTag(ExponentialGrowthModelParser.POPULATION_SIZE);
                        writer.writeOpenTag(ExponentialGrowthModelParser.GROWTH_RATE);
                        writeParameter("exponentialReference.growthRate", "exponential.growthRate", beautiOptions.logFileName, (int) (options.mleChainLength * 0.10), writer);
                        writer.writeCloseTag(ExponentialGrowthModelParser.GROWTH_RATE);
                        writer.writeCloseTag(ExponentialGrowthModelParser.EXPONENTIAL_GROWTH_MODEL);

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

                        writer.writeComment("A working prior for the logistic growth model.");
                        writer.writeOpenTag(
                                LogisticGrowthModelParser.LOGISTIC_GROWTH_MODEL,
                                new Attribute[]{
                                        new Attribute.Default<String>(XMLParser.ID, modelPrefix + "logisticReference"),
                                        new Attribute.Default<String>("units", Units.Utils.getDefaultUnitName(beautiOptions.units))
                                }
                        );

                        writer.writeOpenTag(LogisticGrowthModelParser.POPULATION_SIZE);
                        writeParameter("logisticReference.popSize", "logistic.popSize", beautiOptions.logFileName, (int) (options.mleChainLength * 0.10), writer);
                        writer.writeCloseTag(LogisticGrowthModelParser.POPULATION_SIZE);
                        writer.writeOpenTag(LogisticGrowthModelParser.GROWTH_RATE);
                        writeParameter("logisticReference.growthRate", "logistic.growthRate", beautiOptions.logFileName, (int) (options.mleChainLength * 0.10), writer);
                        writer.writeCloseTag(LogisticGrowthModelParser.GROWTH_RATE);
                        writer.writeOpenTag(LogisticGrowthModelParser.TIME_50);
                        writeParameter("logisticReference.t50", "logistic.t50", beautiOptions.logFileName, (int) (options.mleChainLength * 0.10), writer);
                        writer.writeCloseTag(LogisticGrowthModelParser.TIME_50);
                        writer.writeCloseTag(LogisticGrowthModelParser.LOGISTIC_GROWTH_MODEL);

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

                        writer.writeComment("A working prior for the expansion growth model.");
                        writer.writeOpenTag(
                                ExpansionModelParser.EXPANSION_MODEL,
                                new Attribute[]{
                                        new Attribute.Default<String>(XMLParser.ID, modelPrefix + "expansionReference"),
                                        new Attribute.Default<String>("units", Units.Utils.getDefaultUnitName(beautiOptions.units))
                                }
                        );

                        writer.writeOpenTag(ExpansionModelParser.POPULATION_SIZE);
                        writeParameter("expansionReference.popSize", "expansion.popSize", beautiOptions.logFileName, (int) (options.mleChainLength * 0.10), writer);
                        writer.writeCloseTag(ExpansionModelParser.POPULATION_SIZE);
                        writer.writeOpenTag(ExpansionModelParser.GROWTH_RATE);
                        writeParameter("expansionReference.growthRate", "expansion.growthRate", beautiOptions.logFileName, (int) (options.mleChainLength * 0.10), writer);
                        writer.writeCloseTag(ExpansionModelParser.GROWTH_RATE);
                        writer.writeOpenTag(ExpansionModelParser.ANCESTRAL_POPULATION_PROPORTION);
                        writeParameter("expansionReference.ancestralProportion", "expansion.ancestralProportion", beautiOptions.logFileName, (int) (options.mleChainLength * 0.10), writer);
                        writer.writeCloseTag(ExpansionModelParser.ANCESTRAL_POPULATION_PROPORTION);
                        writer.writeCloseTag(ExpansionModelParser.EXPANSION_MODEL);

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

                        //Do not switch to product of exponentials as the coalescentEventsStatistic has not been logged
                        //TODO: show menu that explains mismatch between prior and working prior?
                        //TODO: but show it when the MCM option is wrongfully being selected, don't do anything here

                }

            }

            writer.writeComment("Define marginal likelihood estimator (GSS) settings");

            List<Attribute> attributes = new ArrayList<Attribute>();
            attributes.add(new Attribute.Default<Integer>(MarginalLikelihoodEstimator.CHAIN_LENGTH, options.mleChainLength));
            attributes.add(new Attribute.Default<Integer>(MarginalLikelihoodEstimator.PATH_STEPS, options.pathSteps));
            attributes.add(new Attribute.Default<String>(MarginalLikelihoodEstimator.PATH_SCHEME, options.pathScheme));
            if (!options.pathScheme.equals(MarginalLikelihoodEstimator.LINEAR)) {
                attributes.add(new Attribute.Default<Double>(MarginalLikelihoodEstimator.ALPHA, options.schemeParameter));
            }
            if (options.printOperatorAnalysis) {
                attributes.add(new Attribute.Default<Boolean>(MarginalLikelihoodEstimator.PRINT_OPERATOR_ANALYSIS, true));
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

            //Start with providing working priors for the substitution model(s)
            for (PartitionSubstitutionModel model : beautiOptions.getPartitionSubstitutionModels()) {
                int codonPartitionCount = model.getCodonPartitionCount();

                switch (model.getDataType().getType()) {
                    case DataType.NUCLEOTIDES:

                        switch (model.getNucSubstitutionModel()) {

                            case HKY:
                                if (codonPartitionCount > 1 && model.isUnlinkedSubstitutionModel()) {
                                    for (int i = 1; i <= codonPartitionCount; i++) {
                                        writer.writeOpenTag(WorkingPriorParsers.LOG_TRANSFORMED_NORMAL_REFERENCE_PRIOR,
                                                new Attribute[]{
                                                        new Attribute.Default<String>("fileName", beautiOptions.logFileName),
                                                        new Attribute.Default<String>("parameterColumn", model.getPrefix(i) + "kappa"),
                                                        new Attribute.Default<String>("burnin", "" + beautiOptions.chainLength*0.10)
                                                });
                                        writer.writeIDref(ParameterParser.PARAMETER, model.getPrefix(i) + "kappa");
                                        writer.writeCloseTag(WorkingPriorParsers.LOG_TRANSFORMED_NORMAL_REFERENCE_PRIOR);
                                    }
                                } else {
                                    writer.writeOpenTag(WorkingPriorParsers.LOG_TRANSFORMED_NORMAL_REFERENCE_PRIOR,
                                            new Attribute[]{
                                                    new Attribute.Default<String>("fileName", beautiOptions.logFileName),
                                                    new Attribute.Default<String>("parameterColumn", model.getPrefix() + "kappa"),
                                                    new Attribute.Default<String>("burnin", "" + beautiOptions.chainLength*0.10)
                                            });
                                    writer.writeIDref(ParameterParser.PARAMETER, model.getPrefix() + "kappa");
                                    writer.writeCloseTag(WorkingPriorParsers.LOG_TRANSFORMED_NORMAL_REFERENCE_PRIOR);
                                }
                                break;

                            case TN93:
                                if (codonPartitionCount > 1 && model.isUnlinkedSubstitutionModel()) {
                                    for (int i = 1; i <= codonPartitionCount; i++) {
                                        writer.writeOpenTag(WorkingPriorParsers.LOG_TRANSFORMED_NORMAL_REFERENCE_PRIOR,
                                                new Attribute[]{
                                                        new Attribute.Default<String>("fileName", beautiOptions.logFileName),
                                                        new Attribute.Default<String>("parameterColumn", model.getPrefix(i) + "kappa1"),
                                                        new Attribute.Default<String>("burnin", "" + beautiOptions.chainLength*0.10)
                                                });
                                        writer.writeIDref(ParameterParser.PARAMETER, model.getPrefix(i) + "kappa1");
                                        writer.writeCloseTag(WorkingPriorParsers.LOG_TRANSFORMED_NORMAL_REFERENCE_PRIOR);
                                        writer.writeOpenTag(WorkingPriorParsers.LOG_TRANSFORMED_NORMAL_REFERENCE_PRIOR,
                                                new Attribute[]{
                                                        new Attribute.Default<String>("fileName", beautiOptions.logFileName),
                                                        new Attribute.Default<String>("parameterColumn", model.getPrefix(i) + "kappa2"),
                                                        new Attribute.Default<String>("burnin", "" + beautiOptions.chainLength*0.10)
                                                });
                                        writer.writeIDref(ParameterParser.PARAMETER, model.getPrefix(i) + "kappa2");
                                        writer.writeCloseTag(WorkingPriorParsers.LOG_TRANSFORMED_NORMAL_REFERENCE_PRIOR);
                                    }
                                } else {
                                    writer.writeOpenTag(WorkingPriorParsers.LOG_TRANSFORMED_NORMAL_REFERENCE_PRIOR,
                                            new Attribute[]{
                                                    new Attribute.Default<String>("fileName", beautiOptions.logFileName),
                                                    new Attribute.Default<String>("parameterColumn", model.getPrefix() + "kappa1"),
                                                    new Attribute.Default<String>("burnin", "" + beautiOptions.chainLength*0.10)
                                            });
                                    writer.writeIDref(ParameterParser.PARAMETER, model.getPrefix() + "kappa1");
                                    writer.writeCloseTag(WorkingPriorParsers.LOG_TRANSFORMED_NORMAL_REFERENCE_PRIOR);
                                    writer.writeOpenTag(WorkingPriorParsers.LOG_TRANSFORMED_NORMAL_REFERENCE_PRIOR,
                                            new Attribute[]{
                                                    new Attribute.Default<String>("fileName", beautiOptions.logFileName),
                                                    new Attribute.Default<String>("parameterColumn", model.getPrefix() + "kappa2"),
                                                    new Attribute.Default<String>("burnin", "" + beautiOptions.chainLength*0.10)
                                            });
                                    writer.writeIDref(ParameterParser.PARAMETER, model.getPrefix() + "kappa2");
                                    writer.writeCloseTag(WorkingPriorParsers.LOG_TRANSFORMED_NORMAL_REFERENCE_PRIOR);
                                }
                                break;

                            case GTR:
                                if (codonPartitionCount > 1 && model.isUnlinkedSubstitutionModel()) {
                                    for (int i = 1; i <= codonPartitionCount; i++) {
                                        for (String rateName : PartitionSubstitutionModel.GTR_RATE_NAMES) {
                                            writer.writeOpenTag(WorkingPriorParsers.LOG_TRANSFORMED_NORMAL_REFERENCE_PRIOR,
                                                    new Attribute[]{
                                                            new Attribute.Default<String>("fileName", beautiOptions.logFileName),
                                                            new Attribute.Default<String>("parameterColumn", model.getPrefix(i) + rateName),
                                                            new Attribute.Default<String>("burnin", "" + beautiOptions.chainLength * 0.10)
                                                    });
                                            writer.writeIDref(ParameterParser.PARAMETER, model.getPrefix(i) + rateName);
                                            writer.writeCloseTag(WorkingPriorParsers.LOG_TRANSFORMED_NORMAL_REFERENCE_PRIOR);
                                        }
                                    }
                                } else {
                                    for (String rateName : PartitionSubstitutionModel.GTR_RATE_NAMES) {
                                        writer.writeOpenTag(WorkingPriorParsers.LOG_TRANSFORMED_NORMAL_REFERENCE_PRIOR,
                                                new Attribute[]{
                                                        new Attribute.Default<String>("fileName", beautiOptions.logFileName),
                                                        new Attribute.Default<String>("parameterColumn", model.getPrefix() + rateName),
                                                        new Attribute.Default<String>("burnin", "" + beautiOptions.chainLength * 0.10)
                                                });
                                        writer.writeIDref(ParameterParser.PARAMETER, model.getPrefix() + rateName);
                                        writer.writeCloseTag(WorkingPriorParsers.LOG_TRANSFORMED_NORMAL_REFERENCE_PRIOR);
                                    }
                                }
                                break;

                        }

                        if (model.getFrequencyPolicy() == FrequencyPolicyType.ESTIMATED) {
                            if (codonPartitionCount > 1 && model.isUnlinkedSubstitutionModel() && model.isUnlinkedFrequencyModel()) {
                                for (int i = 1; i <= codonPartitionCount; i++) {
                                    writer.writeOpenTag(WorkingPriorParsers.LOGIT_TRANSFORMED_NORMAL_REFERENCE_PRIOR,
                                            new Attribute[]{
                                                    new Attribute.Default<String>("fileName", beautiOptions.logFileName),
                                                    new Attribute.Default<String>("parameterColumn", model.getPrefix(i) + "frequencies"),
                                                    new Attribute.Default<Integer>("dimension", 4),
                                                    new Attribute.Default<String>("burnin", "" + beautiOptions.chainLength*0.10)
                                            });
                                    writer.writeIDref(ParameterParser.PARAMETER, model.getPrefix(i) + "frequencies");
                                    writer.writeCloseTag(WorkingPriorParsers.LOGIT_TRANSFORMED_NORMAL_REFERENCE_PRIOR);
                                }
                            } else {
                                writer.writeOpenTag(WorkingPriorParsers.LOGIT_TRANSFORMED_NORMAL_REFERENCE_PRIOR,
                                        new Attribute[]{
                                                new Attribute.Default<String>("fileName", beautiOptions.logFileName),
                                                new Attribute.Default<String>("parameterColumn", model.getPrefix() + "frequencies"),
                                                new Attribute.Default<Integer>("dimension", 4),
                                                new Attribute.Default<String>("burnin", "" + beautiOptions.chainLength*0.10)
                                        });
                                writer.writeIDref(ParameterParser.PARAMETER, model.getPrefix() + "frequencies");
                                writer.writeCloseTag(WorkingPriorParsers.LOGIT_TRANSFORMED_NORMAL_REFERENCE_PRIOR);
                            }
                        }

                        break;//NUCLEOTIDES

                    default:
                        throw new IllegalArgumentException("Unknown data type");
                }

                if (model.isGammaHetero()) {
                    if (codonPartitionCount > 1 && model.isUnlinkedHeterogeneityModel()) {
                        for (int i = 1; i <= codonPartitionCount; i++) {
                            writer.writeOpenTag(WorkingPriorParsers.LOG_TRANSFORMED_NORMAL_REFERENCE_PRIOR,
                                    new Attribute[]{
                                            new Attribute.Default<String>("fileName", beautiOptions.logFileName),
                                            new Attribute.Default<String>("parameterColumn", model.getPrefix(i) + "alpha"),
                                            new Attribute.Default<String>("burnin", "" + beautiOptions.chainLength*0.10)
                                    });
                            writer.writeIDref(ParameterParser.PARAMETER, model.getPrefix(i) + "alpha");
                            writer.writeCloseTag(WorkingPriorParsers.LOG_TRANSFORMED_NORMAL_REFERENCE_PRIOR);
                        }
                    } else {
                        writer.writeOpenTag(WorkingPriorParsers.LOG_TRANSFORMED_NORMAL_REFERENCE_PRIOR,
                                new Attribute[]{
                                        new Attribute.Default<String>("fileName", beautiOptions.logFileName),
                                        new Attribute.Default<String>("parameterColumn", model.getPrefix() + "alpha"),
                                        new Attribute.Default<String>("burnin", "" + beautiOptions.chainLength*0.10)
                                });
                        writer.writeIDref(ParameterParser.PARAMETER, model.getPrefix() + "alpha");
                        writer.writeCloseTag(WorkingPriorParsers.LOG_TRANSFORMED_NORMAL_REFERENCE_PRIOR);
                    }
                }

                if (model.isInvarHetero()) {
                    if (codonPartitionCount > 1 && model.isUnlinkedHeterogeneityModel()) {
                        for (int i = 1; i <= codonPartitionCount; i++) {
                            writer.writeOpenTag(WorkingPriorParsers.LOGIT_TRANSFORMED_NORMAL_REFERENCE_PRIOR,
                                    new Attribute[]{
                                            new Attribute.Default<String>("fileName", beautiOptions.logFileName),
                                            new Attribute.Default<String>("parameterColumn", model.getPrefix(i) + "pInv"),
                                            new Attribute.Default<String>("burnin", "" + beautiOptions.chainLength*0.10)
                                    });
                            writer.writeIDref(ParameterParser.PARAMETER, model.getPrefix(i) + "pInv");
                            writer.writeCloseTag(WorkingPriorParsers.LOGIT_TRANSFORMED_NORMAL_REFERENCE_PRIOR);
                        }
                    } else {
                        writer.writeOpenTag(WorkingPriorParsers.LOGIT_TRANSFORMED_NORMAL_REFERENCE_PRIOR,
                                new Attribute[]{
                                        new Attribute.Default<String>("fileName", beautiOptions.logFileName),
                                        new Attribute.Default<String>("parameterColumn", model.getPrefix() + "pInv"),
                                        new Attribute.Default<String>("burnin", "" + beautiOptions.chainLength*0.10)
                                });
                        writer.writeIDref(ParameterParser.PARAMETER, model.getPrefix() + "pInv");
                        writer.writeCloseTag(WorkingPriorParsers.LOGIT_TRANSFORMED_NORMAL_REFERENCE_PRIOR);
                    }
                }
            }

            //Continue with providing working priors for the clock model(s)
            for (PartitionClockModel model : beautiOptions.getPartitionClockModels()) {
                switch(model.getClockType()) {
                    case STRICT_CLOCK:
                        writer.writeOpenTag(WorkingPriorParsers.LOG_TRANSFORMED_NORMAL_REFERENCE_PRIOR,
                                new Attribute[]{
                                        new Attribute.Default<String>("fileName", beautiOptions.logFileName),
                                        new Attribute.Default<String>("parameterColumn", model.getPrefix() + "clock.rate"),
                                        new Attribute.Default<String>("burnin", "" + beautiOptions.chainLength*0.10)
                                });
                        writer.writeIDref(ParameterParser.PARAMETER, model.getPrefix() + "clock.rate");
                        writer.writeCloseTag(WorkingPriorParsers.LOG_TRANSFORMED_NORMAL_REFERENCE_PRIOR);
                        break;

                    case UNCORRELATED:
                        switch (model.getClockDistributionType()) {
                            case LOGNORMAL:
                                if (model.getClockRateParam().isMeanInRealSpace()) {
                                    writer.writeOpenTag(WorkingPriorParsers.LOG_TRANSFORMED_NORMAL_REFERENCE_PRIOR,
                                            new Attribute[]{
                                                    new Attribute.Default<String>("fileName", beautiOptions.logFileName),
                                                    new Attribute.Default<String>("parameterColumn", model.getPrefix() + ClockType.UCLD_MEAN),
                                                    new Attribute.Default<String>("burnin", "" + beautiOptions.chainLength*0.10)
                                            });
                                    writer.writeIDref(ParameterParser.PARAMETER, model.getPrefix() + ClockType.UCLD_MEAN);
                                    writer.writeCloseTag(WorkingPriorParsers.LOG_TRANSFORMED_NORMAL_REFERENCE_PRIOR);
                                    writer.writeOpenTag(WorkingPriorParsers.LOG_TRANSFORMED_NORMAL_REFERENCE_PRIOR,
                                            new Attribute[]{
                                                    new Attribute.Default<String>("fileName", beautiOptions.logFileName),
                                                    new Attribute.Default<String>("parameterColumn", model.getPrefix() + ClockType.UCLD_STDEV),
                                                    new Attribute.Default<String>("burnin", "" + beautiOptions.chainLength*0.10)
                                            });
                                    writer.writeIDref(ParameterParser.PARAMETER, model.getPrefix() + ClockType.UCLD_STDEV);
                                    writer.writeCloseTag(WorkingPriorParsers.LOG_TRANSFORMED_NORMAL_REFERENCE_PRIOR);
                                } else {
                                    writer.writeOpenTag(WorkingPriorParsers.NORMAL_REFERENCE_PRIOR,
                                            new Attribute[]{
                                                    new Attribute.Default<String>("fileName", beautiOptions.logFileName),
                                                    new Attribute.Default<String>("parameterColumn", model.getPrefix() + ClockType.UCLD_MEAN),
                                                    new Attribute.Default<String>("burnin", "" + beautiOptions.chainLength*0.10)
                                            });
                                    writer.writeIDref(ParameterParser.PARAMETER, model.getPrefix() + ClockType.UCLD_MEAN);
                                    writer.writeCloseTag(WorkingPriorParsers.NORMAL_REFERENCE_PRIOR);
                                    writer.writeOpenTag(WorkingPriorParsers.LOG_TRANSFORMED_NORMAL_REFERENCE_PRIOR,
                                            new Attribute[]{
                                                    new Attribute.Default<String>("fileName", beautiOptions.logFileName),
                                                    new Attribute.Default<String>("parameterColumn", model.getPrefix() + ClockType.UCLD_STDEV),
                                                    new Attribute.Default<String>("burnin", "" + beautiOptions.chainLength*0.10)
                                            });
                                    writer.writeIDref(ParameterParser.PARAMETER, model.getPrefix() + ClockType.UCLD_STDEV);
                                    writer.writeCloseTag(WorkingPriorParsers.LOG_TRANSFORMED_NORMAL_REFERENCE_PRIOR);
                                }
                                break;

                            case EXPONENTIAL:
                                writer.writeOpenTag(WorkingPriorParsers.LOG_TRANSFORMED_NORMAL_REFERENCE_PRIOR,
                                        new Attribute[]{
                                                new Attribute.Default<String>("fileName", beautiOptions.logFileName),
                                                new Attribute.Default<String>("parameterColumn", model.getPrefix() + ClockType.UCED_MEAN),
                                                new Attribute.Default<String>("burnin", "" + beautiOptions.chainLength*0.10)
                                        });
                                writer.writeIDref(ParameterParser.PARAMETER, model.getPrefix() + ClockType.UCED_MEAN);
                                writer.writeCloseTag(WorkingPriorParsers.LOG_TRANSFORMED_NORMAL_REFERENCE_PRIOR);
                                break;
                        }
                        break;

                    default:
                        throw new IllegalArgumentException("Unknown clock model");
                }
            }

            //Provide working priors for the coalescent model(s)
            for (PartitionTreePrior model : beautiOptions.getPartitionTreePriors()) {
                TreePriorType nodeHeightPrior = model.getNodeHeightPrior();
                TreePriorParameterizationType parameterization = model.getParameterization();

                switch (nodeHeightPrior) {
                    case CONSTANT:
                        writer.writeOpenTag(WorkingPriorParsers.LOG_TRANSFORMED_NORMAL_REFERENCE_PRIOR,
                                new Attribute[]{
                                        new Attribute.Default<String>("fileName", beautiOptions.logFileName),
                                        new Attribute.Default<String>("parameterColumn", "constant.popSize"),
                                        new Attribute.Default<String>("burnin", "" + beautiOptions.chainLength*0.10)
                                });
                        writer.writeIDref(ParameterParser.PARAMETER, "constant.popSize");
                        writer.writeCloseTag(WorkingPriorParsers.LOG_TRANSFORMED_NORMAL_REFERENCE_PRIOR);
                        break;

                    case EXPONENTIAL:
                        writer.writeOpenTag(WorkingPriorParsers.LOG_TRANSFORMED_NORMAL_REFERENCE_PRIOR,
                                new Attribute[]{
                                        new Attribute.Default<String>("fileName", beautiOptions.logFileName),
                                        new Attribute.Default<String>("parameterColumn", "exponential.popSize"),
                                        new Attribute.Default<String>("burnin", "" + beautiOptions.chainLength*0.10)
                                });
                        writer.writeIDref(ParameterParser.PARAMETER, "exponential.popSize");
                        writer.writeCloseTag(WorkingPriorParsers.LOG_TRANSFORMED_NORMAL_REFERENCE_PRIOR);

                        if (parameterization == TreePriorParameterizationType.GROWTH_RATE) {
                            writer.writeOpenTag(WorkingPriorParsers.NORMAL_REFERENCE_PRIOR,
                                    new Attribute[]{
                                            new Attribute.Default<String>("fileName", beautiOptions.logFileName),
                                            new Attribute.Default<String>("parameterColumn", "exponential.growthRate"),
                                            new Attribute.Default<String>("burnin", "" + beautiOptions.chainLength * 0.10)
                                    });
                            writer.writeIDref(ParameterParser.PARAMETER, "exponential.growthRate");
                            writer.writeCloseTag(WorkingPriorParsers.NORMAL_REFERENCE_PRIOR);
                        } else {
                            writer.writeOpenTag(WorkingPriorParsers.LOG_TRANSFORMED_NORMAL_REFERENCE_PRIOR,
                                    new Attribute[]{
                                            new Attribute.Default<String>("fileName", beautiOptions.logFileName),
                                            new Attribute.Default<String>("parameterColumn", "exponential.doublingTime"),
                                            new Attribute.Default<String>("burnin", "" + beautiOptions.chainLength * 0.10)
                                    });
                            writer.writeIDref(ParameterParser.PARAMETER, "exponential.doublingTime");
                            writer.writeCloseTag(WorkingPriorParsers.LOG_TRANSFORMED_NORMAL_REFERENCE_PRIOR);
                        }

                        break;

                    case LOGISTIC:
                        writer.writeOpenTag(WorkingPriorParsers.LOG_TRANSFORMED_NORMAL_REFERENCE_PRIOR,
                                new Attribute[]{
                                        new Attribute.Default<String>("fileName", beautiOptions.logFileName),
                                        new Attribute.Default<String>("parameterColumn", "logistic.popSize"),
                                        new Attribute.Default<String>("burnin", "" + beautiOptions.chainLength*0.10)
                                });
                        writer.writeIDref(ParameterParser.PARAMETER, "logistic.popSize");
                        writer.writeCloseTag(WorkingPriorParsers.LOG_TRANSFORMED_NORMAL_REFERENCE_PRIOR);

                        if (parameterization == TreePriorParameterizationType.GROWTH_RATE) {
                            writer.writeOpenTag(WorkingPriorParsers.NORMAL_REFERENCE_PRIOR,
                                    new Attribute[]{
                                            new Attribute.Default<String>("fileName", beautiOptions.logFileName),
                                            new Attribute.Default<String>("parameterColumn", "logistic.growthRate"),
                                            new Attribute.Default<String>("burnin", "" + beautiOptions.chainLength * 0.10)
                                    });
                            writer.writeIDref(ParameterParser.PARAMETER, "logistic.growthRate");
                            writer.writeCloseTag(WorkingPriorParsers.NORMAL_REFERENCE_PRIOR);
                        } else {
                            writer.writeOpenTag(WorkingPriorParsers.LOG_TRANSFORMED_NORMAL_REFERENCE_PRIOR,
                                    new Attribute[]{
                                            new Attribute.Default<String>("fileName", beautiOptions.logFileName),
                                            new Attribute.Default<String>("parameterColumn", "logistic.doublingTime"),
                                            new Attribute.Default<String>("burnin", "" + beautiOptions.chainLength * 0.10)
                                    });
                            writer.writeIDref(ParameterParser.PARAMETER, "logistic.doublingTime");
                            writer.writeCloseTag(WorkingPriorParsers.LOG_TRANSFORMED_NORMAL_REFERENCE_PRIOR);
                        }

                        writer.writeOpenTag(WorkingPriorParsers.LOG_TRANSFORMED_NORMAL_REFERENCE_PRIOR,
                                new Attribute[]{
                                        new Attribute.Default<String>("fileName", beautiOptions.logFileName),
                                        new Attribute.Default<String>("parameterColumn", "logistic.t50"),
                                        new Attribute.Default<String>("burnin", "" + beautiOptions.chainLength*0.10)
                                });
                        writer.writeIDref(ParameterParser.PARAMETER, "logistic.t50");
                        writer.writeCloseTag(WorkingPriorParsers.LOG_TRANSFORMED_NORMAL_REFERENCE_PRIOR);

                        break;

                    case EXPANSION:
                        writer.writeOpenTag(WorkingPriorParsers.LOG_TRANSFORMED_NORMAL_REFERENCE_PRIOR,
                                new Attribute[]{
                                        new Attribute.Default<String>("fileName", beautiOptions.logFileName),
                                        new Attribute.Default<String>("parameterColumn", "expansion.popSize"),
                                        new Attribute.Default<String>("burnin", "" + beautiOptions.chainLength*0.10)
                                });
                        writer.writeIDref(ParameterParser.PARAMETER, "expansion.popSize");
                        writer.writeCloseTag(WorkingPriorParsers.LOG_TRANSFORMED_NORMAL_REFERENCE_PRIOR);

                        if (parameterization == TreePriorParameterizationType.GROWTH_RATE) {
                            writer.writeOpenTag(WorkingPriorParsers.NORMAL_REFERENCE_PRIOR,
                                    new Attribute[]{
                                            new Attribute.Default<String>("fileName", beautiOptions.logFileName),
                                            new Attribute.Default<String>("parameterColumn", "expansion.growthRate"),
                                            new Attribute.Default<String>("burnin", "" + beautiOptions.chainLength * 0.10)
                                    });
                            writer.writeIDref(ParameterParser.PARAMETER, "expansion.growthRate");
                            writer.writeCloseTag(WorkingPriorParsers.NORMAL_REFERENCE_PRIOR);
                        } else {
                            writer.writeOpenTag(WorkingPriorParsers.LOG_TRANSFORMED_NORMAL_REFERENCE_PRIOR,
                                    new Attribute[]{
                                            new Attribute.Default<String>("fileName", beautiOptions.logFileName),
                                            new Attribute.Default<String>("parameterColumn", "expansion.doublingTime"),
                                            new Attribute.Default<String>("burnin", "" + beautiOptions.chainLength * 0.10)
                                    });
                            writer.writeIDref(ParameterParser.PARAMETER, "expansion.doublingTime");
                            writer.writeCloseTag(WorkingPriorParsers.LOG_TRANSFORMED_NORMAL_REFERENCE_PRIOR);
                        }

                        writer.writeOpenTag(WorkingPriorParsers.LOGIT_TRANSFORMED_NORMAL_REFERENCE_PRIOR,
                                new Attribute[]{
                                        new Attribute.Default<String>("fileName", beautiOptions.logFileName),
                                        new Attribute.Default<String>("parameterColumn", "expansion.ancestralProportion"),
                                        new Attribute.Default<String>("burnin", "" + beautiOptions.chainLength*0.10)
                                });
                        writer.writeIDref(ParameterParser.PARAMETER, "expansion.ancestralProportion");
                        writer.writeCloseTag(WorkingPriorParsers.LOGIT_TRANSFORMED_NORMAL_REFERENCE_PRIOR);

                        break;
                }
            }

            //TODO: take care of anything else I missed



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
