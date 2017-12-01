/*
 * MarginalLikelihoodEstimationGenerator.java
 *
 * Copyright (c) 2002-2017 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

import dr.app.beauti.BeautiFrame;
import dr.app.beauti.generator.BaseComponentGenerator;
import dr.app.beauti.generator.TreePriorGenerator;
import dr.app.beauti.options.*;
import dr.app.beauti.types.*;
import dr.app.beauti.util.XMLWriter;
import dr.evolution.datatype.DataType;
import dr.evolution.util.Taxa;
import dr.evolution.util.Units;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.tree.TreeModel;
import dr.evomodelxml.TreeWorkingPriorParsers;
import dr.evomodelxml.branchratemodel.*;
import dr.evomodelxml.coalescent.*;
import dr.evomodelxml.speciation.SpeciationLikelihoodParser;
import dr.evomodelxml.speciation.SpeciesTreeModelParser;
import dr.evomodelxml.speciation.YuleModelParser;
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
import java.util.EnumSet;
import java.util.List;

/**
 * @author Andrew Rambaut
 * @author Guy Baele
 * @version $Id$
 */
public class MarginalLikelihoodEstimationGenerator extends BaseComponentGenerator {

    public static final boolean DEBUG = false;

    private BeautiOptions beautiOptions = null;

    MarginalLikelihoodEstimationGenerator(final BeautiOptions options) {
        super(options);
        this.beautiOptions = options;
    }

    @Override
    public void checkOptions() throws GeneratorException {

        MarginalLikelihoodEstimationOptions mleOptions = (MarginalLikelihoodEstimationOptions)options.getComponentOptions(MarginalLikelihoodEstimationOptions.class);

        if (DEBUG) {
            System.out.println("mleOptions.performMLE: " + mleOptions.performMLE);
            System.out.println("mleOptions.performMLEGSS: " + mleOptions.performMLEGSS);
        }

        //++++++++++++++++ Improper priors ++++++++++++++++++
        if (mleOptions.performMLE) {
            for (Parameter param : options.selectParameters()) {
                if (param.isPriorImproper() || (param.priorType == PriorType.ONE_OVER_X_PRIOR && !param.getBaseName().contains("popSize"))) {
                    throw new GeneratorException("Parameter \"" + param.getName() + "\":" +
                            "\nhas an improper prior and will not sample correctly when estimating " +
                            "the marginal likelihood. " +
                            "\nPlease check the Prior panel.", BeautiFrame.PRIORS);
                }
            }
        }

        //++++++++++++++++ Coalescent Events available for GSS ++++++++++++++++++
        if (mleOptions.performMLEGSS) {
            EnumSet<TreePriorType> allowedTypes = EnumSet.of(
                    TreePriorType.CONSTANT, TreePriorType.EXPONENTIAL, TreePriorType.LOGISTIC, TreePriorType.EXPANSION, TreePriorType.SKYGRID, TreePriorType.GMRF_SKYRIDE, TreePriorType.YULE
            );
            EnumSet<TreePriorType> allowedMCMTypes = EnumSet.of(TreePriorType.CONSTANT, TreePriorType.EXPONENTIAL, TreePriorType.LOGISTIC, TreePriorType.EXPANSION);
            for (PartitionTreeModel model : options.getPartitionTreeModels()) {
                PartitionTreePrior prior = model.getPartitionTreePrior();
                if (!allowedTypes.contains(prior.getNodeHeightPrior())) {
                    throw new GeneratorException("Generalized stepping stone sampling can only be performed\n" +
                            "on standard parameteric coalescent tree priors and the Skyride and Skygrid models. " +
                            "\nPlease check the Trees panel.", BeautiFrame.TREES);
                }
                if (mleOptions.choiceTreeWorkingPrior.equals("Matching coalescent model") && !allowedMCMTypes.contains(prior.getNodeHeightPrior())) {
                    throw new GeneratorException("A Matching Coalescent Model cannot be constructed for\n" +
                            "the Skyride and Skygrid models. Please check the Marginal Likelihood\n" +
                            "Estimation settings via the MCMC panel.");
                }
            }

            // Shouldn't get here as the MLE switch in the MCMC tab already checks.
            for (AbstractPartitionData partition : options.getDataPartitions()) {
                if (partition.getDataType().getType() != DataType.NUCLEOTIDES) {
                    throw new GeneratorException(
                            "Generalized stepping-stone sampling is not currently\n" +
                                    "compatible with substitution models other than those\n" +
                                    "for nucleotide data. \n\n" +
                                    BeautiFrame.MCMC);
                }
            }

        }
    }

    public boolean usesInsertionPoint(final InsertionPoint point) {
        MarginalLikelihoodEstimationOptions component = (MarginalLikelihoodEstimationOptions) options.getComponentOptions(MarginalLikelihoodEstimationOptions.class);

        if (!component.performMLE && !component.performMLEGSS) {
            return false;
        }

        switch (point) {
            case AFTER_MCMC:
                return true;
            case IN_FILE_LOG_PARAMETERS:
                return options.logCoalescentEventsStatistic;
        }
        return false;
    }

    protected void generate(final InsertionPoint point, final Object item, final String prefix, final XMLWriter writer) {
        MarginalLikelihoodEstimationOptions component = (MarginalLikelihoodEstimationOptions) options.getComponentOptions(MarginalLikelihoodEstimationOptions.class);

        /*System.err.println("generate component: " + component);
        System.err.println("options.pathSteps: " + component.pathSteps);
        System.err.println("options.mleChainLength: " + component.mleChainLength);
        System.err.println("options.mleLogEvery: " + component.mleLogEvery);*/

        switch (point) {
            case AFTER_MCMC:
                writeMLE(writer, component);
                break;
            case IN_FILE_LOG_PARAMETERS:
                if (options.logCoalescentEventsStatistic) {
                    writeCoalescentEventsStatistic(writer);
                }
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
            attributes.add(new Attribute.Default<String>("resultsFileName", options.mleResultFileName));
            writer.writeOpenTag(PathSamplingAnalysis.PATH_SAMPLING_ANALYSIS, attributes);
            writer.writeTag("likelihoodColumn", new Attribute.Default<String>("name", "pathLikelihood.delta"), true);
            writer.writeTag("thetaColumn", new Attribute.Default<String>("name", "pathLikelihood.theta"), true);
            writer.writeCloseTag(PathSamplingAnalysis.PATH_SAMPLING_ANALYSIS);

            writer.writeComment("Stepping-stone sampling estimator from collected samples");
            attributes = new ArrayList<Attribute>();
            attributes.add(new Attribute.Default<String>("fileName", options.mleFileName));
            attributes.add(new Attribute.Default<String>("resultsFileName", options.mleResultFileName));
            writer.writeOpenTag(SteppingStoneSamplingAnalysis.STEPPING_STONE_SAMPLING_ANALYSIS, attributes);
            writer.writeTag("likelihoodColumn", new Attribute.Default<String>("name", "pathLikelihood.delta"), true);
            writer.writeTag("thetaColumn", new Attribute.Default<String>("name", "pathLikelihood.theta"), true);
            writer.writeCloseTag(SteppingStoneSamplingAnalysis.STEPPING_STONE_SAMPLING_ANALYSIS);

        } else if (options.performMLEGSS) {

            // TODO: does this need a prefix? I.e., will there ever be more than one of these?
            String modelPrefix = "";

            //First define necessary components for the tree working prior
            if (options.choiceTreeWorkingPrior.equals("Product of exponential distributions")) {
                //more general product of exponentials needs to be constructed

                if (DEBUG) {
                    System.err.println("productOfExponentials selected: " + options.choiceTreeWorkingPrior);
                }

                List<Attribute> attributes = new ArrayList<Attribute>();
                attributes.add(new Attribute.Default<String>(XMLParser.ID, "exponentials"));
                attributes.add(new Attribute.Default<String>("fileName", beautiOptions.logFileName));
                attributes.add(new Attribute.Default<String>("burnin", "" + (int)(beautiOptions.chainLength*0.10)));
                attributes.add(new Attribute.Default<String>("parameterColumn", "coalescentEventsStatistic"));
                attributes.add(new Attribute.Default<String>("dimension", "" + (beautiOptions.taxonList.getTaxonCount()-1)));

                writer.writeOpenTag(TreeWorkingPriorParsers.PRODUCT_OF_EXPONENTIALS_POSTERIOR_MEANS_LOESS, attributes);
                writer.writeIDref(TreeModel.TREE_MODEL, TreeModel.TREE_MODEL);
                writer.writeCloseTag(TreeWorkingPriorParsers.PRODUCT_OF_EXPONENTIALS_POSTERIOR_MEANS_LOESS);

            } else if (options.choiceTreeWorkingPrior.equals("Matching coalescent model")) {
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

                if (DEBUG) {
                    System.err.println("nodeHeightPrior: " + nodeHeightPrior);
                }

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

            } else {
                //matching speciation model has to be constructed
                //getting the speciation model
                if (DEBUG) {
                    System.err.println("matching speciation model selected: " + options.choiceTreeWorkingPrior);
                    System.err.println(beautiOptions.getPartitionTreePriors().get(0).getNodeHeightPrior());
                }
                TreePriorType nodeHeightPrior = beautiOptions.getPartitionTreePriors().get(0).getNodeHeightPrior();

                switch (nodeHeightPrior) {
                    case YULE:

                        writer.writeComment("A working prior for the Yule pure birth model.");
                        writer.writeOpenTag(
                                YuleModelParser.YULE_MODEL,
                                new Attribute[]{
                                        new Attribute.Default<String>(XMLParser.ID, modelPrefix + "yuleReference"),
                                        new Attribute.Default<String>("units", Units.Utils.getDefaultUnitName(beautiOptions.units))
                                }
                        );

                        writer.writeOpenTag(YuleModelParser.BIRTH_RATE);
                        writeParameter("yuleReference.birthRate", "yule.birthRate", beautiOptions.logFileName, (int) (options.mleChainLength * 0.10), writer);
                        writer.writeCloseTag(YuleModelParser.BIRTH_RATE);
                        writer.writeCloseTag(YuleModelParser.YULE_MODEL);

                        writer.writeComment("A working prior for the speciation process.");
                        writer.writeOpenTag(
                                SpeciationLikelihoodParser.SPECIATION_LIKELIHOOD,
                                new Attribute[]{
                                        new Attribute.Default<String>(XMLParser.ID, modelPrefix + "speciationReference")
                                }
                        );
                        writer.writeOpenTag(SpeciationLikelihoodParser.MODEL);
                        writer.writeIDref(YuleModelParser.YULE_MODEL, beautiOptions.getPartitionTreePriors().get(0).getPrefix() + "yuleReference");
                        writer.writeCloseTag(SpeciationLikelihoodParser.MODEL);
                        writer.writeOpenTag(SpeciesTreeModelParser.SPECIES_TREE);
                        writer.writeIDref(TreeModel.TREE_MODEL, modelPrefix + TreeModel.TREE_MODEL);
                        writer.writeCloseTag(SpeciesTreeModelParser.SPECIES_TREE);
                        writer.writeCloseTag(SpeciationLikelihoodParser.SPECIATION_LIKELIHOOD);

                        break;

                    default:

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

                            case JC:

                                if (codonPartitionCount > 1) {
                                    //write working priors for relative rates
                                    writeRelativeRates(writer, model, codonPartitionCount);
                                }
                                break;

                            case HKY:
                                if (codonPartitionCount > 1 && model.isUnlinkedSubstitutionModel()) {
                                    for (int i = 1; i <= codonPartitionCount; i++) {
                                        writer.writeOpenTag(WorkingPriorParsers.LOG_TRANSFORMED_NORMAL_REFERENCE_PRIOR,
                                                new Attribute[]{
                                                        new Attribute.Default<String>("fileName", beautiOptions.logFileName),
                                                        new Attribute.Default<String>("parameterColumn", model.getPrefix(i) + "kappa"),
                                                        new Attribute.Default<String>("burnin", "" + (int) (beautiOptions.chainLength * 0.10))
                                                });
                                        writer.writeIDref(ParameterParser.PARAMETER, model.getPrefix(i) + "kappa");
                                        writer.writeCloseTag(WorkingPriorParsers.LOG_TRANSFORMED_NORMAL_REFERENCE_PRIOR);
                                    }
                                } else {
                                    writer.writeOpenTag(WorkingPriorParsers.LOG_TRANSFORMED_NORMAL_REFERENCE_PRIOR,
                                            new Attribute[]{
                                                    new Attribute.Default<String>("fileName", beautiOptions.logFileName),
                                                    new Attribute.Default<String>("parameterColumn", model.getPrefix() + "kappa"),
                                                    new Attribute.Default<String>("burnin", "" + (int) (beautiOptions.chainLength * 0.10))
                                            });
                                    writer.writeIDref(ParameterParser.PARAMETER, model.getPrefix() + "kappa");
                                    writer.writeCloseTag(WorkingPriorParsers.LOG_TRANSFORMED_NORMAL_REFERENCE_PRIOR);
                                }
                                if (codonPartitionCount > 1) {
                                    //write working priors for relative rates
                                    writeRelativeRates(writer, model, codonPartitionCount);
                                }
                                break;

                            case TN93:
                                if (codonPartitionCount > 1 && model.isUnlinkedSubstitutionModel()) {
                                    for (int i = 1; i <= codonPartitionCount; i++) {
                                        writer.writeOpenTag(WorkingPriorParsers.LOG_TRANSFORMED_NORMAL_REFERENCE_PRIOR,
                                                new Attribute[]{
                                                        new Attribute.Default<String>("fileName", beautiOptions.logFileName),
                                                        new Attribute.Default<String>("parameterColumn", model.getPrefix(i) + "kappa1"),
                                                        new Attribute.Default<String>("burnin", "" + (int) (beautiOptions.chainLength * 0.10))
                                                });
                                        writer.writeIDref(ParameterParser.PARAMETER, model.getPrefix(i) + "kappa1");
                                        writer.writeCloseTag(WorkingPriorParsers.LOG_TRANSFORMED_NORMAL_REFERENCE_PRIOR);
                                        writer.writeOpenTag(WorkingPriorParsers.LOG_TRANSFORMED_NORMAL_REFERENCE_PRIOR,
                                                new Attribute[]{
                                                        new Attribute.Default<String>("fileName", beautiOptions.logFileName),
                                                        new Attribute.Default<String>("parameterColumn", model.getPrefix(i) + "kappa2"),
                                                        new Attribute.Default<String>("burnin", "" + (int) (beautiOptions.chainLength * 0.10))
                                                });
                                        writer.writeIDref(ParameterParser.PARAMETER, model.getPrefix(i) + "kappa2");
                                        writer.writeCloseTag(WorkingPriorParsers.LOG_TRANSFORMED_NORMAL_REFERENCE_PRIOR);
                                    }
                                } else {
                                    writer.writeOpenTag(WorkingPriorParsers.LOG_TRANSFORMED_NORMAL_REFERENCE_PRIOR,
                                            new Attribute[]{
                                                    new Attribute.Default<String>("fileName", beautiOptions.logFileName),
                                                    new Attribute.Default<String>("parameterColumn", model.getPrefix() + "kappa1"),
                                                    new Attribute.Default<String>("burnin", "" + (int) (beautiOptions.chainLength * 0.10))
                                            });
                                    writer.writeIDref(ParameterParser.PARAMETER, model.getPrefix() + "kappa1");
                                    writer.writeCloseTag(WorkingPriorParsers.LOG_TRANSFORMED_NORMAL_REFERENCE_PRIOR);
                                    writer.writeOpenTag(WorkingPriorParsers.LOG_TRANSFORMED_NORMAL_REFERENCE_PRIOR,
                                            new Attribute[]{
                                                    new Attribute.Default<String>("fileName", beautiOptions.logFileName),
                                                    new Attribute.Default<String>("parameterColumn", model.getPrefix() + "kappa2"),
                                                    new Attribute.Default<String>("burnin", "" + (int) (beautiOptions.chainLength * 0.10))
                                            });
                                    writer.writeIDref(ParameterParser.PARAMETER, model.getPrefix() + "kappa2");
                                    writer.writeCloseTag(WorkingPriorParsers.LOG_TRANSFORMED_NORMAL_REFERENCE_PRIOR);
                                }
                                if (codonPartitionCount > 1) {
                                    //write working priors for relative rates
                                    writeRelativeRates(writer, model, codonPartitionCount);
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
                                                            new Attribute.Default<String>("burnin", "" + (int) (beautiOptions.chainLength * 0.10))
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
                                                        new Attribute.Default<String>("burnin", "" + (int) (beautiOptions.chainLength * 0.10))
                                                });
                                        writer.writeIDref(ParameterParser.PARAMETER, model.getPrefix() + rateName);
                                        writer.writeCloseTag(WorkingPriorParsers.LOG_TRANSFORMED_NORMAL_REFERENCE_PRIOR);
                                    }
                                }
                                if (codonPartitionCount > 1) {
                                    //write working priors for relative rates
                                    writeRelativeRates(writer, model, codonPartitionCount);
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
                                                    new Attribute.Default<String>("burnin", "" + (int) (beautiOptions.chainLength * 0.10))
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
                                                new Attribute.Default<String>("burnin", "" + (int) (beautiOptions.chainLength * 0.10))
                                        });
                                writer.writeIDref(ParameterParser.PARAMETER, model.getPrefix() + "frequencies");
                                writer.writeCloseTag(WorkingPriorParsers.LOGIT_TRANSFORMED_NORMAL_REFERENCE_PRIOR);
                            }
                        }

                        break;//NUCLEOTIDES

                    case DataType.AMINO_ACIDS:

                    case DataType.TWO_STATES:

                    case DataType.COVARION:

                    case DataType.GENERAL:

                    case DataType.CONTINUOUS:

                    case DataType.MICRO_SAT:

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
                                            new Attribute.Default<String>("burnin", "" + (int) (beautiOptions.chainLength * 0.10))
                                    });
                            writer.writeIDref(ParameterParser.PARAMETER, model.getPrefix(i) + "alpha");
                            writer.writeCloseTag(WorkingPriorParsers.LOG_TRANSFORMED_NORMAL_REFERENCE_PRIOR);
                        }
                    } else {
                        writer.writeOpenTag(WorkingPriorParsers.LOG_TRANSFORMED_NORMAL_REFERENCE_PRIOR,
                                new Attribute[]{
                                        new Attribute.Default<String>("fileName", beautiOptions.logFileName),
                                        new Attribute.Default<String>("parameterColumn", model.getPrefix() + "alpha"),
                                        new Attribute.Default<String>("burnin", "" + (int) (beautiOptions.chainLength * 0.10))
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
                                            new Attribute.Default<String>("burnin", "" + (int) (beautiOptions.chainLength * 0.10))
                                    });
                            writer.writeIDref(ParameterParser.PARAMETER, model.getPrefix(i) + "pInv");
                            writer.writeCloseTag(WorkingPriorParsers.LOGIT_TRANSFORMED_NORMAL_REFERENCE_PRIOR);
                        }
                    } else {
                        writer.writeOpenTag(WorkingPriorParsers.LOGIT_TRANSFORMED_NORMAL_REFERENCE_PRIOR,
                                new Attribute[]{
                                        new Attribute.Default<String>("fileName", beautiOptions.logFileName),
                                        new Attribute.Default<String>("parameterColumn", model.getPrefix() + "pInv"),
                                        new Attribute.Default<String>("burnin", "" + (int) (beautiOptions.chainLength * 0.10))
                                });
                        writer.writeIDref(ParameterParser.PARAMETER, model.getPrefix() + "pInv");
                        writer.writeCloseTag(WorkingPriorParsers.LOGIT_TRANSFORMED_NORMAL_REFERENCE_PRIOR);
                    }
                }
            }

            //Continue with providing working priors for the clock model(s)
            for (PartitionClockModel model : beautiOptions.getPartitionClockModels()) {
                switch (model.getClockType()) {
                    case STRICT_CLOCK:
                        writer.writeOpenTag(WorkingPriorParsers.LOG_TRANSFORMED_NORMAL_REFERENCE_PRIOR,
                                new Attribute[]{
                                        new Attribute.Default<String>("fileName", beautiOptions.logFileName),
                                        new Attribute.Default<String>("parameterColumn", model.getPrefix() + "clock.rate"),
                                        new Attribute.Default<String>("burnin", "" + (int) (beautiOptions.chainLength * 0.10))
                                });
                        writer.writeIDref(ParameterParser.PARAMETER, model.getPrefix() + "clock.rate");
                        writer.writeCloseTag(WorkingPriorParsers.LOG_TRANSFORMED_NORMAL_REFERENCE_PRIOR);
                        writer.writeIDref(StrictClockBranchRatesParser.STRICT_CLOCK_BRANCH_RATES, model.getPrefix() + BranchRateModel.BRANCH_RATES);
                        break;

                    case UNCORRELATED:

                        if (model.isContinuousQuantile()) {
                            writer.writeIDref(ContinuousBranchRatesParser.CONTINUOUS_BRANCH_RATES, model.getPrefix() + BranchRateModel.BRANCH_RATES);
                        } else {
                            writer.writeIDref(DiscretizedBranchRatesParser.DISCRETIZED_BRANCH_RATES, model.getPrefix() + BranchRateModel.BRANCH_RATES);
                        }

                        switch (model.getClockDistributionType()) {
                            case GAMMA:
                                writer.writeOpenTag(WorkingPriorParsers.LOG_TRANSFORMED_NORMAL_REFERENCE_PRIOR,
                                        new Attribute[]{
                                                new Attribute.Default<String>("fileName", beautiOptions.logFileName),
                                                new Attribute.Default<String>("parameterColumn", model.getPrefix() + ClockType.UCGD_MEAN),
                                                new Attribute.Default<String>("burnin", "" + (int) (beautiOptions.chainLength * 0.10))
                                        });
                                writer.writeIDref(ParameterParser.PARAMETER, model.getPrefix() + ClockType.UCGD_MEAN);
                                writer.writeCloseTag(WorkingPriorParsers.LOG_TRANSFORMED_NORMAL_REFERENCE_PRIOR);
                                writer.writeOpenTag(WorkingPriorParsers.LOG_TRANSFORMED_NORMAL_REFERENCE_PRIOR,
                                        new Attribute[]{
                                                new Attribute.Default<String>("fileName", beautiOptions.logFileName),
                                                new Attribute.Default<String>("parameterColumn", model.getPrefix() + ClockType.UCGD_SHAPE),
                                                new Attribute.Default<String>("burnin", "" + (int) (beautiOptions.chainLength * 0.10))
                                        });
                                writer.writeIDref(ParameterParser.PARAMETER, model.getPrefix() + ClockType.UCGD_SHAPE);
                                writer.writeCloseTag(WorkingPriorParsers.LOG_TRANSFORMED_NORMAL_REFERENCE_PRIOR);
                                break;

                            case LOGNORMAL:
                                if (!model.getClockRateParameter().isInRealSpace()) {
                                    writer.writeOpenTag(WorkingPriorParsers.LOG_TRANSFORMED_NORMAL_REFERENCE_PRIOR,
                                            new Attribute[]{
                                                    new Attribute.Default<String>("fileName", beautiOptions.logFileName),
                                                    new Attribute.Default<String>("parameterColumn", model.getPrefix() + ClockType.UCLD_MEAN),
                                                    new Attribute.Default<String>("burnin", "" + (int) (beautiOptions.chainLength * 0.10))
                                            });
                                    writer.writeIDref(ParameterParser.PARAMETER, model.getPrefix() + ClockType.UCLD_MEAN);
                                    writer.writeCloseTag(WorkingPriorParsers.LOG_TRANSFORMED_NORMAL_REFERENCE_PRIOR);
                                    writer.writeOpenTag(WorkingPriorParsers.LOG_TRANSFORMED_NORMAL_REFERENCE_PRIOR,
                                            new Attribute[]{
                                                    new Attribute.Default<String>("fileName", beautiOptions.logFileName),
                                                    new Attribute.Default<String>("parameterColumn", model.getPrefix() + ClockType.UCLD_STDEV),
                                                    new Attribute.Default<String>("burnin", "" + (int) (beautiOptions.chainLength * 0.10))
                                            });
                                    writer.writeIDref(ParameterParser.PARAMETER, model.getPrefix() + ClockType.UCLD_STDEV);
                                    writer.writeCloseTag(WorkingPriorParsers.LOG_TRANSFORMED_NORMAL_REFERENCE_PRIOR);
                                } else {
                                    writer.writeOpenTag(WorkingPriorParsers.NORMAL_REFERENCE_PRIOR,
                                            new Attribute[]{
                                                    new Attribute.Default<String>("fileName", beautiOptions.logFileName),
                                                    new Attribute.Default<String>("parameterColumn", model.getPrefix() + ClockType.UCLD_MEAN),
                                                    new Attribute.Default<String>("burnin", "" + (int) (beautiOptions.chainLength * 0.10))
                                            });
                                    writer.writeIDref(ParameterParser.PARAMETER, model.getPrefix() + ClockType.UCLD_MEAN);
                                    writer.writeCloseTag(WorkingPriorParsers.NORMAL_REFERENCE_PRIOR);
                                    writer.writeOpenTag(WorkingPriorParsers.LOG_TRANSFORMED_NORMAL_REFERENCE_PRIOR,
                                            new Attribute[]{
                                                    new Attribute.Default<String>("fileName", beautiOptions.logFileName),
                                                    new Attribute.Default<String>("parameterColumn", model.getPrefix() + ClockType.UCLD_STDEV),
                                                    new Attribute.Default<String>("burnin", "" + (int) (beautiOptions.chainLength * 0.10))
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
                                                new Attribute.Default<String>("burnin", "" + (int) (beautiOptions.chainLength * 0.10))
                                        });
                                writer.writeIDref(ParameterParser.PARAMETER, model.getPrefix() + ClockType.UCED_MEAN);
                                writer.writeCloseTag(WorkingPriorParsers.LOG_TRANSFORMED_NORMAL_REFERENCE_PRIOR);
                                break;
                            case MODEL_AVERAGING:
                                throw new RuntimeException("Marginal likelihood estimation cannot be performed on a clock model that performs model averaging.");
                        }
                        break;

                    case FIXED_LOCAL_CLOCK:
                        writer.writeOpenTag(WorkingPriorParsers.LOG_TRANSFORMED_NORMAL_REFERENCE_PRIOR,
                                new Attribute[]{
                                        new Attribute.Default<String>("fileName", beautiOptions.logFileName),
                                        new Attribute.Default<String>("parameterColumn", model.getPrefix() + "clock.rate"),
                                        new Attribute.Default<String>("burnin", "" + (int) (beautiOptions.chainLength * 0.10))
                                });
                        writer.writeIDref(ParameterParser.PARAMETER, model.getPrefix() + "clock.rate");
                        writer.writeCloseTag(WorkingPriorParsers.LOG_TRANSFORMED_NORMAL_REFERENCE_PRIOR);

                        for (Taxa taxonSet : beautiOptions.taxonSets) {
                            if (beautiOptions.taxonSetsMono.get(taxonSet)) {
                                String parameterName = taxonSet.getId() + ".rate";

                                writer.writeOpenTag(WorkingPriorParsers.LOG_TRANSFORMED_NORMAL_REFERENCE_PRIOR,
                                        new Attribute[]{
                                                new Attribute.Default<String>("fileName", beautiOptions.logFileName),
                                                new Attribute.Default<String>("parameterColumn", model.getPrefix() + parameterName),
                                                new Attribute.Default<String>("burnin", "" + (int) (beautiOptions.chainLength * 0.10))
                                        });
                                writer.writeIDref(ParameterParser.PARAMETER, model.getPrefix() + parameterName);
                                writer.writeCloseTag(WorkingPriorParsers.LOG_TRANSFORMED_NORMAL_REFERENCE_PRIOR);
                            }
                        }

                        writer.writeIDref(LocalClockModelParser.LOCAL_CLOCK_MODEL, model.getPrefix() + BranchRateModel.BRANCH_RATES);
                        break;

                    case RANDOM_LOCAL_CLOCK:

                        //TODO



                        writer.writeIDref(RandomLocalClockModelParser.LOCAL_BRANCH_RATES, model.getPrefix() + BranchRateModel.BRANCH_RATES);
                        break;

                    default:
                        throw new IllegalArgumentException("Unknown clock model");
                }
            }

            //Provide working priors for the coalescent model(s)
            for (PartitionTreePrior model : beautiOptions.getPartitionTreePriors()) {
                TreePriorType nodeHeightPrior = model.getNodeHeightPrior();
                TreePriorParameterizationType parameterization = model.getParameterization();

                if (DEBUG) {
                    System.err.println("nodeHeightPrior: " + nodeHeightPrior);
                }

                switch (nodeHeightPrior) {
                    case CONSTANT:
                        writer.writeOpenTag(WorkingPriorParsers.LOG_TRANSFORMED_NORMAL_REFERENCE_PRIOR,
                                new Attribute[]{
                                        new Attribute.Default<String>("fileName", beautiOptions.logFileName),
                                        new Attribute.Default<String>("parameterColumn", "constant.popSize"),
                                        new Attribute.Default<String>("burnin", "" + (int) (beautiOptions.chainLength * 0.10))
                                });
                        writer.writeIDref(ParameterParser.PARAMETER, "constant.popSize");
                        writer.writeCloseTag(WorkingPriorParsers.LOG_TRANSFORMED_NORMAL_REFERENCE_PRIOR);
                        break;

                    case EXPONENTIAL:
                        writer.writeOpenTag(WorkingPriorParsers.LOG_TRANSFORMED_NORMAL_REFERENCE_PRIOR,
                                new Attribute[]{
                                        new Attribute.Default<String>("fileName", beautiOptions.logFileName),
                                        new Attribute.Default<String>("parameterColumn", "exponential.popSize"),
                                        new Attribute.Default<String>("burnin", "" + (int) (beautiOptions.chainLength * 0.10))
                                });
                        writer.writeIDref(ParameterParser.PARAMETER, "exponential.popSize");
                        writer.writeCloseTag(WorkingPriorParsers.LOG_TRANSFORMED_NORMAL_REFERENCE_PRIOR);

                        if (parameterization == TreePriorParameterizationType.GROWTH_RATE) {
                            writer.writeOpenTag(WorkingPriorParsers.NORMAL_REFERENCE_PRIOR,
                                    new Attribute[]{
                                            new Attribute.Default<String>("fileName", beautiOptions.logFileName),
                                            new Attribute.Default<String>("parameterColumn", "exponential.growthRate"),
                                            new Attribute.Default<String>("burnin", "" + (int) (beautiOptions.chainLength * 0.10))
                                    });
                            writer.writeIDref(ParameterParser.PARAMETER, "exponential.growthRate");
                            writer.writeCloseTag(WorkingPriorParsers.NORMAL_REFERENCE_PRIOR);
                        } else {
                            writer.writeOpenTag(WorkingPriorParsers.LOG_TRANSFORMED_NORMAL_REFERENCE_PRIOR,
                                    new Attribute[]{
                                            new Attribute.Default<String>("fileName", beautiOptions.logFileName),
                                            new Attribute.Default<String>("parameterColumn", "exponential.doublingTime"),
                                            new Attribute.Default<String>("burnin", "" + (int) (beautiOptions.chainLength * 0.10))
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
                                        new Attribute.Default<String>("burnin", "" + (int) (beautiOptions.chainLength * 0.10))
                                });
                        writer.writeIDref(ParameterParser.PARAMETER, "logistic.popSize");
                        writer.writeCloseTag(WorkingPriorParsers.LOG_TRANSFORMED_NORMAL_REFERENCE_PRIOR);

                        if (parameterization == TreePriorParameterizationType.GROWTH_RATE) {
                            writer.writeOpenTag(WorkingPriorParsers.NORMAL_REFERENCE_PRIOR,
                                    new Attribute[]{
                                            new Attribute.Default<String>("fileName", beautiOptions.logFileName),
                                            new Attribute.Default<String>("parameterColumn", "logistic.growthRate"),
                                            new Attribute.Default<String>("burnin", "" + (int) (beautiOptions.chainLength * 0.10))
                                    });
                            writer.writeIDref(ParameterParser.PARAMETER, "logistic.growthRate");
                            writer.writeCloseTag(WorkingPriorParsers.NORMAL_REFERENCE_PRIOR);
                        } else {
                            writer.writeOpenTag(WorkingPriorParsers.LOG_TRANSFORMED_NORMAL_REFERENCE_PRIOR,
                                    new Attribute[]{
                                            new Attribute.Default<String>("fileName", beautiOptions.logFileName),
                                            new Attribute.Default<String>("parameterColumn", "logistic.doublingTime"),
                                            new Attribute.Default<String>("burnin", "" + (int) (beautiOptions.chainLength * 0.10))
                                    });
                            writer.writeIDref(ParameterParser.PARAMETER, "logistic.doublingTime");
                            writer.writeCloseTag(WorkingPriorParsers.LOG_TRANSFORMED_NORMAL_REFERENCE_PRIOR);
                        }

                        writer.writeOpenTag(WorkingPriorParsers.LOG_TRANSFORMED_NORMAL_REFERENCE_PRIOR,
                                new Attribute[]{
                                        new Attribute.Default<String>("fileName", beautiOptions.logFileName),
                                        new Attribute.Default<String>("parameterColumn", "logistic.t50"),
                                        new Attribute.Default<String>("burnin", "" + (int) (beautiOptions.chainLength * 0.10))
                                });
                        writer.writeIDref(ParameterParser.PARAMETER, "logistic.t50");
                        writer.writeCloseTag(WorkingPriorParsers.LOG_TRANSFORMED_NORMAL_REFERENCE_PRIOR);

                        break;

                    case EXPANSION:
                        writer.writeOpenTag(WorkingPriorParsers.LOG_TRANSFORMED_NORMAL_REFERENCE_PRIOR,
                                new Attribute[]{
                                        new Attribute.Default<String>("fileName", beautiOptions.logFileName),
                                        new Attribute.Default<String>("parameterColumn", "expansion.popSize"),
                                        new Attribute.Default<String>("burnin", "" + (int) (beautiOptions.chainLength * 0.10))
                                });
                        writer.writeIDref(ParameterParser.PARAMETER, "expansion.popSize");
                        writer.writeCloseTag(WorkingPriorParsers.LOG_TRANSFORMED_NORMAL_REFERENCE_PRIOR);

                        if (parameterization == TreePriorParameterizationType.GROWTH_RATE) {
                            writer.writeOpenTag(WorkingPriorParsers.NORMAL_REFERENCE_PRIOR,
                                    new Attribute[]{
                                            new Attribute.Default<String>("fileName", beautiOptions.logFileName),
                                            new Attribute.Default<String>("parameterColumn", "expansion.growthRate"),
                                            new Attribute.Default<String>("burnin", "" + (int) (beautiOptions.chainLength * 0.10))
                                    });
                            writer.writeIDref(ParameterParser.PARAMETER, "expansion.growthRate");
                            writer.writeCloseTag(WorkingPriorParsers.NORMAL_REFERENCE_PRIOR);
                        } else {
                            writer.writeOpenTag(WorkingPriorParsers.LOG_TRANSFORMED_NORMAL_REFERENCE_PRIOR,
                                    new Attribute[]{
                                            new Attribute.Default<String>("fileName", beautiOptions.logFileName),
                                            new Attribute.Default<String>("parameterColumn", "expansion.doublingTime"),
                                            new Attribute.Default<String>("burnin", "" + (int) (beautiOptions.chainLength * 0.10))
                                    });
                            writer.writeIDref(ParameterParser.PARAMETER, "expansion.doublingTime");
                            writer.writeCloseTag(WorkingPriorParsers.LOG_TRANSFORMED_NORMAL_REFERENCE_PRIOR);
                        }

                        writer.writeOpenTag(WorkingPriorParsers.LOGIT_TRANSFORMED_NORMAL_REFERENCE_PRIOR,
                                new Attribute[]{
                                        new Attribute.Default<String>("fileName", beautiOptions.logFileName),
                                        new Attribute.Default<String>("parameterColumn", "expansion.ancestralProportion"),
                                        new Attribute.Default<String>("burnin", "" + (int) (beautiOptions.chainLength * 0.10))
                                });
                        writer.writeIDref(ParameterParser.PARAMETER, "expansion.ancestralProportion");
                        writer.writeCloseTag(WorkingPriorParsers.LOGIT_TRANSFORMED_NORMAL_REFERENCE_PRIOR);

                        break;

                    case GMRF_SKYRIDE:
                        writer.writeOpenTag(WorkingPriorParsers.NORMAL_REFERENCE_PRIOR,
                                new Attribute[]{
                                        new Attribute.Default<String>("fileName", beautiOptions.logFileName),
                                        new Attribute.Default<String>("parameterColumn", "skyride.logPopSize"),
                                        new Attribute.Default<Integer>("dimension", beautiOptions.taxonList.getTaxonCount() - 1),
                                        new Attribute.Default<String>("burnin", "" + (int) (beautiOptions.chainLength * 0.10))
                                });
                        writer.writeIDref(ParameterParser.PARAMETER, "skyride.logPopSize");
                        writer.writeCloseTag(WorkingPriorParsers.NORMAL_REFERENCE_PRIOR);

                        writer.writeOpenTag(WorkingPriorParsers.LOG_TRANSFORMED_NORMAL_REFERENCE_PRIOR,
                                new Attribute[]{
                                        new Attribute.Default<String>("fileName", beautiOptions.logFileName),
                                        new Attribute.Default<String>("parameterColumn", "skyride.precision"),
                                        new Attribute.Default<String>("burnin", "" + (int) (beautiOptions.chainLength * 0.10))
                                });
                        writer.writeIDref(ParameterParser.PARAMETER, "skyride.precision");
                        writer.writeCloseTag(WorkingPriorParsers.LOG_TRANSFORMED_NORMAL_REFERENCE_PRIOR);

                        break;

                    case SKYGRID:
                        writer.writeOpenTag(WorkingPriorParsers.NORMAL_REFERENCE_PRIOR,
                                new Attribute[]{
                                        new Attribute.Default<String>("fileName", beautiOptions.logFileName),
                                        new Attribute.Default<String>("parameterColumn", "skygrid.logPopSize"),
                                        new Attribute.Default<Integer>("dimension", model.getSkyGridCount()),
                                        new Attribute.Default<String>("burnin", "" + (int) (beautiOptions.chainLength * 0.10))
                                });
                        writer.writeIDref(ParameterParser.PARAMETER, "skygrid.logPopSize");
                        writer.writeCloseTag(WorkingPriorParsers.NORMAL_REFERENCE_PRIOR);

                        writer.writeOpenTag(WorkingPriorParsers.LOG_TRANSFORMED_NORMAL_REFERENCE_PRIOR,
                                new Attribute[]{
                                        new Attribute.Default<String>("fileName", beautiOptions.logFileName),
                                        new Attribute.Default<String>("parameterColumn", "skygrid.precision"),
                                        new Attribute.Default<String>("burnin", "" + (int) (beautiOptions.chainLength * 0.10))
                                });
                        writer.writeIDref(ParameterParser.PARAMETER, "skygrid.precision");
                        writer.writeCloseTag(WorkingPriorParsers.LOG_TRANSFORMED_NORMAL_REFERENCE_PRIOR);

                        break;

                    case YULE:
                        writer.writeOpenTag(WorkingPriorParsers.LOG_TRANSFORMED_NORMAL_REFERENCE_PRIOR,
                                new Attribute[]{
                                        new Attribute.Default<String>("fileName", beautiOptions.logFileName),
                                        new Attribute.Default<String>("parameterColumn", "yule.birthRate"),
                                        new Attribute.Default<String>("burnin", "" + (int) (beautiOptions.chainLength * 0.10))
                                });
                        writer.writeIDref(ParameterParser.PARAMETER, "yule.birthRate");
                        writer.writeCloseTag(WorkingPriorParsers.LOG_TRANSFORMED_NORMAL_REFERENCE_PRIOR);

                        break;

                }
            }

            //TODO: take care of anything else I missed

            if (options.choiceTreeWorkingPrior.equals("Product of exponential distributions")) {
                writer.writeIDref("productOfExponentialsPosteriorMeansLoess", "exponentials");
            } else if (options.choiceTreeWorkingPrior.equals("Matching coalescent model")) {
                writer.writeIDref(CoalescentLikelihoodParser.COALESCENT_LIKELIHOOD, "coalescentReference");
            } else {
                writer.writeIDref(YuleModelParser.YULE_MODEL, "yuleReference");
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
            attributes.add(new Attribute.Default<String>("resultsFileName", options.mleResultFileName));
            writer.writeOpenTag(GeneralizedSteppingStoneSamplingAnalysis.GENERALIZED_STEPPING_STONE_SAMPLING_ANALYSIS, attributes);
            writer.writeTag("sourceColumn", new Attribute.Default<String>("name", "pathLikelihood.source"), true);
            writer.writeTag("destinationColumn", new Attribute.Default<String>("name", "pathLikelihood.destination"), true);
            writer.writeTag("thetaColumn", new Attribute.Default<String>("name", "pathLikelihood.theta"), true);
            writer.writeCloseTag(GeneralizedSteppingStoneSamplingAnalysis.GENERALIZED_STEPPING_STONE_SAMPLING_ANALYSIS);

        }

    }

    private void writeRelativeRates(XMLWriter writer, PartitionSubstitutionModel model, int codonPartitionCount) {
        for (int i = 1; i <= codonPartitionCount; i++) {
            writer.writeOpenTag(WorkingPriorParsers.LOG_TRANSFORMED_NORMAL_REFERENCE_PRIOR,
                    new Attribute[]{
                            new Attribute.Default<String>("fileName", beautiOptions.logFileName),
                            new Attribute.Default<String>("parameterColumn", model.getPrefix(i) + "mu"),
                            new Attribute.Default<String>("burnin", "" + (int)(beautiOptions.chainLength*0.10)),
                            new Attribute.Default<String>("upperLimit", "" + (double)(codonPartitionCount))
                    });
            writer.writeIDref(ParameterParser.PARAMETER, model.getPrefix(i) + "mu");
            writer.writeCloseTag(WorkingPriorParsers.LOG_TRANSFORMED_NORMAL_REFERENCE_PRIOR);
        }
    }

    private void writeCoalescentEventsStatistic(XMLWriter writer) {
        writer.writeOpenTag("coalescentEventsStatistic");
        // coalescentLikelihood
        for (PartitionTreeModel model : options.getPartitionTreeModels()) {
            PartitionTreePrior prior = model.getPartitionTreePrior();
            TreePriorGenerator.writePriorLikelihoodReferenceLog(prior, model, writer);
            writer.writeText("");
        }

            /*for (PartitionTreePrior prior : options.getPartitionTreePriors()) {
                if (prior.getNodeHeightPrior() == TreePriorType.EXTENDED_SKYLINE) {
                    writer.writeIDref(CoalescentLikelihoodParser.COALESCENT_LIKELIHOOD, prior.getPrefix() + COALESCENT); // only 1 coalescent
                } else if (prior.getNodeHeightPrior() == TreePriorType.SKYGRID) {
                    writer.writeIDref(GMRFSkyrideLikelihoodParser.SKYGRID_LIKELIHOOD, prior.getPrefix() + "skygrid");
                }
            }*/
        writer.writeCloseTag("coalescentEventsStatistic");
    }


    private void writeParameterIdref(XMLWriter writer, Parameter parameter) {
        if (parameter.isStatistic) {
            writer.writeIDref("statistic", parameter.getName());
        } else {
            writer.writeIDref(ParameterParser.PARAMETER, parameter.getName());
        }
    }

}
