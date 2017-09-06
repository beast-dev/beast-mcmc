/*
 * TreePriorGenerator.java
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
import dr.app.beauti.types.StartingTreeType;
import dr.app.beauti.types.TreePriorParameterizationType;
import dr.app.beauti.types.TreePriorType;
import dr.app.beauti.util.XMLWriter;
import dr.evolution.util.Taxa;
import dr.evolution.util.Units;
import dr.evomodel.tree.TreeModel;
import dr.evomodelxml.CSVExporterParser;
import dr.evomodelxml.coalescent.*;
import dr.evomodelxml.speciation.*;
import dr.evoxml.TaxaParser;
import dr.inference.distribution.ExponentialDistributionModel;
import dr.inference.distribution.ExponentialMarkovModel;
import dr.inference.model.ParameterParser;
import dr.inferencexml.distribution.DistributionModelParser;
import dr.inferencexml.distribution.ExponentialMarkovModelParser;
import dr.inferencexml.distribution.MixedDistributionLikelihoodParser;
import dr.inferencexml.model.SumStatisticParser;
import dr.util.Attribute;
import dr.xml.XMLParser;

/**
 * @author Alexei Drummond
 * @author Andrew Rambaut
 */
public class TreePriorGenerator extends Generator {

    public TreePriorGenerator(BeautiOptions options, ComponentFactory[] components) {
        super(options, components);
    }

//    void writeTreePrior(PartitionTreePrior prior, PartitionTreeModel model, XMLWriter writer) {    // for species, partitionName.treeModel
//        setModelPrefix(prior.getPrefix()); // only has prefix, if (options.getPartitionTreePriors().size() > 1)
//
//        writePriorLikelihood(prior, model, writer);
//    }

    /**
     * Write a tree prior (coalescent or speciational) model
     *
     * @param prior  the partition tree prior
     * @param writer the writer
     */
    void writeTreePriorModel(PartitionTreePrior prior, XMLWriter writer) {

        String prefix = prior.getPrefix();

        String initialPopSize = null;

        TreePriorType nodeHeightPrior = prior.getNodeHeightPrior();
        Units.Type units = options.units;
        TreePriorParameterizationType parameterization = prior.getParameterization();

        switch (nodeHeightPrior) {
            case CONSTANT:
                writer.writeComment("A prior assumption that the population size has remained constant",
                        "throughout the time spanned by the genealogy.");
                writer.writeOpenTag(
                        ConstantPopulationModelParser.CONSTANT_POPULATION_MODEL,
                        new Attribute[]{
                                new Attribute.Default<String>(XMLParser.ID, prefix + "constant"),
                                new Attribute.Default<String>("units", Units.Utils.getDefaultUnitName(options.units))
                        }
                );

                writer.writeOpenTag(ConstantPopulationModelParser.POPULATION_SIZE);
                writeParameter("constant.popSize", prior, writer);
                writer.writeCloseTag(ConstantPopulationModelParser.POPULATION_SIZE);
                writer.writeCloseTag(ConstantPopulationModelParser.CONSTANT_POPULATION_MODEL);

                break;

            case EXPONENTIAL:
                // generate an exponential prior tree

                writer.writeComment("A prior assumption that the population size has grown exponentially",
                        "throughout the time spanned by the genealogy.");
                writer.writeOpenTag(
                        ExponentialGrowthModelParser.EXPONENTIAL_GROWTH_MODEL,
                        new Attribute[]{
                                new Attribute.Default<String>(XMLParser.ID, prefix + "exponential"),
                                new Attribute.Default<String>("units", Units.Utils.getDefaultUnitName(options.units))
                        }
                );

                // write pop size socket
                writer.writeOpenTag(ExponentialGrowthModelParser.POPULATION_SIZE);
                writeParameter("exponential.popSize", prior, writer);
                writer.writeCloseTag(ExponentialGrowthModelParser.POPULATION_SIZE);

                if (parameterization == TreePriorParameterizationType.GROWTH_RATE) {
                    // write growth rate socket
                    writer.writeOpenTag(ExponentialGrowthModelParser.GROWTH_RATE);
                    writeParameter("exponential.growthRate", prior, writer);
                    writer.writeCloseTag(ExponentialGrowthModelParser.GROWTH_RATE);
                } else {
                    // write doubling time socket
                    writer.writeOpenTag(ExponentialGrowthModelParser.DOUBLING_TIME);
                    writeParameter("exponential.doublingTime", prior, writer);
                    writer.writeCloseTag(ExponentialGrowthModelParser.DOUBLING_TIME);
                }

                writer.writeCloseTag(ExponentialGrowthModelParser.EXPONENTIAL_GROWTH_MODEL);

                break;

            case LOGISTIC:
                // generate an exponential prior tree
                writer.writeComment("A prior assumption that the population size has grown logistically",
                        "throughout the time spanned by the genealogy.");
                writer.writeOpenTag(
                        LogisticGrowthModelParser.LOGISTIC_GROWTH_MODEL,
                        new Attribute[]{
                                new Attribute.Default<String>(XMLParser.ID, prefix + "logistic"),
                                new Attribute.Default<String>("units", Units.Utils.getDefaultUnitName(options.units))
                        }
                );

                // write pop size socket
                writer.writeOpenTag(LogisticGrowthModelParser.POPULATION_SIZE);
                writeParameter("logistic.popSize", prior, writer);
                writer.writeCloseTag(LogisticGrowthModelParser.POPULATION_SIZE);

                if (parameterization == TreePriorParameterizationType.GROWTH_RATE) {
                    // write growth rate socket
                    writer.writeOpenTag(LogisticGrowthModelParser.GROWTH_RATE);
                    writeParameter("logistic.growthRate", prior, writer);
                    writer.writeCloseTag(LogisticGrowthModelParser.GROWTH_RATE);
                } else {
                    // write doubling time socket
                    writer.writeOpenTag(LogisticGrowthModelParser.DOUBLING_TIME);
                    writeParameter("logistic.doublingTime", prior, writer);
                    writer.writeCloseTag(LogisticGrowthModelParser.DOUBLING_TIME);
                }

                // write logistic t50 socket
                writer.writeOpenTag(LogisticGrowthModelParser.TIME_50);

//	            if (options.clockModelOptions.getRateOptionClockModel() == FixRateType.FIX_MEAN
//	        			|| options.clockModelOptions.getRateOptionClockModel() == FixRateType.RELATIVE_TO) {
//            		writer.writeComment("No calibration");
//                writer.writeComment("logistic.t50 initial always has to < treeRootHeight initial");
//                dr.app.beauti.options.Parameter priorPara = prior.getParameter("logistic.t50");
//
//                double initRootHeight;
//                if (options.isShareSameTreePrior()) {
//                    initRootHeight = priorPara.initial;
//                    for (PartitionTreeModel tree : options.getPartitionTreeModels()) {
//                        double tmpRootHeight = tree.getParameter("treeModel.rootHeight").initial;
//                        if (initRootHeight > tmpRootHeight) { // take min
//                            initRootHeight = tmpRootHeight;
//                        }
//                    }
//                } else {
//                    initRootHeight = prior.getTreeModel().getParameter("treeModel.rootHeight").initial;
//                }
//                // logistic.t50 initial always has to < treeRootHeight initial
//                if (priorPara.initial >= initRootHeight) {
//                    priorPara.initial = initRootHeight / 2; // tree prior.initial has to < treeRootHeight.initial
//                }
//	            } else {
//	            	writer.writeComment("Has calibration");
//
//	            	throw new IllegalArgumentException("This function is not available in this release !");
//	            }

                writeParameter("logistic.t50", prior, writer);
                writer.writeCloseTag(LogisticGrowthModelParser.TIME_50);

                writer.writeCloseTag(LogisticGrowthModelParser.LOGISTIC_GROWTH_MODEL);

                initialPopSize = "logistic.popSize";

                break;

            case EXPANSION:
                // generate an exponential prior tree
                writer.writeComment("A prior assumption that the population size has grown exponentially",
                        "from some ancestral population size in the past.");
                writer.writeOpenTag(
                        ExpansionModelParser.EXPANSION_MODEL,
                        new Attribute[]{
                                new Attribute.Default<String>(XMLParser.ID, prefix + "expansion"),
                                new Attribute.Default<String>("units", Units.Utils.getDefaultUnitName(options.units))
                        }
                );

                // write pop size socket
                writeParameter(ExpansionModelParser.POPULATION_SIZE, "expansion.popSize", prior, writer);

                if (parameterization == TreePriorParameterizationType.GROWTH_RATE) {
                    // write growth rate socket
                    writeParameter(ExpansionModelParser.GROWTH_RATE, "expansion.growthRate", prior, writer);
                } else {
                    // write doubling time socket
                    writeParameter(ExpansionModelParser.DOUBLING_TIME, "expansion.doublingTime", prior, writer);
                }

                // write ancestral proportion socket
                writeParameter(ExpansionModelParser.ANCESTRAL_POPULATION_PROPORTION, "expansion.ancestralProportion", prior, writer);

                writer.writeCloseTag(ExpansionModelParser.EXPANSION_MODEL);

                initialPopSize = "expansion.popSize";

                break;

            case YULE:
            case YULE_CALIBRATION:
                if (nodeHeightPrior == TreePriorType.YULE_CALIBRATION) {
                    writer.writeComment("Calibrated Yule: Heled J, Drummond AJ (2011), Syst Biol, doi: " +
                            "10.1093/sysbio/syr087");
                } else {
                    writer.writeComment("A prior on the distribution node heights defined given",
                            "a Yule speciation process (a pure birth process).");
                }
                writer.writeOpenTag(
                        YuleModelParser.YULE_MODEL,
                        new Attribute[]{
                                new Attribute.Default<String>(XMLParser.ID, prefix + YuleModelParser.YULE),
                                new Attribute.Default<String>("units", Units.Utils.getDefaultUnitName(units))
                        }
                );

                writeParameter(YuleModelParser.BIRTH_RATE, "yule.birthRate", prior, writer);
                writer.writeCloseTag(YuleModelParser.YULE_MODEL);

                break;

            case BIRTH_DEATH:
            case BIRTH_DEATH_INCOMPLETE_SAMPLING:
                writer.writeComment("A prior on the distribution node heights defined given");
                writer.writeComment(nodeHeightPrior == TreePriorType.BIRTH_DEATH_INCOMPLETE_SAMPLING ?
                        BirthDeathModelParser.getCitationRHO() : BirthDeathModelParser.getCitation());
                writer.writeOpenTag(
                        BirthDeathModelParser.BIRTH_DEATH_MODEL,
                        new Attribute[]{
                                new Attribute.Default<String>(XMLParser.ID, prefix + BirthDeathModelParser.BIRTH_DEATH),
                                new Attribute.Default<String>("units", Units.Utils.getDefaultUnitName(units))
                        }
                );

                writeParameter(BirthDeathModelParser.BIRTHDIFF_RATE, BirthDeathModelParser.MEAN_GROWTH_RATE_PARAM_NAME, prior, writer);
                writeParameter(BirthDeathModelParser.RELATIVE_DEATH_RATE, BirthDeathModelParser.RELATIVE_DEATH_RATE_PARAM_NAME, prior, writer);

                if (nodeHeightPrior == TreePriorType.BIRTH_DEATH_INCOMPLETE_SAMPLING) {
                    writeParameter(BirthDeathModelParser.SAMPLE_PROB,
                            BirthDeathModelParser.BIRTH_DEATH + "." + BirthDeathModelParser.SAMPLE_PROB, prior, writer);
                }

                writer.writeCloseTag(BirthDeathModelParser.BIRTH_DEATH_MODEL);

                break;

            case BIRTH_DEATH_SERIAL_SAMPLING:
                writer.writeComment(BirthDeathSerialSamplingModelParser.getCitationPsiOrg());

                writer.writeOpenTag(
                        BirthDeathSerialSamplingModelParser.BIRTH_DEATH_SERIAL_MODEL,
                        new Attribute[]{
                                new Attribute.Default<String>(XMLParser.ID, prefix + BirthDeathSerialSamplingModelParser.BDSS),
                                new Attribute.Default<String>("units", Units.Utils.getDefaultUnitName(units)),
                                new Attribute.Default<Boolean>(BirthDeathSerialSamplingModelParser.HAS_FINAL_SAMPLE, false)
                        }
                );

                writeParameter(BirthDeathSerialSamplingModelParser.LAMBDA,
                        BirthDeathSerialSamplingModelParser.BDSS + "." + BirthDeathSerialSamplingModelParser.LAMBDA, prior, writer);
                writeParameter(BirthDeathSerialSamplingModelParser.RELATIVE_MU,
                        BirthDeathSerialSamplingModelParser.BDSS + "." + BirthDeathSerialSamplingModelParser.RELATIVE_MU, prior, writer);
//                writeParameter(BirthDeathSerialSamplingModelParser.SAMPLE_PROBABILITY,
//                        BirthDeathSerialSamplingModelParser.BDSS + "." + BirthDeathSerialSamplingModelParser.SAMPLE_PROBABILITY, prior, writer);
                writeParameter(BirthDeathSerialSamplingModelParser.PSI,
                        BirthDeathSerialSamplingModelParser.BDSS + "." + BirthDeathSerialSamplingModelParser.PSI, prior, writer);
                writeParameter(BirthDeathSerialSamplingModelParser.ORIGIN,
                        BirthDeathSerialSamplingModelParser.BDSS + "." + BirthDeathSerialSamplingModelParser.ORIGIN, prior, writer);

                writer.writeCloseTag(BirthDeathSerialSamplingModelParser.BIRTH_DEATH_SERIAL_MODEL);

                break;

            case BIRTH_DEATH_BASIC_REPRODUCTIVE_NUMBER:
                writer.writeComment(BirthDeathSerialSamplingModelParser.getCitationRT());

                writer.writeOpenTag(
                        BirthDeathEpidemiologyModelParser.BIRTH_DEATH_EPIDEMIOLOGY,
                        new Attribute[]{
                                new Attribute.Default<String>(XMLParser.ID, prefix + BirthDeathEpidemiologyModelParser.BIRTH_DEATH_EPIDEMIOLOGY),
                                new Attribute.Default<String>("units", Units.Utils.getDefaultUnitName(units))
                        }
                );

                writeParameter(BirthDeathEpidemiologyModelParser.R0,
                        BirthDeathEpidemiologyModelParser.R0, prior, writer);
                writeParameter(BirthDeathEpidemiologyModelParser.RECOVERY_RATE,
                        BirthDeathEpidemiologyModelParser.RECOVERY_RATE, prior, writer);
                writeParameter(BirthDeathEpidemiologyModelParser.SAMPLING_PROBABILITY,
                        BirthDeathEpidemiologyModelParser.SAMPLING_PROBABILITY, prior, writer);
                writeParameter(BirthDeathEpidemiologyModelParser.ORIGIN,
                        BirthDeathEpidemiologyModelParser.ORIGIN, prior, writer);

                writer.writeCloseTag(BirthDeathEpidemiologyModelParser.BIRTH_DEATH_EPIDEMIOLOGY);

                break;

            case SPECIES_BIRTH_DEATH:
            case SPECIES_YULE:
            case SPECIES_YULE_CALIBRATION:
                writer.writeComment("A prior assumption that the population size has remained constant");
                writer.writeComment("throughout the time spanned by the genealogy.");
                if (nodeHeightPrior == TreePriorType.SPECIES_YULE_CALIBRATION)
                    writer.writeComment("Calibrated Yule: Heled J, Drummond AJ (2011), Syst Biol, doi: " +
                            "10.1093/sysbio/syr087");

                writer.writeOpenTag(
                        ConstantPopulationModelParser.CONSTANT_POPULATION_MODEL,
                        new Attribute[]{
                                new Attribute.Default<String>(XMLParser.ID, prefix + "constant"),
                                new Attribute.Default<String>("units", Units.Utils.getDefaultUnitName(options.units))
                        }
                );

                // initial value for pop mean is the same as what used to be the value for the population size
                Parameter para = options.starBEASTOptions.getParameter(TraitData.TRAIT_SPECIES + "." + options.starBEASTOptions.POP_MEAN);
                prior.getParameter("constant.popSize").setInitial(para.getInitial());

                writer.writeOpenTag(ConstantPopulationModelParser.POPULATION_SIZE);
                writeParameter("constant.popSize", prior, writer);
                writer.writeCloseTag(ConstantPopulationModelParser.POPULATION_SIZE);
                writer.writeCloseTag(ConstantPopulationModelParser.CONSTANT_POPULATION_MODEL);

                break;
        }

        if ((!options.useStarBEAST) && nodeHeightPrior != TreePriorType.CONSTANT && nodeHeightPrior != TreePriorType.EXPONENTIAL) {
            // If the node height prior is not one of these two then we need to simulate a
            // random starting tree under a constant size coalescent.
            writer.writeComment("This is a simple constant population size coalescent model",
                    "that is used to generate an initial tree for the chain.");
            writer.writeOpenTag(
                    ConstantPopulationModelParser.CONSTANT_POPULATION_MODEL,
                    new Attribute[]{
                            new Attribute.Default<String>(XMLParser.ID, prefix + "initialDemo"),
                            new Attribute.Default<String>("units", Units.Utils.getDefaultUnitName(units))
                    }
            );

            writer.writeOpenTag(ConstantPopulationModelParser.POPULATION_SIZE);
            if (initialPopSize != null) {
                writer.writeIDref(ParameterParser.PARAMETER, prefix + initialPopSize);
            } else {
                writeParameter(prefix + "initialDemo.popSize", 1, 100.0, Double.NaN, Double.NaN, writer);
            }
            writer.writeCloseTag(ConstantPopulationModelParser.POPULATION_SIZE);
            writer.writeCloseTag(ConstantPopulationModelParser.CONSTANT_POPULATION_MODEL);
        }

//        if (nodeHeightPrior == TreePriorType.BIRTH_DEATH_BASIC_REPRODUCTIVE_NUMBER) {
//            writer.writeComment("R0 = b/(b*d+s*r)");
//            writer.writeOpenTag(RPNcalculatorStatisticParser.RPN_STATISTIC,
//                    new Attribute[]{
//                            new Attribute.Default<String>(XMLParser.ID, modelPrefix + "R0")
//                    });
//
//            writer.writeOpenTag(RPNcalculatorStatisticParser.VARIABLE,
//                    new Attribute[]{
//                            new Attribute.Default<String>(Statistic.NAME, modelPrefix + "b")
//                    });
//            writeParameterRef(modelPrefix + BirthDeathSerialSamplingModelParser.BDSS + "." + BirthDeathSerialSamplingModelParser.LAMBDA, writer);
//            writer.writeCloseTag(RPNcalculatorStatisticParser.VARIABLE);
//
//            writer.writeOpenTag(RPNcalculatorStatisticParser.VARIABLE,
//                    new Attribute[]{
//                            new Attribute.Default<String>(Statistic.NAME, modelPrefix + "d")
//                    });
//            writeParameterRef(modelPrefix + BirthDeathSerialSamplingModelParser.BDSS + "." + BirthDeathSerialSamplingModelParser.RELATIVE_MU, writer);
//            writer.writeCloseTag(RPNcalculatorStatisticParser.VARIABLE);
//
//            writer.writeOpenTag(RPNcalculatorStatisticParser.VARIABLE,
//                    new Attribute[]{
//                            new Attribute.Default<String>(Statistic.NAME, modelPrefix + "s")
//                    });
//            writeParameterRef(modelPrefix + BirthDeathSerialSamplingModelParser.BDSS + "." + BirthDeathSerialSamplingModelParser.PSI, writer);
//            writer.writeCloseTag(RPNcalculatorStatisticParser.VARIABLE);
//
//            writer.writeOpenTag(RPNcalculatorStatisticParser.VARIABLE,
//                    new Attribute[]{
//                            new Attribute.Default<String>(Statistic.NAME, modelPrefix + "r")
//                    });
//            writeParameterRef(modelPrefix + BirthDeathSerialSamplingModelParser.BDSS + "." + BirthDeathSerialSamplingModelParser.R, writer);
//            writer.writeCloseTag(RPNcalculatorStatisticParser.VARIABLE);
//
//            writer.writeOpenTag(RPNcalculatorStatisticParser.EXPRESSION,
//                    new Attribute[]{
//                            new Attribute.Default<String>(Statistic.NAME, modelPrefix + "R0")
//                    });
//            writer.writeText(modelPrefix + "b " + modelPrefix + "b " + modelPrefix + "d " + "* " + modelPrefix + "s " + modelPrefix + "r " + "* + /");
//            writer.writeCloseTag(RPNcalculatorStatisticParser.EXPRESSION);
//
//            writer.writeCloseTag(RPNcalculatorStatisticParser.RPN_STATISTIC);
//        }
    }

    /**
     * Write the prior on node heights (coalescent or speciational models)
     *
     * @param model  PartitionTreeModel
     * @param writer the writer
     */
    void writePriorLikelihood(PartitionTreeModel model, XMLWriter writer) {

        //tree model prefix
        String prefix = model.getPrefix();

        PartitionTreePrior prior = model.getPartitionTreePrior();
        TreePriorType treePrior = prior.getNodeHeightPrior();

//        String priorPrefix = prior.getPrefix();

        switch (treePrior) {
            case YULE:
            case BIRTH_DEATH:
            case BIRTH_DEATH_INCOMPLETE_SAMPLING:
            case BIRTH_DEATH_SERIAL_SAMPLING:
            case BIRTH_DEATH_BASIC_REPRODUCTIVE_NUMBER:
            case YULE_CALIBRATION:
                // generate a speciational process
                writer.writeComment("Generate a speciation likelihood for Yule or Birth Death");
                writer.writeOpenTag(
                        SpeciationLikelihoodParser.SPECIATION_LIKELIHOOD,
                        new Attribute[]{
                                new Attribute.Default<String>(XMLParser.ID, prefix + "speciation")
                        }
                );

                // write pop size socket
                writer.writeOpenTag(SpeciationLikelihoodParser.MODEL);
                writeNodeHeightPriorModelRef(prior, writer);
                writer.writeCloseTag(SpeciationLikelihoodParser.MODEL);
                writer.writeOpenTag(SpeciationLikelihoodParser.TREE);
                writer.writeIDref(TreeModel.TREE_MODEL, prefix + TreeModel.TREE_MODEL);
                writer.writeCloseTag(SpeciationLikelihoodParser.TREE);

                if (treePrior == TreePriorType.YULE_CALIBRATION) {

                    if (options.treeModelOptions.isNodeCalibrated(model) == 0) {
                        writer.writeOpenTag(SpeciationLikelihoodParser.CALIBRATION,
                                new Attribute[]{
                                        new Attribute.Default<String>(SpeciationLikelihoodParser.CORRECTION, prior.getCalibCorrectionType().toString())
                                });
                        writer.writeOpenTag(SpeciationLikelihoodParser.POINT);

                        String taxaId;
                        if (options.hasIdenticalTaxa()) {
                            taxaId = TaxaParser.TAXA;
                        } else {
                            taxaId = options.getDataPartitions(model).get(0).getPrefix() + TaxaParser.TAXA;
                        }
                        writer.writeIDref(TaxaParser.TAXA, taxaId);

                        writeDistribution(model.getParameter("treeModel.rootHeight"), true, writer);

                        writer.writeCloseTag(SpeciationLikelihoodParser.POINT);
                        writer.writeCloseTag(SpeciationLikelihoodParser.CALIBRATION);

                    } else if (options.treeModelOptions.isNodeCalibrated(model) == 1) {
                        // should be only 1 calibrated internal node with monophyletic for each tree at moment
                        Taxa t = (Taxa) options.getKeysFromValue(options.taxonSetsTreeModel, model).get(0);
                        Parameter nodeCalib = options.getStatistic(t);

                        writer.writeOpenTag(SpeciationLikelihoodParser.CALIBRATION,
                                new Attribute[]{
                                        new Attribute.Default<String>(SpeciationLikelihoodParser.CORRECTION, prior.getCalibCorrectionType().toString())
                                });
                        writer.writeOpenTag(SpeciationLikelihoodParser.POINT);

                        writer.writeIDref(TaxaParser.TAXA, t.getId());
                        writeDistribution(nodeCalib, true, writer);

                        writer.writeCloseTag(SpeciationLikelihoodParser.POINT);
                        writer.writeCloseTag(SpeciationLikelihoodParser.CALIBRATION);

                        if (!options.treeModelOptions.isNodeCalibrated(nodeCalib)) {
                            throw new IllegalArgumentException("Calibrated Yule model requires a calibration to be specified for node, " +
                                    nodeCalib.getName() + ".");
                        }
                    }
                }

                writer.writeCloseTag(SpeciationLikelihoodParser.SPECIATION_LIKELIHOOD);

                break;

            // AR - removed this special case for Logistic - it causes terrible problems.
            // Need to put informative priors on Logistic parameters.
//    		case LOGISTIC:
//    			writer.writeComment("Generate a boolean likelihood for Coalescent: Logistic Growth");
//    			writer.writeOpenTag(
//    	                BooleanLikelihoodParser.BOOLEAN_LIKELIHOOD,
//    	                new Attribute[]{new Attribute.Default<String>(XMLParser.ID, modelPrefix + "booleanLikelihood1")}
//    	        );
//    	        writer.writeOpenTag(
//    	                TestStatisticParser.TEST_STATISTIC,
//    	                new Attribute[]{
//    	                        new Attribute.Default<String>(XMLParser.ID, "test1"),
//    	                        new Attribute.Default<String>(Statistic.NAME, "test1")
//    	                }
//    	        );
//    	        writer.writeIDref(ParameterParser.PARAMETER, priorPrefix + "logistic.t50");
//    	        writer.writeOpenTag("lessThan");
//    	        writer.writeIDref(ParameterParser.PARAMETER, modelPrefix + "treeModel.rootHeight");
//    	        writer.writeCloseTag("lessThan");
//    	        writer.writeCloseTag(TestStatisticParser.TEST_STATISTIC);
//    	        writer.writeCloseTag(BooleanLikelihoodParser.BOOLEAN_LIKELIHOOD);
//
//    	        writer.writeOpenTag(
//	                    CoalescentLikelihoodParser.COALESCENT_LIKELIHOOD,
//	                    new Attribute[]{new Attribute.Default<String>(XMLParser.ID, modelPrefix + COALESCENT)}
//	            );
//	            writer.writeOpenTag(CoalescentLikelihoodParser.MODEL);
//	            writeNodeHeightPriorModelRef(prior, writer);
//	            writer.writeCloseTag(CoalescentLikelihoodParser.MODEL);
//	            writer.writeOpenTag(CoalescentLikelihoodParser.POPULATION_TREE);
//	            writer.writeIDref(TreeModel.TREE_MODEL, modelPrefix + TreeModel.TREE_MODEL);
//	            writer.writeCloseTag(CoalescentLikelihoodParser.POPULATION_TREE);
//	            writer.writeCloseTag(CoalescentLikelihoodParser.COALESCENT_LIKELIHOOD);
//
//    	        break;

            case SKYLINE:
                // generate a Bayesian skyline plot
                writer.writeComment("Generate a generalizedSkyLineLikelihood for Bayesian Skyline");
                writer.writeOpenTag(
                        BayesianSkylineLikelihoodParser.SKYLINE_LIKELIHOOD,
                        new Attribute[]{
                                new Attribute.Default<String>(XMLParser.ID, prefix + "skyline"),
                                new Attribute.Default<String>("linear",
                                        prior.getSkylineModel() == TreePriorParameterizationType.LINEAR_SKYLINE ? "true" : "false")
                        }
                );

                // write pop size socket
                writer.writeOpenTag(BayesianSkylineLikelihoodParser.POPULATION_SIZES);
                if (prior.getSkylineModel() == TreePriorParameterizationType.LINEAR_SKYLINE) {
                    writeParameter(prior.getParameter("skyline.popSize"), prior.getSkylineGroupCount() + 1, writer);
                } else {
                    writeParameter(prior.getParameter("skyline.popSize"), prior.getSkylineGroupCount(), writer);
                }
                writer.writeCloseTag(BayesianSkylineLikelihoodParser.POPULATION_SIZES);

                // write group size socket
                writer.writeOpenTag(BayesianSkylineLikelihoodParser.GROUP_SIZES);
                writeParameter(prior.getParameter("skyline.groupSize"), prior.getSkylineGroupCount(), writer);
                writer.writeCloseTag(BayesianSkylineLikelihoodParser.GROUP_SIZES);

                writer.writeOpenTag(CoalescentLikelihoodParser.POPULATION_TREE);
                writer.writeIDref(TreeModel.TREE_MODEL, prefix + TreeModel.TREE_MODEL);
                writer.writeCloseTag(CoalescentLikelihoodParser.POPULATION_TREE);

                writer.writeCloseTag(BayesianSkylineLikelihoodParser.SKYLINE_LIKELIHOOD);

                writer.writeText("");
                writeExponentialMarkovLikelihood(prior, writer);

                break;

            case EXTENDED_SKYLINE:
                // different format
                break;

            case GMRF_SKYRIDE:
                writer.writeComment("Generate a gmrfSkyrideLikelihood for GMRF Bayesian Skyride process");
                writer.writeOpenTag(
                        GMRFSkyrideLikelihoodParser.SKYLINE_LIKELIHOOD,
                        new Attribute[]{
                                new Attribute.Default<String>(XMLParser.ID, prefix + "skyride"),
                                new Attribute.Default<String>(GMRFSkyrideLikelihoodParser.TIME_AWARE_SMOOTHING,
                                        prior.getSkyrideSmoothing() == TreePriorParameterizationType.TIME_AWARE_SKYRIDE ? "true" : "false"),
                                new Attribute.Default<String>(GMRFSkyrideLikelihoodParser.RANDOMIZE_TREE,
                                        //TODO For GMRF, tree model/tree prior combination not implemented by BEAST yet. The validation is in BeastGenerator.checkOptions()
                                        options.getPartitionTreeModels(prior).get(0).getStartingTreeType() == StartingTreeType.UPGMA ? "true" : "false"),
                        }
                );

                int skyrideIntervalCount = options.taxonList.getTaxonCount() - 1;
                writer.writeOpenTag(GMRFSkyrideLikelihoodParser.POPULATION_PARAMETER);
                writer.writeComment("skyride.logPopSize is in log units unlike other popSize");
                writeParameter(prior.getParameter("skyride.logPopSize"), skyrideIntervalCount, writer);
                writer.writeCloseTag(GMRFSkyrideLikelihoodParser.POPULATION_PARAMETER);

                writer.writeOpenTag(GMRFSkyrideLikelihoodParser.GROUP_SIZES);
                writeParameter(prior.getParameter("skyride.groupSize"), skyrideIntervalCount, writer);
                writer.writeCloseTag(GMRFSkyrideLikelihoodParser.GROUP_SIZES);

                writer.writeOpenTag(GMRFSkyrideLikelihoodParser.PRECISION_PARAMETER);
                writeParameter(prior.getParameter("skyride.precision"), 1, writer);
                writer.writeCloseTag(GMRFSkyrideLikelihoodParser.PRECISION_PARAMETER);

                writer.writeOpenTag(GMRFSkyrideLikelihoodParser.POPULATION_TREE);
                writer.writeIDref(TreeModel.TREE_MODEL, prefix + TreeModel.TREE_MODEL);
                writer.writeCloseTag(GMRFSkyrideLikelihoodParser.POPULATION_TREE);

                writer.writeCloseTag(GMRFSkyrideLikelihoodParser.SKYLINE_LIKELIHOOD);

                break;

            case SKYGRID:
                break;
            case SPECIES_YULE:
            case SPECIES_YULE_CALIBRATION:
            case SPECIES_BIRTH_DEATH:
                break;

            default:
                // generate a coalescent process
                writer.writeComment("Generate a coalescent likelihood");
                writer.writeOpenTag(
                        CoalescentLikelihoodParser.COALESCENT_LIKELIHOOD,
                        new Attribute[]{new Attribute.Default<String>(XMLParser.ID, prefix + COALESCENT)}
                );
                writer.writeOpenTag(CoalescentLikelihoodParser.MODEL);
                writeNodeHeightPriorModelRef(prior, writer);
                writer.writeCloseTag(CoalescentLikelihoodParser.MODEL);
                writer.writeOpenTag(CoalescentLikelihoodParser.POPULATION_TREE);
                writer.writeIDref(TreeModel.TREE_MODEL, prefix + TreeModel.TREE_MODEL);
                writer.writeCloseTag(CoalescentLikelihoodParser.POPULATION_TREE);
                writer.writeCloseTag(CoalescentLikelihoodParser.COALESCENT_LIKELIHOOD);
        }
    }

    void writeNodeHeightPriorModelRef(PartitionTreePrior prior, XMLWriter writer) {
        TreePriorType treePrior = prior.getNodeHeightPrior();
        String priorPrefix = prior.getPrefix();

        switch (treePrior) {
            case CONSTANT:
            case SPECIES_YULE:
            case SPECIES_BIRTH_DEATH:
            case SPECIES_YULE_CALIBRATION:
                writer.writeIDref(ConstantPopulationModelParser.CONSTANT_POPULATION_MODEL, priorPrefix + "constant");
                break;
            case EXPONENTIAL:
                writer.writeIDref(ExponentialGrowthModelParser.EXPONENTIAL_GROWTH_MODEL, priorPrefix + "exponential");
                break;
            case LOGISTIC:
                writer.writeIDref(LogisticGrowthModelParser.LOGISTIC_GROWTH_MODEL, priorPrefix + "logistic");
                break;
            case EXPANSION:
                writer.writeIDref(ExpansionModelParser.EXPANSION_MODEL, priorPrefix + "expansion");
                break;
            case SKYLINE:
                writer.writeIDref(BayesianSkylineLikelihoodParser.SKYLINE_LIKELIHOOD, priorPrefix + "skyline");
                break;
            case GMRF_SKYRIDE:
                writer.writeIDref(GMRFSkyrideLikelihoodParser.SKYLINE_LIKELIHOOD, priorPrefix + "skyride");
                break;
            case SKYGRID:
                writer.writeIDref(GMRFSkyrideLikelihoodParser.SKYGRID_LIKELIHOOD, priorPrefix + "skygrid");
                break;
            case YULE:
            case YULE_CALIBRATION:
                writer.writeIDref(YuleModelParser.YULE_MODEL, priorPrefix + YuleModelParser.YULE);
                break;
            case BIRTH_DEATH:
            case BIRTH_DEATH_INCOMPLETE_SAMPLING:
                writer.writeIDref(BirthDeathModelParser.BIRTH_DEATH_MODEL, priorPrefix + BirthDeathModelParser.BIRTH_DEATH);
                break;
            case BIRTH_DEATH_SERIAL_SAMPLING:
                writer.writeIDref(BirthDeathSerialSamplingModelParser.BIRTH_DEATH_SERIAL_MODEL,
                        priorPrefix + BirthDeathSerialSamplingModelParser.BDSS);
                break;
            case BIRTH_DEATH_BASIC_REPRODUCTIVE_NUMBER:
                writer.writeIDref(BirthDeathEpidemiologyModelParser.BIRTH_DEATH_EPIDEMIOLOGY,
                        priorPrefix + BirthDeathEpidemiologyModelParser.BIRTH_DEATH_EPIDEMIOLOGY);
                break;
            default:
                throw new IllegalArgumentException("No tree prior has been specified so cannot refer to it");
        }
    }

    void writeMultiLociTreePriors(PartitionTreePrior prior, XMLWriter writer) {

        String priorPrefix = prior.getPrefix();

        if (prior.getNodeHeightPrior() == TreePriorType.SKYGRID) {

            writer.writeComment("Generate a gmrfSkyGridLikelihood for the Bayesian SkyGrid process");
            writer.writeOpenTag(
                    GMRFSkyrideLikelihoodParser.SKYGRID_LIKELIHOOD,
                    new Attribute[]{
                            new Attribute.Default<String>(XMLParser.ID, priorPrefix + "skygrid"),
                    }
            );

            int skyGridIntervalCount = prior.getSkyGridCount();
            double skyGridInterval = prior.getSkyGridInterval();

            writer.writeOpenTag(GMRFSkyrideLikelihoodParser.POPULATION_PARAMETER);
            writer.writeComment("skygrid.logPopSize is in log units unlike other popSize");
            writeParameter(prior.getParameter("skygrid.logPopSize"), skyGridIntervalCount, writer);
            writer.writeCloseTag(GMRFSkyrideLikelihoodParser.POPULATION_PARAMETER);

            writer.writeOpenTag(GMRFSkyrideLikelihoodParser.PRECISION_PARAMETER);
            writeParameter(prior.getParameter("skygrid.precision"), 1, writer);
            writer.writeCloseTag(GMRFSkyrideLikelihoodParser.PRECISION_PARAMETER);

            writer.writeOpenTag(GMRFSkyrideLikelihoodParser.NUM_GRID_POINTS);
            Parameter numGridPoint = prior.getParameter("skygrid.numGridPoints");
            numGridPoint.setInitial(skyGridIntervalCount - 1);
            writeParameter(numGridPoint, 1, writer);
            writer.writeCloseTag(GMRFSkyrideLikelihoodParser.NUM_GRID_POINTS);

            writer.writeOpenTag(GMRFSkyrideLikelihoodParser.CUT_OFF);
            Parameter cutOff = prior.getParameter("skygrid.cutOff");
            cutOff.setInitial(skyGridInterval);
            writeParameter(cutOff, 1, writer);
            writer.writeCloseTag(GMRFSkyrideLikelihoodParser.CUT_OFF);

            writer.writeOpenTag(GMRFSkyrideLikelihoodParser.POPULATION_TREE);
            // TODO Add all linked trees
            if (options.isShareSameTreePrior()) {
                for (PartitionTreeModel thisModel : options.getPartitionTreeModels()) {
                    writer.writeIDref(TreeModel.TREE_MODEL, thisModel.getPrefix() + TreeModel.TREE_MODEL);
                }
            } else {
                writer.writeIDref(TreeModel.TREE_MODEL, options.getPartitionTreeModels(prior).get(0).getPrefix() + TreeModel.TREE_MODEL);
            }
            writer.writeCloseTag(GMRFSkyrideLikelihoodParser.POPULATION_TREE);

            writer.writeCloseTag(GMRFSkyrideLikelihoodParser.SKYGRID_LIKELIHOOD);

        } else if (prior.getNodeHeightPrior() == TreePriorType.EXTENDED_SKYLINE) {

            final String tagName = VariableDemographicModelParser.MODEL_NAME;
            writer.writeComment("Generate a variableDemographic for extended Bayesian skyline process");
            writer.writeOpenTag(tagName, new Attribute[]{
                    new Attribute.Default<String>(XMLParser.ID, priorPrefix + VariableDemographicModelParser.demoElementName),
                    new Attribute.Default<String>(VariableDemographicModelParser.TYPE, prior.getExtendedSkylineModel().toString()),
                    // use midpoint by default (todo) would be nice to have a user 'tickable' option
                    new Attribute.Default<String>(VariableDemographicModelParser.USE_MIDPOINTS, "true")
            }
            );

            Parameter popSize = prior.getParameter(VariableDemographicModelParser.demoElementName + ".popSize");
            Parameter populationMean = prior.getParameter(VariableDemographicModelParser.demoElementName + ".populationMean");
            popSize.setInitial(populationMean.getInitial());

            writer.writeOpenTag(VariableDemographicModelParser.POPULATION_SIZES);
            writer.writeComment("popSize value = populationMean value");
            writer.writeTag(ParameterParser.PARAMETER,
                    new Attribute[]{
                            new Attribute.Default<String>(XMLParser.ID, priorPrefix + VariableDemographicModelParser.demoElementName + ".popSize"),
                            new Attribute.Default<String>(ParameterParser.VALUE, Double.toString(popSize.getInitial()))}, true);
//	        writeParameter(popSize, -1, writer);
            writer.writeCloseTag(VariableDemographicModelParser.POPULATION_SIZES);

//            Parameter indicators = prior.getParameter(VariableDemographicModelParser.demoElementName + ".indicators");
            writer.writeOpenTag(VariableDemographicModelParser.INDICATOR_PARAMETER);
            writer.writeTag(ParameterParser.PARAMETER,
                    new Attribute[]{
                            new Attribute.Default<String>(XMLParser.ID, priorPrefix + VariableDemographicModelParser.demoElementName + ".indicators"),
                            new Attribute.Default<String>(ParameterParser.VALUE, Double.toString(0.0))}, true); // also 0.0
//	        writeParameter(prior.getParameter(VariableDemographicModelParser.demoElementName + ".indicators"), -1, writer); // not need dimension
            writer.writeCloseTag(VariableDemographicModelParser.INDICATOR_PARAMETER);

            writer.writeOpenTag(VariableDemographicModelParser.POPULATION_TREES);

            if (options.isShareSameTreePrior()) {
                for (PartitionTreeModel model : options.getPartitionTreeModels()) {
                    writer.writeOpenTag(VariableDemographicModelParser.POP_TREE, new Attribute[]{
                            new Attribute.Default<String>(SpeciesBindingsParser.PLOIDY, Double.toString(model.getPloidyType().getValue()))
                    }
                    );
                    writer.writeIDref(TreeModel.TREE_MODEL, model.getPrefix() + TreeModel.TREE_MODEL);
                    writer.writeCloseTag(VariableDemographicModelParser.POP_TREE);
                }
            } else {//TODO correct for not sharing same prior?
                writer.writeOpenTag(VariableDemographicModelParser.POP_TREE, new Attribute[]{
                        new Attribute.Default<String>(SpeciesBindingsParser.PLOIDY, Double.toString(options.getPartitionTreeModels(prior).get(0).getPloidyType().getValue()))
                }
                );
                writer.writeIDref(TreeModel.TREE_MODEL, options.getPartitionTreeModels(prior).get(0).getPrefix() + TreeModel.TREE_MODEL);
                writer.writeCloseTag(VariableDemographicModelParser.POP_TREE);
            }

            writer.writeCloseTag(VariableDemographicModelParser.POPULATION_TREES);

            writer.writeCloseTag(tagName);

            writer.writeOpenTag(CoalescentLikelihoodParser.COALESCENT_LIKELIHOOD, new Attribute.Default<String>(XMLParser.ID, priorPrefix + COALESCENT));
            writer.writeOpenTag(CoalescentLikelihoodParser.MODEL);
            writer.writeIDref(tagName, priorPrefix + VariableDemographicModelParser.demoElementName);
            writer.writeCloseTag(CoalescentLikelihoodParser.MODEL);
            writer.writeComment("Take population Tree from demographic");
            writer.writeCloseTag(CoalescentLikelihoodParser.COALESCENT_LIKELIHOOD);

            writer.writeOpenTag(SumStatisticParser.SUM_STATISTIC,
                    new Attribute[]{
                            new Attribute.Default<String>(XMLParser.ID, priorPrefix + VariableDemographicModelParser.demoElementName + ".populationSizeChanges"),
                            new Attribute.Default<String>("elementwise", "true")
                    });
            writer.writeIDref(ParameterParser.PARAMETER, priorPrefix + VariableDemographicModelParser.demoElementName + ".indicators");
            writer.writeCloseTag(SumStatisticParser.SUM_STATISTIC);
            writer.writeOpenTag(ExponentialDistributionModel.EXPONENTIAL_DISTRIBUTION_MODEL,
                    new Attribute[]{
                            new Attribute.Default<String>(XMLParser.ID, priorPrefix + VariableDemographicModelParser.demoElementName + ".populationMeanDist")
                            //,new Attribute.Default<String>("elementwise", "true")
                    });
            writer.writeOpenTag(DistributionModelParser.MEAN);

            writer.writeComment("prefer populationMean value = 1");
            populationMean = prior.getParameter(VariableDemographicModelParser.demoElementName + ".populationMean");
            writer.writeTag(ParameterParser.PARAMETER,
                    new Attribute[]{
                            new Attribute.Default<String>(XMLParser.ID, priorPrefix + VariableDemographicModelParser.demoElementName + ".populationMean"),
                            new Attribute.Default<String>(ParameterParser.VALUE, Double.toString(populationMean.getInitial()))}, true);

            writer.writeCloseTag(DistributionModelParser.MEAN);
            writer.writeCloseTag(ExponentialDistributionModel.EXPONENTIAL_DISTRIBUTION_MODEL);
        }
    }

    void writeParameterLog(PartitionTreePrior prior, XMLWriter writer) {

        String priorPrefix = prior.getPrefix();

        switch (prior.getNodeHeightPrior()) {

            case CONSTANT:
                writeParameterRef(priorPrefix + "constant.popSize", writer);
                break;
            case EXPONENTIAL:
                writeParameterRef(priorPrefix + "exponential.popSize", writer);
                if (prior.getParameterization() == TreePriorParameterizationType.GROWTH_RATE) {
                    writeParameterRef(priorPrefix + "exponential.growthRate", writer);
                } else {
                    writeParameterRef(priorPrefix + "exponential.doublingTime", writer);
                }
                break;
            case LOGISTIC:
                writeParameterRef(priorPrefix + "logistic.popSize", writer);
                if (prior.getParameterization() == TreePriorParameterizationType.GROWTH_RATE) {
                    writeParameterRef(priorPrefix + "logistic.growthRate", writer);
                } else {
                    writeParameterRef(priorPrefix + "logistic.doublingTime", writer);
                }
                writeParameterRef(priorPrefix + "logistic.t50", writer);
                break;
            case EXPANSION:
                writeParameterRef(priorPrefix + "expansion.popSize", writer);
                if (prior.getParameterization() == TreePriorParameterizationType.GROWTH_RATE) {
                    writeParameterRef(priorPrefix + "expansion.growthRate", writer);
                } else {
                    writeParameterRef(priorPrefix + "expansion.doublingTime", writer);
                }
                writeParameterRef(priorPrefix + "expansion.ancestralProportion", writer);
                break;
            case SKYLINE:
                writeParameterRef(priorPrefix + "skyline.popSize", writer);
                writeParameterRef(priorPrefix + "skyline.groupSize", writer);
                break;
            case EXTENDED_SKYLINE:
                writer.writeIDref(SumStatisticParser.SUM_STATISTIC, "demographic.populationSizeChanges");
                writeParameterRef(priorPrefix + "demographic.populationMean", writer);
                writeParameterRef(priorPrefix + "demographic.popSize", writer);
                writeParameterRef(priorPrefix + "demographic.indicators", writer);
                break;
            case SKYGRID:
                writeParameterRef(priorPrefix + "skygrid.precision", writer);
                writeParameterRef(priorPrefix + "skygrid.logPopSize", writer);
                writeParameterRef(priorPrefix + "skygrid.cutOff", writer);
//                writeParameterRef(priorPrefix + "skygrid.groupSize", writer);
                break;
            case GMRF_SKYRIDE:
                writeParameterRef(priorPrefix + "skyride.precision", writer);
                writeParameterRef(priorPrefix + "skyride.logPopSize", writer);
                writeParameterRef(priorPrefix + "skyride.groupSize", writer);
                break;
            case YULE:
            case YULE_CALIBRATION:
                writeParameterRef(priorPrefix + "yule.birthRate", writer);
                break;
            case BIRTH_DEATH:
            case BIRTH_DEATH_INCOMPLETE_SAMPLING:
                writeParameterRef(priorPrefix + BirthDeathModelParser.MEAN_GROWTH_RATE_PARAM_NAME, writer);
                writeParameterRef(priorPrefix + BirthDeathModelParser.RELATIVE_DEATH_RATE_PARAM_NAME, writer);
                if (prior.getNodeHeightPrior() == TreePriorType.BIRTH_DEATH_INCOMPLETE_SAMPLING)
                    writeParameterRef(priorPrefix + BirthDeathModelParser.BIRTH_DEATH + "."
                            + BirthDeathModelParser.SAMPLE_PROB, writer);
                break;
            case BIRTH_DEATH_SERIAL_SAMPLING:
                writeParameterRef(priorPrefix + BirthDeathSerialSamplingModelParser.BDSS + "."
                        + BirthDeathSerialSamplingModelParser.LAMBDA, writer);
                writeParameterRef(priorPrefix + BirthDeathSerialSamplingModelParser.BDSS + "."
                        + BirthDeathSerialSamplingModelParser.RELATIVE_MU, writer);
                //Issue 656
//                writeParameterRef(modelPrefix + BirthDeathSerialSamplingModelParser.BDSS + "."
//                        + BirthDeathSerialSamplingModelParser.SAMPLE_PROBABILITY, writer);
                writeParameterRef(priorPrefix + BirthDeathSerialSamplingModelParser.BDSS + "."
                        + BirthDeathSerialSamplingModelParser.PSI, writer);
                writeParameterRef(priorPrefix + BirthDeathSerialSamplingModelParser.BDSS + "."
                        + BirthDeathSerialSamplingModelParser.ORIGIN, writer);
                break;
            case BIRTH_DEATH_BASIC_REPRODUCTIVE_NUMBER:
                writeParameterRef(priorPrefix + BirthDeathEpidemiologyModelParser.R0, writer);
                writeParameterRef(priorPrefix + BirthDeathEpidemiologyModelParser.RECOVERY_RATE, writer);
                writeParameterRef(priorPrefix + BirthDeathEpidemiologyModelParser.SAMPLING_PROBABILITY, writer);
                writeParameterRef(priorPrefix + BirthDeathEpidemiologyModelParser.ORIGIN, writer);
            case SPECIES_YULE:
            case SPECIES_BIRTH_DEATH:
            case SPECIES_YULE_CALIBRATION:
                break;
            default:
                throw new IllegalArgumentException("No tree prior has been specified so cannot refer to it");
        }

    }

    void writeEBSPAnalysisToCSVfile(PartitionTreePrior prior, XMLWriter writer) {

        String priorPrefix = prior.getPrefix();

        String logFileName = options.logFileName;

        if (prior.getNodeHeightPrior() == TreePriorType.EXTENDED_SKYLINE) {
            writer.writeOpenTag(EBSPAnalysisParser.VD_ANALYSIS, new Attribute[]{
                    new Attribute.Default<String>(XMLParser.ID, priorPrefix + "demographic.analysis"),
                    new Attribute.Default<Double>(EBSPAnalysisParser.BURN_IN, 0.1),
                    new Attribute.Default<Boolean>(VariableDemographicModelParser.USE_MIDPOINTS, true)}
            );

            writer.writeOpenTag(EBSPAnalysisParser.LOG_FILE_NAME);
            writer.writeText(logFileName);
            writer.writeCloseTag(EBSPAnalysisParser.LOG_FILE_NAME);

            writer.writeOpenTag(EBSPAnalysisParser.TREE_FILE_NAMES);
            for (String treeFN : options.treeFileName) {
                writer.writeOpenTag(EBSPAnalysisParser.TREE_LOG);
                writer.writeText(treeFN);
                writer.writeCloseTag(EBSPAnalysisParser.TREE_LOG);
            }
            writer.writeCloseTag(EBSPAnalysisParser.TREE_FILE_NAMES);

            writer.writeOpenTag(EBSPAnalysisParser.MODEL_TYPE);
            writer.writeText(prior.getExtendedSkylineModel().toString());
            writer.writeCloseTag(EBSPAnalysisParser.MODEL_TYPE);

            writer.writeOpenTag(EBSPAnalysisParser.POPULATION_FIRST_COLUMN);
            writer.writeText(VariableDemographicModelParser.demoElementName + ".popSize1");
            writer.writeCloseTag(EBSPAnalysisParser.POPULATION_FIRST_COLUMN);

            writer.writeOpenTag(EBSPAnalysisParser.INDICATORS_FIRST_COLUMN);
            writer.writeText(VariableDemographicModelParser.demoElementName + ".indicators1");
            writer.writeCloseTag(EBSPAnalysisParser.INDICATORS_FIRST_COLUMN);

            writer.writeCloseTag(EBSPAnalysisParser.VD_ANALYSIS);

            writer.writeOpenTag(CSVExporterParser.CSV_EXPORT,
                    new Attribute[]{
                            new Attribute.Default<String>(CSVExporterParser.FILE_NAME,
                                    logFileName.subSequence(0, logFileName.length() - 4) + ".csv"), //.log
                            new Attribute.Default<String>(CSVExporterParser.SEPARATOR, ",")
                    });
            writer.writeOpenTag(CSVExporterParser.COLUMNS);
            writer.writeIDref(EBSPAnalysisParser.VD_ANALYSIS, priorPrefix + "demographic.analysis");
            writer.writeCloseTag(CSVExporterParser.COLUMNS);
            writer.writeCloseTag(CSVExporterParser.CSV_EXPORT);
        }
    }

    private void writeExponentialMarkovLikelihood(PartitionTreePrior prior, XMLWriter writer) {

        writer.writeOpenTag(
                ExponentialMarkovModel.EXPONENTIAL_MARKOV_MODEL,
                new Attribute[]{new Attribute.Default<String>(XMLParser.ID, prior.getPrefix() + "eml1"),
                        new Attribute.Default<String>("jeffreys", "true")}
        );

        writeParameterRef(ExponentialMarkovModelParser.CHAIN_PARAMETER, prior.getPrefix() + "skyline.popSize", writer);

        writer.writeCloseTag(ExponentialMarkovModel.EXPONENTIAL_MARKOV_MODEL);
    }

    public static void writePriorLikelihoodReferenceLog(PartitionTreePrior prior, PartitionTreeModel model, XMLWriter writer) {

        String prefix = model.getPrefix();

        switch (prior.getNodeHeightPrior()) {

            case YULE:
            case YULE_CALIBRATION:
            case BIRTH_DEATH:
            case BIRTH_DEATH_INCOMPLETE_SAMPLING:
            case BIRTH_DEATH_SERIAL_SAMPLING:
            case BIRTH_DEATH_BASIC_REPRODUCTIVE_NUMBER:
                writer.writeIDref(SpeciationLikelihoodParser.SPECIATION_LIKELIHOOD, prefix + "speciation");
                break;
            case SKYLINE:
                writer.writeIDref(BayesianSkylineLikelihoodParser.SKYLINE_LIKELIHOOD, prefix + "skyline");
//                writer.writeIDref(ExponentialMarkovModel.EXPONENTIAL_MARKOV_MODEL, modelPrefix + "eml1");
                break;
            case GMRF_SKYRIDE:
                writer.writeIDref(GMRFSkyrideLikelihoodParser.SKYLINE_LIKELIHOOD, prefix + "skyride");
                break;
            case SKYGRID:
                writer.writeIDref(GMRFSkyrideLikelihoodParser.SKYLINE_LIKELIHOOD, prefix + "skygrid");
                // only 1 coalescent, so write it separately after this method
                break;
            case LOGISTIC:
//                writer.writeIDref(BooleanLikelihoodParser.BOOLEAN_LIKELIHOOD, modelPrefix + "booleanLikelihood1");
                writer.writeIDref(CoalescentLikelihoodParser.COALESCENT_LIKELIHOOD, prefix + COALESCENT);
                break;
            case EXTENDED_SKYLINE:
                // only 1 coalescent, so write it separately after this method
            case SPECIES_YULE:
            case SPECIES_YULE_CALIBRATION:
            case SPECIES_BIRTH_DEATH:
                // do not need
                break;
            default:
                writer.writeIDref(CoalescentLikelihoodParser.COALESCENT_LIKELIHOOD, prefix + COALESCENT);
        }
    }

    // id is written in writePriorLikelihood (PartitionTreePrior prior, PartitionTreeModel model, XMLWriter writer)
    public void writePriorLikelihoodReference(PartitionTreePrior prior, PartitionTreeModel model, XMLWriter writer) {

        String prefix = model.getPrefix();

        switch (prior.getNodeHeightPrior()) {

            case YULE:
            case YULE_CALIBRATION:
            case BIRTH_DEATH:
            case BIRTH_DEATH_INCOMPLETE_SAMPLING:
            case BIRTH_DEATH_SERIAL_SAMPLING:
            case BIRTH_DEATH_BASIC_REPRODUCTIVE_NUMBER:
                writer.writeIDref(SpeciationLikelihoodParser.SPECIATION_LIKELIHOOD, prefix + "speciation");
                break;
            case SKYLINE:
                writer.writeIDref(BayesianSkylineLikelihoodParser.SKYLINE_LIKELIHOOD, prefix + "skyline");
                writer.writeIDref(ExponentialMarkovModel.EXPONENTIAL_MARKOV_MODEL, prefix + "eml1");
                break;
            case GMRF_SKYRIDE:
                writer.writeIDref(GMRFSkyrideLikelihoodParser.SKYLINE_LIKELIHOOD, prefix + "skyride");
                break;
            case SKYGRID:
//                writer.writeIDref(GMRFSkyrideLikelihoodParser.SKYLINE_LIKELIHOOD, prefix + "skygrid");
                // only 1 coalescent, so write it separately after this method
                break;
//            case LOGISTIC:
//                writer.writeIDref(BooleanLikelihoodParser.BOOLEAN_LIKELIHOOD, prefix + "booleanLikelihood1");
//                writer.writeIDref(CoalescentLikelihoodParser.COALESCENT_LIKELIHOOD, modelPrefix + COALESCENT);
//                break;
            case EXTENDED_SKYLINE:
                // only 1 coalescent, so write it in writeEBSPVariableDemographicReference
            case SPECIES_YULE:
            case SPECIES_YULE_CALIBRATION:
            case SPECIES_BIRTH_DEATH:
                // do not need
                break;
            default:
                writer.writeIDref(CoalescentLikelihoodParser.COALESCENT_LIKELIHOOD, prefix + COALESCENT);
        }
    }

    public void writeMultiLociLikelihoodReference(PartitionTreePrior prior, XMLWriter writer) {

        String prefix = prior.getPrefix();

        //TODO: make suitable for *BEAST
        if (prior.getNodeHeightPrior() == TreePriorType.EXTENDED_SKYLINE) {

            writer.writeIDref(CoalescentLikelihoodParser.COALESCENT_LIKELIHOOD, prefix + COALESCENT); // only 1 coalescent

            writer.writeOpenTag(MixedDistributionLikelihoodParser.DISTRIBUTION_LIKELIHOOD);

            writer.writeOpenTag(MixedDistributionLikelihoodParser.DISTRIBUTION0);
            writer.writeIDref(ExponentialDistributionModel.EXPONENTIAL_DISTRIBUTION_MODEL, prefix + "demographic.populationMeanDist");
            writer.writeCloseTag(MixedDistributionLikelihoodParser.DISTRIBUTION0);

            writer.writeOpenTag(MixedDistributionLikelihoodParser.DISTRIBUTION1);
            writer.writeIDref(ExponentialDistributionModel.EXPONENTIAL_DISTRIBUTION_MODEL, prefix + "demographic.populationMeanDist");
            writer.writeCloseTag(MixedDistributionLikelihoodParser.DISTRIBUTION1);

            writeParameterRef(MixedDistributionLikelihoodParser.DATA, prefix + "demographic.popSize", writer);

            writeParameterRef(MixedDistributionLikelihoodParser.INDICATORS, prefix + "demographic.indicators", writer);

            writer.writeCloseTag(MixedDistributionLikelihoodParser.DISTRIBUTION_LIKELIHOOD);
        } else if (prior.getNodeHeightPrior() == TreePriorType.SKYGRID) {
            writer.writeIDref(GMRFSkyrideLikelihoodParser.SKYGRID_LIKELIHOOD, prefix + "skygrid");
        }

    }
}
