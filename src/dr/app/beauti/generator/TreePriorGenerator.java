/*
 * TreePriorGenerator.java
 *
 * Copyright (C) 2002-2009 Alexei Drummond and Andrew Rambaut
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
 * BEAST is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package dr.app.beauti.generator;

import dr.app.beauti.components.ComponentFactory;
import dr.app.beauti.enumTypes.StartingTreeType;
import dr.app.beauti.enumTypes.TreePriorParameterizationType;
import dr.app.beauti.enumTypes.TreePriorType;
import dr.app.beauti.options.*;
import dr.app.beauti.util.XMLWriter;
import dr.evolution.util.Units;
import dr.evomodel.speciation.BirthDeathGernhard08Model;
import dr.evomodel.tree.TreeModel;
import dr.evomodelxml.CSVExporterParser;
import dr.evomodelxml.coalescent.*;
import dr.evomodelxml.speciation.BirthDeathModelParser;
import dr.evomodelxml.speciation.SpeciationLikelihoodParser;
import dr.evomodelxml.speciation.SpeciesBindingsParser;
import dr.evomodelxml.speciation.YuleModelParser;
import dr.inference.distribution.ExponentialDistributionModel;
import dr.inference.distribution.ExponentialMarkovModel;
import dr.inference.model.ParameterParser;
import dr.inferencexml.distribution.DistributionModelParser;
import dr.inferencexml.distribution.ExponentialMarkovModelParser;
import dr.inferencexml.distribution.MixedDistributionLikelihoodParser;
import dr.inferencexml.model.BooleanLikelihoodParser;
import dr.inferencexml.model.SumStatisticParser;
import dr.inferencexml.model.TestStatisticParser;
import dr.util.Attribute;
import dr.xml.XMLParser;

/**
 * @author Alexei Drummond
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

        setModelPrefix(prior.getPrefix()); // only has prefix, if (options.getPartitionTreePriors().size() > 1)

        String initialPopSize = null;

        TreePriorType nodeHeightPrior = prior.getNodeHeightPrior();
        Units.Type units = options.units;
        TreePriorParameterizationType parameterization = prior.getParameterization();

        switch (nodeHeightPrior) {
        	case CONSTANT:  
	            writer.writeComment("A prior assumption that the population size has remained constant");
	            writer.writeComment("throughout the time spanned by the genealogy.");
	            writer.writeOpenTag(
	                    ConstantPopulationModelParser.CONSTANT_POPULATION_MODEL,
	                    new Attribute[]{
	                            new Attribute.Default<String>(XMLParser.ID, modelPrefix + "constant"),
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
	
	            writer.writeComment("A prior assumption that the population size has grown exponentially");
	            writer.writeComment("throughout the time spanned by the genealogy.");
	            writer.writeOpenTag(
	                    ExponentialGrowthModelParser.EXPONENTIAL_GROWTH_MODEL,
	                    new Attribute[]{
	                            new Attribute.Default<String>(XMLParser.ID, modelPrefix + "exponential"),
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
	
	            writer.writeComment("A prior assumption that the population size has grown logistically");
	            writer.writeComment("throughout the time spanned by the genealogy.");
	            writer.writeOpenTag(
	                    LogisticGrowthModelParser.LOGISTIC_GROWTH_MODEL,
	                    new Attribute[]{
	                            new Attribute.Default<String>(XMLParser.ID, modelPrefix + "logistic"),
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
                writer.writeComment("logistic.t50 initial always has to < treeRootHeight initial");
                dr.app.beauti.options.Parameter priorPara = prior.getParameter("logistic.t50");

                double initRootHeight;
                if (options.isShareSameTreePrior()) {
                    initRootHeight = priorPara.initial;
                    for (PartitionTreeModel tree : options.getPartitionTreeModels()) {
                        double tmpRootHeight = tree.getParameter("treeModel.rootHeight").initial;
                        if (initRootHeight > tmpRootHeight) { // take min
                            initRootHeight = tmpRootHeight;
                        }
                    }
                } else {
                    initRootHeight = prior.getTreeModel().getParameter("treeModel.rootHeight").initial;
                }
                // logistic.t50 initial always has to < treeRootHeight initial
                if (priorPara.initial >= initRootHeight) {
                    priorPara.initial = initRootHeight / 2; // tree prior.initial has to < treeRootHeight.initial
                }
//	            } else {
//	            	writer.writeComment("Has calibration");
//	            	//TODO
//	            	throw new IllegalArgumentException("This function is not available in this release !");
//	            }
	            
	            writeParameter("logistic.t50", prior, writer);
	            writer.writeCloseTag(LogisticGrowthModelParser.TIME_50);
	
	            writer.writeCloseTag(LogisticGrowthModelParser.LOGISTIC_GROWTH_MODEL);
	
	            initialPopSize = "logistic.popSize";
	            
	            break;

        	case EXPANSION:
        		// generate an exponential prior tree
	            writer.writeComment("A prior assumption that the population size has grown exponentially");
	            writer.writeComment("from some ancestral population size in the past.");
	            writer.writeOpenTag(
	                    ExpansionModelParser.EXPANSION_MODEL,
	                    new Attribute[]{
	                            new Attribute.Default<String>(XMLParser.ID, modelPrefix + "expansion"),
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
	            writer.writeComment("A prior on the distribution node heights defined given");
	            writer.writeComment("a Yule speciation process (a pure birth process).");
	            writer.writeOpenTag(
	                    YuleModelParser.YULE_MODEL,
	                    new Attribute[]{
	                            new Attribute.Default<String>(XMLParser.ID, modelPrefix + "yule"),
	                            new Attribute.Default<String>("units", Units.Utils.getDefaultUnitName(units))
	                    }
	            );
	
	            writeParameter(YuleModelParser.BIRTH_RATE, "yule.birthRate", prior, writer);
	            writer.writeCloseTag(YuleModelParser.YULE_MODEL);
	            
	            break;
	            
        	case BIRTH_DEATH:
	            writer.writeComment("A prior on the distribution node heights defined given");
	            writer.writeComment("a Birth-Death speciation process (Gernhard 2008).");
	            writer.writeOpenTag(
	                    BirthDeathGernhard08Model.BIRTH_DEATH_MODEL,
	                    new Attribute[]{
	                            new Attribute.Default<String>(XMLParser.ID, modelPrefix + "birthDeath"),
	                            new Attribute.Default<String>("units", Units.Utils.getDefaultUnitName(units))
	                    }
	            );
	
	            writeParameter(BirthDeathModelParser.BIRTHDIFF_RATE, BirthDeathModelParser.MEAN_GROWTH_RATE_PARAM_NAME, prior, writer);
	            writeParameter(BirthDeathModelParser.RELATIVE_DEATH_RATE, BirthDeathModelParser.RELATIVE_DEATH_RATE_PARAM_NAME, prior, writer);
	            writer.writeCloseTag(BirthDeathGernhard08Model.BIRTH_DEATH_MODEL);
	            
	            break;
	            
        	case SPECIES_BIRTH_DEATH:
        	case SPECIES_YULE:
	            writer.writeComment("A prior assumption that the population size has remained constant");
	            writer.writeComment("throughout the time spanned by the genealogy.");
	            writer.writeOpenTag(
	                    ConstantPopulationModelParser.CONSTANT_POPULATION_MODEL,
	                    new Attribute[]{
	                            new Attribute.Default<String>(XMLParser.ID, modelPrefix + "constant"),
	                            new Attribute.Default<String>("units", Units.Utils.getDefaultUnitName(options.units))
	                    }
	            );
	            
	            // initial value for pop mean is the same as what used to be the value for the population size
	            Parameter para = options.starBEASTOptions.getParameter(TraitData.Traits.TRAIT_SPECIES + "." + options.starBEASTOptions.POP_MEAN);
	            prior.getParameter("constant.popSize").initial = para.initial;
	
	            writer.writeOpenTag(ConstantPopulationModelParser.POPULATION_SIZE);
	            writeParameter("constant.popSize", prior, writer);
	            writer.writeCloseTag(ConstantPopulationModelParser.POPULATION_SIZE);
	            writer.writeCloseTag(ConstantPopulationModelParser.CONSTANT_POPULATION_MODEL);
	            
	            break;	            
        }

        if ((!options.starBEASTOptions.isSpeciesAnalysis()) && nodeHeightPrior != TreePriorType.CONSTANT && nodeHeightPrior != TreePriorType.EXPONENTIAL) {
            // If the node height prior is not one of these two then we need to simulate a
            // random starting tree under a constant size coalescent.

            writer.writeComment("This is a simple constant population size coalescent model");
            writer.writeComment("that is used to generate an initial tree for the chain.");
            writer.writeOpenTag(
                    ConstantPopulationModelParser.CONSTANT_POPULATION_MODEL,
                    new Attribute[]{
                            new Attribute.Default<String>(XMLParser.ID, modelPrefix + "initialDemo"),
                            new Attribute.Default<String>("units", Units.Utils.getDefaultUnitName(units))
                    }
            );

            writer.writeOpenTag(ConstantPopulationModelParser.POPULATION_SIZE);
            if (initialPopSize != null) {
            	writer.writeIDref(ParameterParser.PARAMETER, modelPrefix + initialPopSize);
            } else {
                writeParameter(modelPrefix + "initialDemo.popSize", 1, 100.0, Double.NaN, Double.NaN, writer);
            }
            writer.writeCloseTag(ConstantPopulationModelParser.POPULATION_SIZE);
            writer.writeCloseTag(ConstantPopulationModelParser.CONSTANT_POPULATION_MODEL);
        }
    }

    /**
     * Write the prior on node heights (coalescent or speciational models)
     *
     * @param prior  the partition tree prior
     * @param model   PartitionTreeModel
     * @param writer the writer
     */
    void writePriorLikelihood (PartitionTreePrior prior, PartitionTreeModel model, XMLWriter writer) {
    	
    	//tree model prefix
    	setModelPrefix(model.getPrefix()); // only has prefix, if (options.getPartitionTreePriors().size() > 1)
    	String priorPrefix = prior.getPrefix();

        TreePriorType treePrior = prior.getNodeHeightPrior();
        
        switch (treePrior) {
    		case YULE:
    		case BIRTH_DEATH:
	            // generate a speciational process
	            writer.writeComment("Generate a speciation likelihood for Yule or Birth Death");
	            writer.writeOpenTag(
	                    SpeciationLikelihoodParser.SPECIATION_LIKELIHOOD,
	                    new Attribute[]{
	                            new Attribute.Default<String>(XMLParser.ID, modelPrefix + "speciation")
	                    }
	            );
	
	            // write pop size socket
	            writer.writeOpenTag(SpeciationLikelihoodParser.MODEL);
	            writeNodeHeightPriorModelRef(prior, writer);
	            writer.writeCloseTag(SpeciationLikelihoodParser.MODEL);
	            writer.writeOpenTag(SpeciationLikelihoodParser.TREE);
	            writer.writeIDref(TreeModel.TREE_MODEL, modelPrefix + TreeModel.TREE_MODEL);
	            writer.writeCloseTag(SpeciationLikelihoodParser.TREE);
	
	            writer.writeCloseTag(SpeciationLikelihoodParser.SPECIATION_LIKELIHOOD);
	            
	            break;
	            
    		case LOGISTIC:
    			writer.writeComment("Generate a boolean likelihood for Coalescent: Logistic Growth");
    			writer.writeOpenTag(
    	                BooleanLikelihoodParser.BOOLEAN_LIKELIHOOD,
    	                new Attribute[]{new Attribute.Default<String>(XMLParser.ID, modelPrefix + "booleanLikelihood1")}
    	        );
    	        writer.writeOpenTag(
    	                TestStatisticParser.TEST_STATISTIC,
    	                new Attribute[]{
    	                        new Attribute.Default<String>(XMLParser.ID, "test1"),
    	                        new Attribute.Default<String>("name", "test1")
    	                }
    	        );
    	        writer.writeIDref(ParameterParser.PARAMETER, priorPrefix + "logistic.t50"); //TODO correct?
    	        writer.writeOpenTag("lessThan");
    	        writer.writeIDref(ParameterParser.PARAMETER, modelPrefix + "treeModel.rootHeight");
    	        writer.writeCloseTag("lessThan");
    	        writer.writeCloseTag(TestStatisticParser.TEST_STATISTIC);
    	        writer.writeCloseTag(BooleanLikelihoodParser.BOOLEAN_LIKELIHOOD);
    	        
    	        writer.writeOpenTag(
	                    CoalescentLikelihoodParser.COALESCENT_LIKELIHOOD,
	                    new Attribute[]{new Attribute.Default<String>(XMLParser.ID, modelPrefix + COALESCENT)}
	            );
	            writer.writeOpenTag(CoalescentLikelihoodParser.MODEL);
	            writeNodeHeightPriorModelRef(prior, writer);
	            writer.writeCloseTag(CoalescentLikelihoodParser.MODEL);
	            writer.writeOpenTag(CoalescentLikelihoodParser.POPULATION_TREE);
	            writer.writeIDref(TreeModel.TREE_MODEL, modelPrefix + TreeModel.TREE_MODEL);
	            writer.writeCloseTag(CoalescentLikelihoodParser.POPULATION_TREE);
	            writer.writeCloseTag(CoalescentLikelihoodParser.COALESCENT_LIKELIHOOD);
	            
    	        break;
	            
    		case SKYLINE:
	            // generate a Bayesian skyline plot
	            writer.writeComment("Generate a generalizedSkyLineLikelihood for Bayesian Skyline");
	            writer.writeOpenTag(
	                    BayesianSkylineLikelihoodParser.SKYLINE_LIKELIHOOD,
	                    new Attribute[]{
	                            new Attribute.Default<String>(XMLParser.ID, modelPrefix + "skyline"),
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
	            writer.writeIDref(TreeModel.TREE_MODEL, modelPrefix + TreeModel.TREE_MODEL);
	            writer.writeCloseTag(CoalescentLikelihoodParser.POPULATION_TREE);
	
	            writer.writeCloseTag(BayesianSkylineLikelihoodParser.SKYLINE_LIKELIHOOD);
	            
	            writer.writeText("");
	            writeExponentialMarkovLikelihood(writer);
	            
	            break;
	            
    		case EXTENDED_SKYLINE:        	
	            // different format	            
	            break;
            
    		case GMRF_SKYRIDE:
	            writer.writeComment("Generate a gmrfSkyrideLikelihood for GMRF Bayesian Skyride process");
	            writer.writeOpenTag(
	                    GMRFSkyrideLikelihoodParser.SKYLINE_LIKELIHOOD,
	                    new Attribute[]{
	                            new Attribute.Default<String>(XMLParser.ID, modelPrefix + "skyride"),
	                            new Attribute.Default<String>(GMRFSkyrideLikelihoodParser.TIME_AWARE_SMOOTHING,
	                                    prior.getSkyrideSmoothing() == TreePriorParameterizationType.TIME_AWARE_SKYRIDE ? "true" : "false"),
	                            new Attribute.Default<String>(GMRFSkyrideLikelihoodParser.RANDOMIZE_TREE,
	                                    //TODO For GMRF, tree model/tree prior combination not implemented by BEAST yet. The validation is in BeastGenerator.checkOptions()
	                                    prior.getTreeModel().getStartingTreeType() == StartingTreeType.UPGMA ? "true" : "false"),
	                    }
	            );
	
	            int skyrideIntervalCount = BeautiOptions.taxonList.getTaxonCount() - 1;
	            writer.writeOpenTag(GMRFSkyrideLikelihoodParser.POPULATION_PARAMETER);
                writer.writeComment("skyride.logPopSize is log unit unlike other popSize");
	            writeParameter(prior.getParameter("skyride.logPopSize"), skyrideIntervalCount, writer);
	            writer.writeCloseTag(GMRFSkyrideLikelihoodParser.POPULATION_PARAMETER);
	
	            writer.writeOpenTag(GMRFSkyrideLikelihoodParser.GROUP_SIZES);
	            writeParameter(prior.getParameter("skyride.groupSize"), skyrideIntervalCount, writer);
	            writer.writeCloseTag(GMRFSkyrideLikelihoodParser.GROUP_SIZES);
	
	            writer.writeOpenTag(GMRFSkyrideLikelihoodParser.PRECISION_PARAMETER);
	            writeParameter(prior.getParameter("skyride.precision"), 1, writer);
	            writer.writeCloseTag(GMRFSkyrideLikelihoodParser.PRECISION_PARAMETER);
	
	            writer.writeOpenTag(GMRFSkyrideLikelihoodParser.POPULATION_TREE);
	            writer.writeIDref(TreeModel.TREE_MODEL, modelPrefix + TreeModel.TREE_MODEL);
	            writer.writeCloseTag(GMRFSkyrideLikelihoodParser.POPULATION_TREE);
	
	            writer.writeCloseTag(GMRFSkyrideLikelihoodParser.SKYLINE_LIKELIHOOD);

	            break;
	        
    		case SPECIES_YULE:
            case SPECIES_BIRTH_DEATH:
            	break;
            	
    		default: 
	            // generate a coalescent process
	            writer.writeComment("Generate a coalescent likelihood");
	            writer.writeOpenTag(
	                    CoalescentLikelihoodParser.COALESCENT_LIKELIHOOD,
	                    new Attribute[]{new Attribute.Default<String>(XMLParser.ID, modelPrefix + COALESCENT)}
	            );
	            writer.writeOpenTag(CoalescentLikelihoodParser.MODEL);
	            writeNodeHeightPriorModelRef(prior, writer);
	            writer.writeCloseTag(CoalescentLikelihoodParser.MODEL);
	            writer.writeOpenTag(CoalescentLikelihoodParser.POPULATION_TREE);
	            writer.writeIDref(TreeModel.TREE_MODEL, modelPrefix + TreeModel.TREE_MODEL);
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
            case YULE:
                writer.writeIDref(YuleModelParser.YULE_MODEL, priorPrefix + "yule");
                break;
            case BIRTH_DEATH:
                writer.writeIDref(BirthDeathGernhard08Model.BIRTH_DEATH_MODEL, priorPrefix + "birthDeath");
                break;
            default:
                throw new RuntimeException("No tree prior has been specified so cannot refer to it");
        }
    }
    
    void writeEBSPVariableDemographic(PartitionTreePrior prior, XMLWriter writer) {
    	
    	if (prior.getNodeHeightPrior() == TreePriorType.EXTENDED_SKYLINE) {
    		
	    	setModelPrefix(prior.getPrefix());
	    	    	
	    	final String tagName = VariableDemographicModelParser.MODEL_NAME;
	        writer.writeComment("Generate a variableDemographic for extended Bayesian skyline process");
	        writer.writeOpenTag(tagName, new Attribute[]{
	                new Attribute.Default<String>(XMLParser.ID, modelPrefix + VariableDemographicModelParser.demoElementName),
	                new Attribute.Default<String>(VariableDemographicModelParser.TYPE, prior.getExtendedSkylineModel().toString()),
	                // use midpoint by default (todo) would be nice to have a user 'tickable' option
	                new Attribute.Default<String>(VariableDemographicModelParser.USE_MIDPOINTS, "true")
	            }
	        );
	        
	        Parameter popSize = prior.getParameter(VariableDemographicModelParser.demoElementName + ".popSize");
	        Parameter populationMean = prior.getParameter(VariableDemographicModelParser.demoElementName + ".populationMean");
	        popSize.initial = populationMean.initial;
	        
	        writer.writeOpenTag(VariableDemographicModelParser.POPULATION_SIZES);
	        writer.writeComment("popSize value = populationMean value");
	        writer.writeTag(ParameterParser.PARAMETER,
                    new Attribute[]{
	                     new Attribute.Default<String>(XMLParser.ID, modelPrefix + VariableDemographicModelParser.demoElementName + ".popSize"), 
                         new Attribute.Default<String>(ParameterParser.VALUE, Double.toString(popSize.initial))}, true);
//	        writeParameter(popSize, -1, writer);
	        writer.writeCloseTag(VariableDemographicModelParser.POPULATION_SIZES);

//            Parameter indicators = prior.getParameter(VariableDemographicModelParser.demoElementName + ".indicators");
	        writer.writeOpenTag(VariableDemographicModelParser.INDICATOR_PARAMETER);
            writer.writeTag(ParameterParser.PARAMETER,
                    new Attribute[]{
	                     new Attribute.Default<String>(XMLParser.ID, modelPrefix + VariableDemographicModelParser.demoElementName + ".indicators"),
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
		                new Attribute.Default<String>(SpeciesBindingsParser.PLOIDY, Double.toString(prior.getTreeModel().getPloidyType().getValue()))
	            	}
	        	);
	        	writer.writeIDref(TreeModel.TREE_MODEL, prior.getTreeModel().getPrefix() + TreeModel.TREE_MODEL);
	            writer.writeCloseTag(VariableDemographicModelParser.POP_TREE);
	        }
	        
	        writer.writeCloseTag(VariableDemographicModelParser.POPULATION_TREES);
	
	        writer.writeCloseTag(tagName);
	
	        writer.writeOpenTag(CoalescentLikelihoodParser.COALESCENT_LIKELIHOOD, new Attribute.Default<String>(XMLParser.ID, modelPrefix + COALESCENT));
	        writer.writeOpenTag(CoalescentLikelihoodParser.MODEL);
	        writer.writeIDref(tagName, modelPrefix + VariableDemographicModelParser.demoElementName);
	        writer.writeCloseTag(CoalescentLikelihoodParser.MODEL);
	        writer.writeComment("Take population Tree from demographic");
	        writer.writeCloseTag(CoalescentLikelihoodParser.COALESCENT_LIKELIHOOD);
	
	        writer.writeOpenTag(SumStatisticParser.SUM_STATISTIC,
	                new Attribute[]{
	                        new Attribute.Default<String>(XMLParser.ID, modelPrefix + VariableDemographicModelParser.demoElementName + ".populationSizeChanges"),
	                        new Attribute.Default<String>("elementwise", "true")
	                });
	        writer.writeIDref(ParameterParser.PARAMETER, modelPrefix + VariableDemographicModelParser.demoElementName + ".indicators");
	        writer.writeCloseTag(SumStatisticParser.SUM_STATISTIC);
	        writer.writeOpenTag(ExponentialDistributionModel.EXPONENTIAL_DISTRIBUTION_MODEL,
	                new Attribute[]{
	                        new Attribute.Default<String>(XMLParser.ID, modelPrefix + VariableDemographicModelParser.demoElementName + ".populationMeanDist")
	                        //,new Attribute.Default<String>("elementwise", "true")
	                });
	        writer.writeOpenTag(DistributionModelParser.MEAN);
	        
	        writer.writeComment("prefer populationMean value = 1");	        
	        populationMean = prior.getParameter(VariableDemographicModelParser.demoElementName + ".populationMean");
	        writer.writeTag(ParameterParser.PARAMETER,
	                new Attribute[]{
	                        new Attribute.Default<String>(XMLParser.ID, modelPrefix + VariableDemographicModelParser.demoElementName + ".populationMean"),
	                        new Attribute.Default<String>(ParameterParser.VALUE, Double.toString(populationMean.initial))}, true);
	        
	        writer.writeCloseTag(DistributionModelParser.MEAN);
	        writer.writeCloseTag(ExponentialDistributionModel.EXPONENTIAL_DISTRIBUTION_MODEL);
    	}
    }

    void writeParameterLog(PartitionTreePrior prior, XMLWriter writer) {

        setModelPrefix(prior.getPrefix());

        switch (prior.getNodeHeightPrior()) {

            case CONSTANT:            
                writer.writeIDref(ParameterParser.PARAMETER, modelPrefix + "constant.popSize");
                break;
            case EXPONENTIAL:
                writer.writeIDref(ParameterParser.PARAMETER, modelPrefix + "exponential.popSize");
                if (prior.getParameterization() == TreePriorParameterizationType.GROWTH_RATE) {
                    writer.writeIDref(ParameterParser.PARAMETER, modelPrefix + "exponential.growthRate");
                } else {
                    writer.writeIDref(ParameterParser.PARAMETER, modelPrefix + "exponential.doublingTime");
                }
                break;
            case LOGISTIC:
                writer.writeIDref(ParameterParser.PARAMETER, modelPrefix + "logistic.popSize");
                if (prior.getParameterization() == TreePriorParameterizationType.GROWTH_RATE) {
                    writer.writeIDref(ParameterParser.PARAMETER, modelPrefix + "logistic.growthRate");
                } else {
                    writer.writeIDref(ParameterParser.PARAMETER, modelPrefix + "logistic.doublingTime");
                }
                writer.writeIDref(ParameterParser.PARAMETER, modelPrefix + "logistic.t50");
                break;
            case EXPANSION:
                writer.writeIDref(ParameterParser.PARAMETER, modelPrefix + "expansion.popSize");
                if (prior.getParameterization() == TreePriorParameterizationType.GROWTH_RATE) {
                    writer.writeIDref(ParameterParser.PARAMETER, modelPrefix + "expansion.growthRate");
                } else {
                    writer.writeIDref(ParameterParser.PARAMETER, modelPrefix + "expansion.doublingTime");
                }
                writer.writeIDref(ParameterParser.PARAMETER, modelPrefix + "expansion.ancestralProportion");
                break;
            case SKYLINE:
                writer.writeIDref(ParameterParser.PARAMETER, modelPrefix + "skyline.popSize");
                writer.writeIDref(ParameterParser.PARAMETER, modelPrefix + "skyline.groupSize");
                break;
            case EXTENDED_SKYLINE:
                writer.writeIDref(SumStatisticParser.SUM_STATISTIC, "demographic.populationSizeChanges");
                writer.writeIDref(ParameterParser.PARAMETER, modelPrefix + "demographic.populationMean");
                writer.writeIDref(ParameterParser.PARAMETER, modelPrefix + "demographic.popSize");
                writer.writeIDref(ParameterParser.PARAMETER, modelPrefix + "demographic.indicators");
                break;
            case GMRF_SKYRIDE:
                writer.writeIDref(ParameterParser.PARAMETER, modelPrefix + "skyride.precision");
                writer.writeIDref(ParameterParser.PARAMETER, modelPrefix + "skyride.logPopSize");
                writer.writeIDref(ParameterParser.PARAMETER, modelPrefix + "skyride.groupSize");
                break;
            case YULE:
                writer.writeIDref(ParameterParser.PARAMETER, modelPrefix + "yule.birthRate");
                break;
            case BIRTH_DEATH:
                writer.writeIDref(ParameterParser.PARAMETER, modelPrefix + BirthDeathModelParser.MEAN_GROWTH_RATE_PARAM_NAME);
                writer.writeIDref(ParameterParser.PARAMETER, modelPrefix + BirthDeathModelParser.RELATIVE_DEATH_RATE_PARAM_NAME);
                break;
            case SPECIES_YULE:
            case SPECIES_BIRTH_DEATH:
            	break;
            default:
                throw new RuntimeException("No tree prior has been specified so cannot refer to it");
        }

    }

    void writeEBSPAnalysisToCSVfile(PartitionTreePrior prior, XMLWriter writer) {
    	
    	setModelPrefix(prior.getPrefix());
    	
        String logFileName = options.logFileName;

        if (prior.getNodeHeightPrior() == TreePriorType.EXTENDED_SKYLINE) {
            writer.writeOpenTag(EBSPAnalysisParser.VD_ANALYSIS, new Attribute[]{
                    new Attribute.Default<String>(XMLParser.ID, modelPrefix + "demographic.analysis"),
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
            writer.writeIDref(EBSPAnalysisParser.VD_ANALYSIS, modelPrefix + "demographic.analysis");
            writer.writeCloseTag(CSVExporterParser.COLUMNS);
            writer.writeCloseTag(CSVExporterParser.CSV_EXPORT);
        }
    }

    private void writeExponentialMarkovLikelihood(XMLWriter writer) {
        writer.writeOpenTag(
                ExponentialMarkovModel.EXPONENTIAL_MARKOV_MODEL,
                new Attribute[]{new Attribute.Default<String>(XMLParser.ID, modelPrefix + "eml1"),
                        new Attribute.Default<String>("jeffreys", "true")}
        );
        writer.writeOpenTag(ExponentialMarkovModelParser.CHAIN_PARAMETER);
        writer.writeIDref(ParameterParser.PARAMETER, modelPrefix + "skyline.popSize");
        writer.writeCloseTag(ExponentialMarkovModelParser.CHAIN_PARAMETER);
        writer.writeCloseTag(ExponentialMarkovModel.EXPONENTIAL_MARKOV_MODEL);
    }

    public void writePriorLikelihoodReferenceLog(PartitionTreePrior prior, PartitionTreeModel model, XMLWriter writer) {
    	//tree model prefix
    	setModelPrefix(model.getPrefix()); // only has prefix, if (options.getPartitionTreePriors().size() > 1)
    	
        switch (prior.getNodeHeightPrior()) {

            case YULE:
            case BIRTH_DEATH:
                writer.writeIDref(SpeciationLikelihoodParser.SPECIATION_LIKELIHOOD, modelPrefix + "speciation");
                break;
            case SKYLINE:
                writer.writeIDref(BayesianSkylineLikelihoodParser.SKYLINE_LIKELIHOOD, modelPrefix + "skyline");
//                writer.writeIDref(ExponentialMarkovModel.EXPONENTIAL_MARKOV_MODEL, modelPrefix + "eml1");
                break;
            case GMRF_SKYRIDE:
                writer.writeIDref(GMRFSkyrideLikelihoodParser.SKYLINE_LIKELIHOOD, modelPrefix + "skyride");
                break;
            case LOGISTIC:
//                writer.writeIDref(BooleanLikelihoodParser.BOOLEAN_LIKELIHOOD, modelPrefix + "booleanLikelihood1");
                writer.writeIDref(CoalescentLikelihoodParser.COALESCENT_LIKELIHOOD, modelPrefix + COALESCENT);
                break;
            case EXTENDED_SKYLINE:
            	// only 1 coalescent, so write it separately after this method
            case SPECIES_YULE:
            case SPECIES_BIRTH_DEATH:
                // do not need
                break;
            default: 
                writer.writeIDref(CoalescentLikelihoodParser.COALESCENT_LIKELIHOOD, modelPrefix + COALESCENT);
        }
    }
    
    // id is written in writePriorLikelihood (PartitionTreePrior prior, PartitionTreeModel model, XMLWriter writer)
    public void writePriorLikelihoodReference(PartitionTreePrior prior, PartitionTreeModel model, XMLWriter writer) {
        //tree model prefix
    	setModelPrefix(model.getPrefix()); // only has prefix, if (options.getPartitionTreePriors().size() > 1)
    	
        switch (prior.getNodeHeightPrior()) {

            case YULE:
            case BIRTH_DEATH:
                writer.writeIDref(SpeciationLikelihoodParser.SPECIATION_LIKELIHOOD, modelPrefix + "speciation");
                break;
            case SKYLINE:
                writer.writeIDref(BayesianSkylineLikelihoodParser.SKYLINE_LIKELIHOOD, modelPrefix + "skyline");
                writer.writeIDref(ExponentialMarkovModel.EXPONENTIAL_MARKOV_MODEL, modelPrefix + "eml1");
                break;
            case GMRF_SKYRIDE:
                writer.writeIDref(GMRFSkyrideLikelihoodParser.SKYLINE_LIKELIHOOD, modelPrefix + "skyride");
                break;
            case LOGISTIC:
                writer.writeIDref(BooleanLikelihoodParser.BOOLEAN_LIKELIHOOD, modelPrefix + "booleanLikelihood1");
                writer.writeIDref(CoalescentLikelihoodParser.COALESCENT_LIKELIHOOD, modelPrefix + COALESCENT);
                break;
            case EXTENDED_SKYLINE:
            	// only 1 coalescent, so write it in writeEBSPVariableDemographicReference
            case SPECIES_YULE:
            case SPECIES_BIRTH_DEATH:
                // do not need
                break;
            default: 
                writer.writeIDref(CoalescentLikelihoodParser.COALESCENT_LIKELIHOOD, modelPrefix + COALESCENT);
        }        
    }
    
    public void writeEBSPVariableDemographicReference(PartitionTreePrior prior, XMLWriter writer) {
        
        setModelPrefix(prior.getPrefix());
        	
        //TODO: make suitable for *BEAST
        if (prior.getNodeHeightPrior() == TreePriorType.EXTENDED_SKYLINE) { 
        	
        	writer.writeIDref(CoalescentLikelihoodParser.COALESCENT_LIKELIHOOD, modelPrefix + COALESCENT); // only 1 coalescent
        	
            writer.writeOpenTag(MixedDistributionLikelihoodParser.DISTRIBUTION_LIKELIHOOD);

            writer.writeOpenTag(MixedDistributionLikelihoodParser.DISTRIBUTION0);
            writer.writeIDref(ExponentialDistributionModel.EXPONENTIAL_DISTRIBUTION_MODEL, modelPrefix + "demographic.populationMeanDist");
            writer.writeCloseTag(MixedDistributionLikelihoodParser.DISTRIBUTION0);

            writer.writeOpenTag(MixedDistributionLikelihoodParser.DISTRIBUTION1);
            writer.writeIDref(ExponentialDistributionModel.EXPONENTIAL_DISTRIBUTION_MODEL, modelPrefix + "demographic.populationMeanDist");
            writer.writeCloseTag(MixedDistributionLikelihoodParser.DISTRIBUTION1);

            writer.writeOpenTag(MixedDistributionLikelihoodParser.DATA);

            writer.writeIDref(ParameterParser.PARAMETER, modelPrefix + "demographic.popSize");

            writer.writeCloseTag(MixedDistributionLikelihoodParser.DATA);

            writer.writeOpenTag(MixedDistributionLikelihoodParser.INDICATORS);

            writer.writeIDref(ParameterParser.PARAMETER, modelPrefix + "demographic.indicators");

            writer.writeCloseTag(MixedDistributionLikelihoodParser.INDICATORS);

            writer.writeCloseTag(MixedDistributionLikelihoodParser.DISTRIBUTION_LIKELIHOOD);
        }

    }
}
