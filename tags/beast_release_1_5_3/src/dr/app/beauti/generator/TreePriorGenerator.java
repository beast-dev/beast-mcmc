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
import dr.evomodel.coalescent.*;
import dr.evomodel.speciation.BirthDeathGernhard08Model;
import dr.evomodel.speciation.SpeciationLikelihood;
import dr.evomodel.speciation.SpeciesBindings;
import dr.evomodel.tree.TreeModel;
import dr.evomodelxml.BirthDeathModelParser;
import dr.evomodelxml.CSVExporterParser;
import dr.evomodelxml.YuleModelParser;
import dr.inference.distribution.ExponentialDistributionModel;
import dr.inference.distribution.ExponentialMarkovModel;
import dr.inference.distribution.MixedDistributionLikelihood;
import dr.inference.model.BooleanLikelihood;
import dr.inference.model.ParameterParser;
import dr.inference.model.SumStatistic;
import dr.inference.model.TestStatistic;
import dr.inferencexml.DistributionModelParser;
import dr.inferencexml.ExponentialMarkovModelParser;
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
	                    ConstantPopulationModel.CONSTANT_POPULATION_MODEL,
	                    new Attribute[]{
	                            new Attribute.Default<String>(XMLParser.ID, modelPrefix + "constant"),
	                            new Attribute.Default<String>("units", Units.Utils.getDefaultUnitName(options.units))
	                    }
	            );
	
	            writer.writeOpenTag(ConstantPopulationModel.POPULATION_SIZE);
	            writeParameter("constant.popSize", prior, writer);
	            writer.writeCloseTag(ConstantPopulationModel.POPULATION_SIZE);
	            writer.writeCloseTag(ConstantPopulationModel.CONSTANT_POPULATION_MODEL);
	            
	            break;

        	case EXPONENTIAL:
	            // generate an exponential prior tree
	
	            writer.writeComment("A prior assumption that the population size has grown exponentially");
	            writer.writeComment("throughout the time spanned by the genealogy.");
	            writer.writeOpenTag(
	                    ExponentialGrowthModel.EXPONENTIAL_GROWTH_MODEL,
	                    new Attribute[]{
	                            new Attribute.Default<String>(XMLParser.ID, modelPrefix + "exponential"),
	                            new Attribute.Default<String>("units", Units.Utils.getDefaultUnitName(options.units))
	                    }
	            );
	
	            // write pop size socket
	            writer.writeOpenTag(ExponentialGrowthModel.POPULATION_SIZE);
	            writeParameter("exponential.popSize", prior, writer);
	            writer.writeCloseTag(ExponentialGrowthModel.POPULATION_SIZE);
	
	            if (parameterization == TreePriorParameterizationType.GROWTH_RATE) {
	                // write growth rate socket
	                writer.writeOpenTag(ExponentialGrowthModel.GROWTH_RATE);
	                writeParameter("exponential.growthRate", prior, writer);
	                writer.writeCloseTag(ExponentialGrowthModel.GROWTH_RATE);
	            } else {
	                // write doubling time socket
	                writer.writeOpenTag(ExponentialGrowthModel.DOUBLING_TIME);
	                writeParameter("exponential.doublingTime", prior, writer);
	                writer.writeCloseTag(ExponentialGrowthModel.DOUBLING_TIME);
	            }
	
	            writer.writeCloseTag(ExponentialGrowthModel.EXPONENTIAL_GROWTH_MODEL);
	            
	            break;
	            
        	case LOGISTIC:
	            // generate an exponential prior tree
	
	            writer.writeComment("A prior assumption that the population size has grown logistically");
	            writer.writeComment("throughout the time spanned by the genealogy.");
	            writer.writeOpenTag(
	                    LogisticGrowthModel.LOGISTIC_GROWTH_MODEL,
	                    new Attribute[]{
	                            new Attribute.Default<String>(XMLParser.ID, modelPrefix + "logistic"),
	                            new Attribute.Default<String>("units", Units.Utils.getDefaultUnitName(options.units))
	                    }
	            );
	
	            // write pop size socket
	            writer.writeOpenTag(LogisticGrowthModel.POPULATION_SIZE);
	            writeParameter("logistic.popSize", prior, writer);
	            writer.writeCloseTag(LogisticGrowthModel.POPULATION_SIZE);
	
	            if (parameterization == TreePriorParameterizationType.GROWTH_RATE) {
	                // write growth rate socket
	                writer.writeOpenTag(LogisticGrowthModel.GROWTH_RATE);
	                writeParameter("logistic.growthRate", prior, writer);
	                writer.writeCloseTag(LogisticGrowthModel.GROWTH_RATE);
	            } else {
	                // write doubling time socket
	                writer.writeOpenTag(LogisticGrowthModel.DOUBLING_TIME);
	                writeParameter("logistic.doublingTime", prior, writer);
	                writer.writeCloseTag(LogisticGrowthModel.DOUBLING_TIME);
	            }
	
	            // write logistic t50 socket
	            writer.writeOpenTag(LogisticGrowthModel.TIME_50);
	            
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
	            writer.writeCloseTag(LogisticGrowthModel.TIME_50);
	
	            writer.writeCloseTag(LogisticGrowthModel.LOGISTIC_GROWTH_MODEL);
	
	            initialPopSize = "logistic.popSize";
	            
	            break;

        	case EXPANSION:
        		// generate an exponential prior tree
	            writer.writeComment("A prior assumption that the population size has grown exponentially");
	            writer.writeComment("from some ancestral population size in the past.");
	            writer.writeOpenTag(
	                    ExpansionModel.EXPANSION_MODEL,
	                    new Attribute[]{
	                            new Attribute.Default<String>(XMLParser.ID, modelPrefix + "expansion"),
	                            new Attribute.Default<String>("units", Units.Utils.getDefaultUnitName(options.units))
	                    }
	            );
	
	            // write pop size socket
	            writeParameter(ExpansionModel.POPULATION_SIZE, "expansion.popSize", prior, writer);
	
	            if (parameterization == TreePriorParameterizationType.GROWTH_RATE) {
	                // write growth rate socket
	                writeParameter(ExpansionModel.GROWTH_RATE, "expansion.growthRate", prior, writer);
	            } else {
	                // write doubling time socket
	                writeParameter(ExpansionModel.DOUBLING_TIME, "expansion.doublingTime", prior, writer);
	            }
	
	            // write ancestral proportion socket
	            writeParameter(ExpansionModel.ANCESTRAL_POPULATION_PROPORTION, "expansion.ancestralProportion", prior, writer);
	
	            writer.writeCloseTag(ExpansionModel.EXPANSION_MODEL);
	
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
	                    ConstantPopulationModel.CONSTANT_POPULATION_MODEL,
	                    new Attribute[]{
	                            new Attribute.Default<String>(XMLParser.ID, modelPrefix + "constant"),
	                            new Attribute.Default<String>("units", Units.Utils.getDefaultUnitName(options.units))
	                    }
	            );
	            
	            // initial value for pop mean is the same as what used to be the value for the population size
	            Parameter para = options.starBEASTOptions.getParameter(TraitGuesser.Traits.TRAIT_SPECIES + "." + options.starBEASTOptions.POP_MEAN);
	            prior.getParameter("constant.popSize").initial = para.initial;
	
	            writer.writeOpenTag(ConstantPopulationModel.POPULATION_SIZE);
	            writeParameter("constant.popSize", prior, writer);
	            writer.writeCloseTag(ConstantPopulationModel.POPULATION_SIZE);
	            writer.writeCloseTag(ConstantPopulationModel.CONSTANT_POPULATION_MODEL);
	            
	            break;	            
        }

        if ((!options.starBEASTOptions.isSpeciesAnalysis()) && nodeHeightPrior != TreePriorType.CONSTANT && nodeHeightPrior != TreePriorType.EXPONENTIAL) {
            // If the node height prior is not one of these two then we need to simulate a
            // random starting tree under a constant size coalescent.

            writer.writeComment("This is a simple constant population size coalescent model");
            writer.writeComment("that is used to generate an initial tree for the chain.");
            writer.writeOpenTag(
                    ConstantPopulationModel.CONSTANT_POPULATION_MODEL,
                    new Attribute[]{
                            new Attribute.Default<String>(XMLParser.ID, modelPrefix + "initialDemo"),
                            new Attribute.Default<String>("units", Units.Utils.getDefaultUnitName(units))
                    }
            );

            writer.writeOpenTag(ConstantPopulationModel.POPULATION_SIZE);
            if (initialPopSize != null) {
            	writer.writeIDref(ParameterParser.PARAMETER, modelPrefix + initialPopSize);
            } else {
                writeParameter(modelPrefix + "initialDemo.popSize", 1, 100.0, Double.NaN, Double.NaN, writer);
            }
            writer.writeCloseTag(ConstantPopulationModel.POPULATION_SIZE);
            writer.writeCloseTag(ConstantPopulationModel.CONSTANT_POPULATION_MODEL);
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
	                    SpeciationLikelihood.SPECIATION_LIKELIHOOD,
	                    new Attribute[]{
	                            new Attribute.Default<String>(XMLParser.ID, modelPrefix + "speciation")
	                    }
	            );
	
	            // write pop size socket
	            writer.writeOpenTag(SpeciationLikelihood.MODEL);
	            writeNodeHeightPriorModelRef(prior, writer);
	            writer.writeCloseTag(SpeciationLikelihood.MODEL);
	            writer.writeOpenTag(SpeciationLikelihood.TREE);
	            writer.writeIDref(TreeModel.TREE_MODEL, modelPrefix + TreeModel.TREE_MODEL);
	            writer.writeCloseTag(SpeciationLikelihood.TREE);
	
	            writer.writeCloseTag(SpeciationLikelihood.SPECIATION_LIKELIHOOD);
	            
	            break;
	            
    		case LOGISTIC:
    			writer.writeComment("Generate a boolean likelihood for Coalescent: Logistic Growth");
    			writer.writeOpenTag(
    	                BooleanLikelihood.BOOLEAN_LIKELIHOOD,
    	                new Attribute[]{new Attribute.Default<String>(XMLParser.ID, modelPrefix + "booleanLikelihood1")}
    	        );
    	        writer.writeOpenTag(
    	                TestStatistic.TEST_STATISTIC,
    	                new Attribute[]{
    	                        new Attribute.Default<String>(XMLParser.ID, "test1"),
    	                        new Attribute.Default<String>("name", "test1")
    	                }
    	        );
    	        writer.writeIDref(ParameterParser.PARAMETER, priorPrefix + "logistic.t50"); //TODO correct?
    	        writer.writeOpenTag("lessThan");
    	        writer.writeIDref(ParameterParser.PARAMETER, modelPrefix + "treeModel.rootHeight");
    	        writer.writeCloseTag("lessThan");
    	        writer.writeCloseTag(TestStatistic.TEST_STATISTIC);
    	        writer.writeCloseTag(BooleanLikelihood.BOOLEAN_LIKELIHOOD);
    	        
    	        writer.writeOpenTag(
	                    CoalescentLikelihood.COALESCENT_LIKELIHOOD,
	                    new Attribute[]{new Attribute.Default<String>(XMLParser.ID, modelPrefix + COALESCENT)}
	            );
	            writer.writeOpenTag(CoalescentLikelihood.MODEL);
	            writeNodeHeightPriorModelRef(prior, writer);
	            writer.writeCloseTag(CoalescentLikelihood.MODEL);
	            writer.writeOpenTag(CoalescentLikelihood.POPULATION_TREE);
	            writer.writeIDref(TreeModel.TREE_MODEL, modelPrefix + TreeModel.TREE_MODEL);
	            writer.writeCloseTag(CoalescentLikelihood.POPULATION_TREE);
	            writer.writeCloseTag(CoalescentLikelihood.COALESCENT_LIKELIHOOD);
	            
    	        break;
	            
    		case SKYLINE:
	            // generate a Bayesian skyline plot
	            writer.writeComment("Generate a generalizedSkyLineLikelihood for Bayesian Skyline");
	            writer.writeOpenTag(
	                    BayesianSkylineLikelihood.SKYLINE_LIKELIHOOD,
	                    new Attribute[]{
	                            new Attribute.Default<String>(XMLParser.ID, modelPrefix + "skyline"),
	                            new Attribute.Default<String>("linear",
	                                    prior.getSkylineModel() == TreePriorParameterizationType.LINEAR_SKYLINE ? "true" : "false")
	                    }
	            );
	
	            // write pop size socket
	            writer.writeOpenTag(BayesianSkylineLikelihood.POPULATION_SIZES);
	            if (prior.getSkylineModel() == TreePriorParameterizationType.LINEAR_SKYLINE) {
	                writeParameter(prior.getParameter("skyline.popSize"), prior.getSkylineGroupCount() + 1, writer);
	            } else {
	                writeParameter(prior.getParameter("skyline.popSize"), prior.getSkylineGroupCount(), writer);
	            }
	            writer.writeCloseTag(BayesianSkylineLikelihood.POPULATION_SIZES);
	
	            // write group size socket
	            writer.writeOpenTag(BayesianSkylineLikelihood.GROUP_SIZES);
	            writeParameter(prior.getParameter("skyline.groupSize"), prior.getSkylineGroupCount(), writer);
	            writer.writeCloseTag(BayesianSkylineLikelihood.GROUP_SIZES);
	
	            writer.writeOpenTag(CoalescentLikelihood.POPULATION_TREE);
	            writer.writeIDref(TreeModel.TREE_MODEL, modelPrefix + TreeModel.TREE_MODEL);
	            writer.writeCloseTag(CoalescentLikelihood.POPULATION_TREE);
	
	            writer.writeCloseTag(BayesianSkylineLikelihood.SKYLINE_LIKELIHOOD);
	            
	            writer.writeText("");
	            writeExponentialMarkovLikelihood(writer);
	            
	            break;
	            
    		case EXTENDED_SKYLINE:        	
	            // different format	            
	            break;
            
    		case GMRF_SKYRIDE:
	            writer.writeComment("Generate a gmrfSkyrideLikelihood for GMRF Bayesian Skyride process");
	            writer.writeOpenTag(
	                    GMRFSkyrideLikelihood.SKYLINE_LIKELIHOOD,
	                    new Attribute[]{
	                            new Attribute.Default<String>(XMLParser.ID, modelPrefix + "skyride"),
	                            new Attribute.Default<String>(GMRFSkyrideLikelihood.TIME_AWARE_SMOOTHING,
	                                    prior.getSkyrideSmoothing() == TreePriorParameterizationType.TIME_AWARE_SKYRIDE ? "true" : "false"),
	                            new Attribute.Default<String>(GMRFSkyrideLikelihood.RANDOMIZE_TREE,
	                                    //TODO For GMRF, tree model/tree prior combination not implemented by BEAST yet. The validation is in BeastGenerator.checkOptions()
	                                    prior.getTreeModel().getStartingTreeType() == StartingTreeType.UPGMA ? "true" : "false"),
	                    }
	            );
	
	            int skyrideIntervalCount = options.taxonList.getTaxonCount() - 1;
	            writer.writeOpenTag(GMRFSkyrideLikelihood.POPULATION_PARAMETER);
                writer.writeComment("skyride.logPopSize is log unit unlike other popSize");
	            writeParameter(prior.getParameter("skyride.logPopSize"), skyrideIntervalCount, writer);
	            writer.writeCloseTag(GMRFSkyrideLikelihood.POPULATION_PARAMETER);
	
	            writer.writeOpenTag(GMRFSkyrideLikelihood.GROUP_SIZES);
	            writeParameter(prior.getParameter("skyride.groupSize"), skyrideIntervalCount, writer);
	            writer.writeCloseTag(GMRFSkyrideLikelihood.GROUP_SIZES);
	
	            writer.writeOpenTag(GMRFSkyrideLikelihood.PRECISION_PARAMETER);
	            writeParameter(prior.getParameter("skyride.precision"), 1, writer);
	            writer.writeCloseTag(GMRFSkyrideLikelihood.PRECISION_PARAMETER);
	
	            writer.writeOpenTag(GMRFSkyrideLikelihood.POPULATION_TREE);
	            writer.writeIDref(TreeModel.TREE_MODEL, modelPrefix + TreeModel.TREE_MODEL);
	            writer.writeCloseTag(GMRFSkyrideLikelihood.POPULATION_TREE);
	
	            writer.writeCloseTag(GMRFSkyrideLikelihood.SKYLINE_LIKELIHOOD);

	            break;
	        
    		case SPECIES_YULE:
            case SPECIES_BIRTH_DEATH:
            	break;
            	
    		default: 
	            // generate a coalescent process
	            writer.writeComment("Generate a coalescent likelihood");
	            writer.writeOpenTag(
	                    CoalescentLikelihood.COALESCENT_LIKELIHOOD,
	                    new Attribute[]{new Attribute.Default<String>(XMLParser.ID, modelPrefix + COALESCENT)}
	            );
	            writer.writeOpenTag(CoalescentLikelihood.MODEL);
	            writeNodeHeightPriorModelRef(prior, writer);
	            writer.writeCloseTag(CoalescentLikelihood.MODEL);
	            writer.writeOpenTag(CoalescentLikelihood.POPULATION_TREE);
	            writer.writeIDref(TreeModel.TREE_MODEL, modelPrefix + TreeModel.TREE_MODEL);
	            writer.writeCloseTag(CoalescentLikelihood.POPULATION_TREE);
	            writer.writeCloseTag(CoalescentLikelihood.COALESCENT_LIKELIHOOD);
        }
    }

    void writeNodeHeightPriorModelRef(PartitionTreePrior prior, XMLWriter writer) {    	
    	TreePriorType treePrior = prior.getNodeHeightPrior();
    	String priorPrefix = prior.getPrefix();
    	
        switch (treePrior) {
            case CONSTANT:
            case SPECIES_YULE:
            case SPECIES_BIRTH_DEATH:
                writer.writeIDref(ConstantPopulationModel.CONSTANT_POPULATION_MODEL, priorPrefix + "constant");
                break;
            case EXPONENTIAL:
                writer.writeIDref(ExponentialGrowthModel.EXPONENTIAL_GROWTH_MODEL, priorPrefix + "exponential");
                break;
            case LOGISTIC:
                writer.writeIDref(LogisticGrowthModel.LOGISTIC_GROWTH_MODEL, priorPrefix + "logistic");
                break;
            case EXPANSION:
                writer.writeIDref(ExpansionModel.EXPANSION_MODEL, priorPrefix + "expansion");
                break;
            case SKYLINE:
                writer.writeIDref(BayesianSkylineLikelihood.SKYLINE_LIKELIHOOD, priorPrefix + "skyline");
                break;
            case GMRF_SKYRIDE:
                writer.writeIDref(GMRFSkyrideLikelihood.SKYLINE_LIKELIHOOD, priorPrefix + "skyride");
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
	    	    	
	    	final String tagName = VariableDemographicModel.PARSER.getParserName();
	        writer.writeComment("Generate a variableDemographic for extended Bayesian skyline process");
	        writer.writeOpenTag(tagName, new Attribute[]{
	                new Attribute.Default<String>(XMLParser.ID, modelPrefix + VariableDemographicModel.demoElementName),
	                new Attribute.Default<String>(VariableDemographicModel.TYPE, prior.getExtendedSkylineModel().toString()),
	                // use midpoint by default (todo) would be nice to have a user 'tickable' option
	                new Attribute.Default<String>(VariableDemographicModel.USE_MIDPOINTS, "true")
	            }
	        );
	        
//	        Parameter popSize = prior.getParameter(VariableDemographicModel.demoElementName + ".popSize");
//	        Parameter populationMean = prior.getParameter(VariableDemographicModel.demoElementName + ".populationMean");
//	        popSize.initial = populationMean.initial;
	        
	        writer.writeOpenTag(VariableDemographicModel.POPULATION_SIZES);
//	        writer.writeComment("popSize value = populationMean value");
	        writer.writeTag(ParameterParser.PARAMETER,
	                        new Attribute.Default<String>(XMLParser.ID, modelPrefix + VariableDemographicModel.demoElementName + ".popSize"), true);
//	        writeParameter(popSize, -1, writer); 
	        writer.writeCloseTag(VariableDemographicModel.POPULATION_SIZES);
	
	        writer.writeOpenTag(VariableDemographicModel.INDICATOR_PARAMETER);
	        writeParameter(prior.getParameter(VariableDemographicModel.demoElementName + ".indicators"), -1, writer); // not need dimension
	        writer.writeCloseTag(VariableDemographicModel.INDICATOR_PARAMETER);
	
	        writer.writeOpenTag(VariableDemographicModel.POPULATION_TREES);
	        
	        if (options.isShareSameTreePrior()) {
	            for (PartitionTreeModel model : options.getPartitionTreeModels()) {
		            writer.writeOpenTag(VariableDemographicModel.POP_TREE, new Attribute[]{
			                new Attribute.Default<String>(SpeciesBindings.PLOIDY, Double.toString(model.getPloidyType().getValue()))
			            }
			        );
		            writer.writeIDref(TreeModel.TREE_MODEL, model.getPrefix() + TreeModel.TREE_MODEL);
		            writer.writeCloseTag(VariableDemographicModel.POP_TREE);
	            }
	        } else {//TODO correct for not sharing same prior?
	        	writer.writeOpenTag(VariableDemographicModel.POP_TREE, new Attribute[]{
		                new Attribute.Default<String>(SpeciesBindings.PLOIDY, Double.toString(prior.getTreeModel().getPloidyType().getValue()))
	            	}
	        	);
	        	writer.writeIDref(TreeModel.TREE_MODEL, prior.getTreeModel().getPrefix() + TreeModel.TREE_MODEL);
	            writer.writeCloseTag(VariableDemographicModel.POP_TREE);
	        }
	        
	        writer.writeCloseTag(VariableDemographicModel.POPULATION_TREES);
	
	        writer.writeCloseTag(tagName);
	
	        writer.writeOpenTag(CoalescentLikelihood.COALESCENT_LIKELIHOOD, new Attribute.Default<String>(XMLParser.ID, modelPrefix + COALESCENT));
	        writer.writeOpenTag(CoalescentLikelihood.MODEL);
	        writer.writeIDref(tagName, modelPrefix + VariableDemographicModel.demoElementName);
	        writer.writeCloseTag(CoalescentLikelihood.MODEL);
	        writer.writeComment("Take population Tree from demographic");
	        writer.writeCloseTag(CoalescentLikelihood.COALESCENT_LIKELIHOOD);
	
	        writer.writeOpenTag(SumStatistic.SUM_STATISTIC,
	                new Attribute[]{
	                        new Attribute.Default<String>(XMLParser.ID, modelPrefix + VariableDemographicModel.demoElementName + ".populationSizeChanges"),
	                        new Attribute.Default<String>("elementwise", "true")
	                });
	        writer.writeIDref(ParameterParser.PARAMETER, modelPrefix + VariableDemographicModel.demoElementName + ".indicators");
	        writer.writeCloseTag(SumStatistic.SUM_STATISTIC);
	        writer.writeOpenTag(ExponentialDistributionModel.EXPONENTIAL_DISTRIBUTION_MODEL,
	                new Attribute[]{
	                        new Attribute.Default<String>(XMLParser.ID, modelPrefix + VariableDemographicModel.demoElementName + ".populationMeanDist")
	                        //,new Attribute.Default<String>("elementwise", "true")
	                });
	        writer.writeOpenTag(DistributionModelParser.MEAN);
	        
	        writer.writeComment("prefer populationMean value = 1");	        
	        Parameter populationMean = prior.getParameter(VariableDemographicModel.demoElementName + ".populationMean");
	        writer.writeTag(ParameterParser.PARAMETER,
	                new Attribute[]{
	                        new Attribute.Default<String>(XMLParser.ID, modelPrefix + VariableDemographicModel.demoElementName + ".populationMean"),
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
                writeSumStatisticColumn(writer, "demographic.populationSizeChanges", "popSize_changes");
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
            writer.writeOpenTag(EBSPAnalysis.VD_ANALYSIS, new Attribute[]{
                    new Attribute.Default<String>(XMLParser.ID, modelPrefix + "demographic.analysis"),
                    new Attribute.Default<Double>(EBSPAnalysis.BURN_IN, 0.1),
                    new Attribute.Default<Boolean>(VariableDemographicModel.USE_MIDPOINTS, true)}
            );

            writer.writeOpenTag(EBSPAnalysis.LOG_FILE_NAME);
            writer.writeText(logFileName);
            writer.writeCloseTag(EBSPAnalysis.LOG_FILE_NAME);

            writer.writeOpenTag(EBSPAnalysis.TREE_FILE_NAMES);
	            for (String treeFN : options.treeFileName) {
		            writer.writeOpenTag(EBSPAnalysis.TREE_LOG);
		            writer.writeText(treeFN);
		            writer.writeCloseTag(EBSPAnalysis.TREE_LOG);
	            }
            writer.writeCloseTag(EBSPAnalysis.TREE_FILE_NAMES);

            writer.writeOpenTag(EBSPAnalysis.MODEL_TYPE);
            writer.writeText(prior.getExtendedSkylineModel().toString());
            writer.writeCloseTag(EBSPAnalysis.MODEL_TYPE);

            writer.writeOpenTag(EBSPAnalysis.POPULATION_FIRST_COLUMN);
            writer.writeText(VariableDemographicModel.demoElementName + ".popSize1");
            writer.writeCloseTag(EBSPAnalysis.POPULATION_FIRST_COLUMN);

            writer.writeOpenTag(EBSPAnalysis.INDICATORS_FIRST_COLUMN);
            writer.writeText(VariableDemographicModel.demoElementName + ".indicators1");
            writer.writeCloseTag(EBSPAnalysis.INDICATORS_FIRST_COLUMN);

            writer.writeCloseTag(EBSPAnalysis.VD_ANALYSIS);

            writer.writeOpenTag(CSVExporterParser.CSV_EXPORT,
                    new Attribute[]{
                            new Attribute.Default<String>(CSVExporterParser.FILE_NAME,
                                    logFileName.subSequence(0, logFileName.length() - 4) + ".csv"), //.log
                            new Attribute.Default<String>(CSVExporterParser.SEPARATOR, ",")
                    });
            writer.writeOpenTag(CSVExporterParser.COLUMNS);
            writer.writeIDref(EBSPAnalysis.VD_ANALYSIS, modelPrefix + "demographic.analysis");
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
                writer.writeIDref(SpeciationLikelihood.SPECIATION_LIKELIHOOD, modelPrefix + "speciation");
                break;
            case SKYLINE:
                writer.writeIDref(BayesianSkylineLikelihood.SKYLINE_LIKELIHOOD, modelPrefix + "skyline");
//                writer.writeIDref(ExponentialMarkovModel.EXPONENTIAL_MARKOV_MODEL, modelPrefix + "eml1");
                break;
            case GMRF_SKYRIDE:
                writer.writeIDref(GMRFSkyrideLikelihood.SKYLINE_LIKELIHOOD, modelPrefix + "skyride");
                break;
            case LOGISTIC:
//                writer.writeIDref(BooleanLikelihood.BOOLEAN_LIKELIHOOD, modelPrefix + "booleanLikelihood1");
                writer.writeIDref(CoalescentLikelihood.COALESCENT_LIKELIHOOD, modelPrefix + COALESCENT);
                break;
            case EXTENDED_SKYLINE:
            	// only 1 coalescent, so write it separately after this method
            case SPECIES_YULE:
            case SPECIES_BIRTH_DEATH:
                // do not need
                break;
            default: 
                writer.writeIDref(CoalescentLikelihood.COALESCENT_LIKELIHOOD, modelPrefix + COALESCENT);
        }
    }
    
    // id is written in writePriorLikelihood (PartitionTreePrior prior, PartitionTreeModel model, XMLWriter writer)
    public void writePriorLikelihoodReference(PartitionTreePrior prior, PartitionTreeModel model, XMLWriter writer) {
        //tree model prefix
    	setModelPrefix(model.getPrefix()); // only has prefix, if (options.getPartitionTreePriors().size() > 1)
    	
        switch (prior.getNodeHeightPrior()) {

            case YULE:
            case BIRTH_DEATH:
                writer.writeIDref(SpeciationLikelihood.SPECIATION_LIKELIHOOD, modelPrefix + "speciation");
                break;
            case SKYLINE:
                writer.writeIDref(BayesianSkylineLikelihood.SKYLINE_LIKELIHOOD, modelPrefix + "skyline");
                writer.writeIDref(ExponentialMarkovModel.EXPONENTIAL_MARKOV_MODEL, modelPrefix + "eml1");
                break;
            case GMRF_SKYRIDE:
                writer.writeIDref(GMRFSkyrideLikelihood.SKYLINE_LIKELIHOOD, modelPrefix + "skyride");
                break;
            case LOGISTIC:
                writer.writeIDref(BooleanLikelihood.BOOLEAN_LIKELIHOOD, modelPrefix + "booleanLikelihood1");
                writer.writeIDref(CoalescentLikelihood.COALESCENT_LIKELIHOOD, modelPrefix + COALESCENT);
                break;
            case EXTENDED_SKYLINE:
            	// only 1 coalescent, so write it in writeEBSPVariableDemographicReference
            case SPECIES_YULE:
            case SPECIES_BIRTH_DEATH:
                // do not need
                break;
            default: 
                writer.writeIDref(CoalescentLikelihood.COALESCENT_LIKELIHOOD, modelPrefix + COALESCENT);
        }        
    }
    
    public void writeEBSPVariableDemographicReference(PartitionTreePrior prior, XMLWriter writer) {
        
        setModelPrefix(prior.getPrefix());
        	
        //TODO: make suitable for *BEAST
        if (prior.getNodeHeightPrior() == TreePriorType.EXTENDED_SKYLINE) { 
        	
        	writer.writeIDref(CoalescentLikelihood.COALESCENT_LIKELIHOOD, modelPrefix + COALESCENT); // only 1 coalescent
        	
            writer.writeOpenTag(MixedDistributionLikelihood.DISTRIBUTION_LIKELIHOOD);

            writer.writeOpenTag(MixedDistributionLikelihood.DISTRIBUTION0);
            writer.writeIDref(ExponentialDistributionModel.EXPONENTIAL_DISTRIBUTION_MODEL, modelPrefix + "demographic.populationMeanDist");
            writer.writeCloseTag(MixedDistributionLikelihood.DISTRIBUTION0);

            writer.writeOpenTag(MixedDistributionLikelihood.DISTRIBUTION1);
            writer.writeIDref(ExponentialDistributionModel.EXPONENTIAL_DISTRIBUTION_MODEL, modelPrefix + "demographic.populationMeanDist");
            writer.writeCloseTag(MixedDistributionLikelihood.DISTRIBUTION1);

            writer.writeOpenTag(MixedDistributionLikelihood.DATA);

            writer.writeIDref(ParameterParser.PARAMETER, modelPrefix + "demographic.popSize");

            writer.writeCloseTag(MixedDistributionLikelihood.DATA);

            writer.writeOpenTag(MixedDistributionLikelihood.INDICATORS);

            writer.writeIDref(ParameterParser.PARAMETER, modelPrefix + "demographic.indicators");

            writer.writeCloseTag(MixedDistributionLikelihood.INDICATORS);

            writer.writeCloseTag(MixedDistributionLikelihood.DISTRIBUTION_LIKELIHOOD);
        }

    }
}
