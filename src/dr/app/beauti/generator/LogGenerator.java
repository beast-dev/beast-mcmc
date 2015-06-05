/*
 * LogGenerator.java
 *
 * Copyright (c) 2002-2013 Alexei Drummond, Andrew Rambaut and Marc Suchard
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
import dr.app.beauti.types.FixRateType;
import dr.app.beauti.types.TreePriorType;
import dr.app.beauti.util.XMLWriter;
import dr.evolution.datatype.DataType;
import dr.evolution.util.Taxa;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.tree.TreeModel;
import dr.evomodelxml.branchratemodel.DiscretizedBranchRatesParser;
import dr.evomodelxml.branchratemodel.LocalClockModelParser;
import dr.evomodelxml.branchratemodel.RandomLocalClockModelParser;
import dr.evomodelxml.branchratemodel.StrictClockBranchRatesParser;
import dr.evomodelxml.clock.ACLikelihoodParser;
import dr.evomodelxml.coalescent.CoalescentLikelihoodParser;
import dr.evomodelxml.coalescent.GMRFSkyrideLikelihoodParser;
import dr.evomodelxml.speciation.*;
import dr.evomodelxml.tree.TMRCAStatisticParser;
import dr.evomodelxml.tree.TreeLoggerParser;
import dr.evomodelxml.tree.TreeModelParser;
import dr.inference.model.ParameterParser;
import dr.inferencexml.distribution.MixedDistributionLikelihoodParser;
import dr.inferencexml.loggers.ColumnsParser;
import dr.inferencexml.loggers.LoggerParser;
import dr.inferencexml.model.CompoundLikelihoodParser;
import dr.inferencexml.model.CompoundParameterParser;
import dr.util.Attribute;
import dr.xml.XMLParser;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @author Walter Xie
 */
public class LogGenerator extends Generator {

    private final static String TREE_FILE_LOG = "treeFileLog";
    private final static String SUB_TREE_FILE_LOG = "substTreeFileLog";

    public LogGenerator(BeautiOptions options, ComponentFactory[] components) {
        super(options, components);
    }

    /**
     * write log to screen
     *
     * @param writer                    XMLWriter
     * @param branchRatesModelGenerator BranchRatesModelGenerator
     */
    public void writeLogToScreen(XMLWriter writer, BranchRatesModelGenerator branchRatesModelGenerator,
                                 SubstitutionModelGenerator substitutionModelGenerator) {
        writer.writeComment("write log to screen");

        writer.writeOpenTag(LoggerParser.LOG,
                new Attribute[]{
                        new Attribute.Default<String>(XMLParser.ID, "screenLog"),
                        new Attribute.Default<String>(LoggerParser.LOG_EVERY, options.echoEvery + "")
                });

        if (options.hasData()) {
            writer.writeOpenTag(ColumnsParser.COLUMN,
                    new Attribute[]{
                            new Attribute.Default<String>(ColumnsParser.LABEL, "Posterior"),
                            new Attribute.Default<String>(ColumnsParser.DECIMAL_PLACES, "4"),
                            new Attribute.Default<String>(ColumnsParser.WIDTH, "12")
                    }
            );
            writer.writeIDref(CompoundLikelihoodParser.POSTERIOR, "posterior");
            writer.writeCloseTag(ColumnsParser.COLUMN);
        }

        writer.writeOpenTag(ColumnsParser.COLUMN,
                new Attribute[]{
                        new Attribute.Default<String>(ColumnsParser.LABEL, "Prior"),
                        new Attribute.Default<String>(ColumnsParser.DECIMAL_PLACES, "4"),
                        new Attribute.Default<String>(ColumnsParser.WIDTH, "12")
                }
        );
        writer.writeIDref(CompoundLikelihoodParser.PRIOR, "prior");
        writer.writeCloseTag(ColumnsParser.COLUMN);

        if (options.hasData()) {
            writer.writeOpenTag(ColumnsParser.COLUMN,
                    new Attribute[]{
                            new Attribute.Default<String>(ColumnsParser.LABEL, "Likelihood"),
                            new Attribute.Default<String>(ColumnsParser.DECIMAL_PLACES, "4"),
                            new Attribute.Default<String>(ColumnsParser.WIDTH, "12")
                    }
            );
            writer.writeIDref(CompoundLikelihoodParser.LIKELIHOOD, "likelihood");
            writer.writeCloseTag(ColumnsParser.COLUMN);
        }

        if (options.useStarBEAST) { // species
            writer.writeOpenTag(ColumnsParser.COLUMN,
                    new Attribute[]{
                            new Attribute.Default<String>(ColumnsParser.LABEL, "PopMean"),
                            new Attribute.Default<String>(ColumnsParser.DECIMAL_PLACES, "4"),
                            new Attribute.Default<String>(ColumnsParser.WIDTH, "12")
                    }
            );
            writer.writeIDref(ParameterParser.PARAMETER, TraitData.TRAIT_SPECIES + "." + options.starBEASTOptions.POP_MEAN);
            writer.writeCloseTag(ColumnsParser.COLUMN);
        }

        for (PartitionTreeModel model : options.getPartitionTreeModels()) {
            writer.writeOpenTag(ColumnsParser.COLUMN,
                    new Attribute[]{
                            new Attribute.Default<String>(ColumnsParser.LABEL, model.getPrefix() + TreeModelParser.ROOT_HEIGHT),
                            new Attribute.Default<String>(ColumnsParser.SIGNIFICANT_FIGURES, "6"),
                            new Attribute.Default<String>(ColumnsParser.WIDTH, "12")
                    }
            );

            writer.writeIDref(ParameterParser.PARAMETER, model.getPrefix() + TreeModel.TREE_MODEL + "." + TreeModelParser.ROOT_HEIGHT);

            writer.writeCloseTag(ColumnsParser.COLUMN);
        }

        for (PartitionClockModel model : options.getPartitionClockModels()) {
            writer.writeOpenTag(ColumnsParser.COLUMN,
                    new Attribute[]{
                            new Attribute.Default<String>(ColumnsParser.LABEL, branchRatesModelGenerator.getClockRateString(model)),
                            new Attribute.Default<String>(ColumnsParser.SIGNIFICANT_FIGURES, "6"),
                            new Attribute.Default<String>(ColumnsParser.WIDTH, "12")
                    }
            );

            branchRatesModelGenerator.writeAllClockRateRefs(model, writer);
//        if (options.clockModelOptions.getRateOptionClockModel() == FixRateType.FIX_MEAN) {
//            writer.writeIDref(ParameterParser.PARAMETER, "allClockRates");
//            for (PartitionClockModel model : options.getPartitionClockModels()) {
//                if (model.getClockType() == ClockType.UNCORRELATED_LOGNORMAL)
//                    writer.writeIDref(ParameterParser.PARAMETER, model.getPrefix() + ClockType.UCLD_STDEV);
//            }
//        } else {
//            for (PartitionClockModel model : options.getPartitionClockModels()) {
//                branchRatesModelGenerator.writeAllClockRateRefs(model, writer);
//            }
//        }
            writer.writeCloseTag(ColumnsParser.COLUMN);
        }

        for (PartitionSubstitutionModel model : options.getPartitionSubstitutionModels()) {
            if (model.getDataType().getType() == DataType.MICRO_SAT)
                substitutionModelGenerator.writeMicrosatSubstModelParameterRef(model, writer);
        }

        generateInsertionPoint(ComponentGenerator.InsertionPoint.IN_SCREEN_LOG, writer);

        writer.writeCloseTag(LoggerParser.LOG);

        generateInsertionPoint(ComponentGenerator.InsertionPoint.AFTER_SCREEN_LOG, writer);
    }

    /**
     * write log to file
     *
     * @param writer                     XMLWriter
     * @param treePriorGenerator         TreePriorGenerator
     * @param branchRatesModelGenerator  BranchRatesModelGenerator
     * @param substitutionModelGenerator SubstitutionModelGenerator
     * @param treeLikelihoodGenerator    TreeLikelihoodGenerator
     */
    public void writeLogToFile(XMLWriter writer,
                               TreePriorGenerator treePriorGenerator,
                               BranchRatesModelGenerator branchRatesModelGenerator,
                               SubstitutionModelGenerator substitutionModelGenerator,
                               TreeLikelihoodGenerator treeLikelihoodGenerator) {
        writer.writeComment("write log to file");

        if (options.logFileName == null) {
            options.logFileName = options.fileNameStem + ".log";
        }
        writer.writeOpenTag(LoggerParser.LOG,
                new Attribute[]{
                        new Attribute.Default<String>(XMLParser.ID, "fileLog"),
                        new Attribute.Default<String>(LoggerParser.LOG_EVERY, options.logEvery + ""),
                        new Attribute.Default<String>(LoggerParser.FILE_NAME, options.logFileName),
                        new Attribute.Default<Boolean>(LoggerParser.ALLOW_OVERWRITE_LOG, options.allowOverwriteLog)
                });

        if (options.hasData()) {
            writer.writeIDref(CompoundLikelihoodParser.POSTERIOR, "posterior");
        }
        writer.writeIDref(CompoundLikelihoodParser.PRIOR, "prior");
        if (options.hasData()) {
            writer.writeIDref(CompoundLikelihoodParser.LIKELIHOOD, "likelihood");
        }

        if (options.useStarBEAST) { // species
            // coalescent prior
            writer.writeIDref(MultiSpeciesCoalescentParser.SPECIES_COALESCENT, TraitData.TRAIT_SPECIES + "." + COALESCENT);
            // prior on population sizes
//            if (options.speciesTreePrior == TreePriorType.SPECIES_YULE) {
            writer.writeIDref(MixedDistributionLikelihoodParser.DISTRIBUTION_LIKELIHOOD, SPOPS);
//            } else {
//                writer.writeIDref(SpeciesTreeBMPrior.STPRIOR, STP);
//            }
            // prior on species tree
            writer.writeIDref(SpeciationLikelihoodParser.SPECIATION_LIKELIHOOD, SPECIATION_LIKE);

            writer.writeIDref(ParameterParser.PARAMETER, TraitData.TRAIT_SPECIES + "." + options.starBEASTOptions.POP_MEAN);
            writer.writeIDref(ParameterParser.PARAMETER, SpeciesTreeModelParser.SPECIES_TREE + "." + SPLIT_POPS);

            if (options.getPartitionTreePriors().get(0).getNodeHeightPrior() == TreePriorType.SPECIES_BIRTH_DEATH) {
                writer.writeIDref(ParameterParser.PARAMETER, TraitData.TRAIT_SPECIES + "." + BirthDeathModelParser.MEAN_GROWTH_RATE_PARAM_NAME);
                writer.writeIDref(ParameterParser.PARAMETER, TraitData.TRAIT_SPECIES + "." + BirthDeathModelParser.RELATIVE_DEATH_RATE_PARAM_NAME);
            } else if (options.getPartitionTreePriors().get(0).getNodeHeightPrior() == TreePriorType.SPECIES_YULE ||
                    options.getPartitionTreePriors().get(0).getNodeHeightPrior() == TreePriorType.SPECIES_YULE_CALIBRATION) {
                writer.writeIDref(ParameterParser.PARAMETER, TraitData.TRAIT_SPECIES + "." + YuleModelParser.YULE + "." + YuleModelParser.BIRTH_RATE);
            } else {
                throw new IllegalArgumentException("Get wrong species tree prior using *BEAST : " + options.getPartitionTreePriors().get(0).getNodeHeightPrior().toString());
            }

            //Species Tree: tmrcaStatistic
            writer.writeIDref(TMRCAStatisticParser.TMRCA_STATISTIC, SpeciesTreeModelParser.SPECIES_TREE + "." + TreeModelParser.ROOT_HEIGHT);
        }

        for (PartitionTreeModel model : options.getPartitionTreeModels()) {
            writer.writeIDref(ParameterParser.PARAMETER, model.getPrefix() + TreeModel.TREE_MODEL + "." + TreeModelParser.ROOT_HEIGHT);
        }

        if (options.useStarBEAST) {
            for (Taxa taxa : options.speciesSets) {
                // make tmrca(tree.name) eay to read in log for Tracer
                writer.writeIDref(TMRCAStatisticParser.TMRCA_STATISTIC, "tmrca(" + taxa.getId() + ")");
            }
        } else {
            for (Taxa taxa : options.taxonSets) {
                // make tmrca(tree.name) eay to read in log for Tracer
                PartitionTreeModel treeModel = options.taxonSetsTreeModel.get(taxa);
                writer.writeIDref(TMRCAStatisticParser.TMRCA_STATISTIC, "tmrca(" + treeModel.getPrefix() + taxa.getId() + ")");
            }
        }

//        if ( options.shareSameTreePrior ) { // Share Same Tree Prior
//	        treePriorGenerator.setModelPrefix("");
//        	treePriorGenerator.writeParameterLog(options.activedSameTreePrior, writer);
//        } else { // no species
        for (PartitionTreePrior prior : options.getPartitionTreePriors()) {
//	        	treePriorGenerator.setModelPrefix(prior.getPrefix()); // priorName.treeModel
            treePriorGenerator.writeParameterLog(prior, writer);
        }
//	    }

        for (PartitionSubstitutionModel model : options.getPartitionSubstitutionModels()) {
            substitutionModelGenerator.writeLog(model, writer);
            if (model.hasCodon()) {
                writer.writeIDref(CompoundParameterParser.COMPOUND_PARAMETER, model.getPrefix() + "allMus");
            }
        }

        for (ClockModelGroup clockModelGroup : options.clockModelOptions.getClockModelGroups()) {
            if (clockModelGroup.getRateTypeOption() == FixRateType.FIX_MEAN) {
                writer.writeIDref(ParameterParser.PARAMETER, clockModelGroup.getName());
                for (PartitionClockModel model : options.getPartitionClockModels(clockModelGroup)) {
                    if (model.getClockType() == ClockType.UNCORRELATED) {
                        switch (model.getClockDistributionType()) {
                            case LOGNORMAL:
                                writer.writeIDref(ParameterParser.PARAMETER, model.getPrefix() + ClockType.UCLD_STDEV);
                                break;
                            case GAMMA:
                                throw new UnsupportedOperationException("Uncorrelated gamma model not implemented yet");
//                            break;
                            case CAUCHY:
                                throw new UnsupportedOperationException("Uncorrelated Cauchy model not implemented yet");
//                            break;
                            case EXPONENTIAL:
                                // nothing required
                                break;
                        }
                    }
                }
            } else {
                for (PartitionClockModel model : options.getPartitionClockModels(clockModelGroup)) {
                    branchRatesModelGenerator.writeLog(model, writer);
                }
            }
        }

        for (PartitionClockModel model : options.getPartitionClockModels()) {
            branchRatesModelGenerator.writeLogStatistic(model, writer);
        }

        if (options.logCoalescentEventsStatistic) {
            writer.writeOpenTag("coalescentEventsStatistic");
            // coalescentLikelihood
            for (PartitionTreeModel model : options.getPartitionTreeModels()) {
                PartitionTreePrior prior = model.getPartitionTreePrior();
                treePriorGenerator.writePriorLikelihoodReferenceLog(prior, model, writer);
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

        generateInsertionPoint(ComponentGenerator.InsertionPoint.IN_FILE_LOG_PARAMETERS, writer);

        treeLikelihoodGenerator.writeTreeLikelihoodReferences(writer);
        branchRatesModelGenerator.writeClockLikelihoodReferences(writer);

        generateInsertionPoint(ComponentGenerator.InsertionPoint.IN_FILE_LOG_LIKELIHOODS, writer);

        // coalescentLikelihood
        for (PartitionTreeModel model : options.getPartitionTreeModels()) {
            PartitionTreePrior prior = model.getPartitionTreePrior();
            treePriorGenerator.writePriorLikelihoodReferenceLog(prior, model, writer);
            writer.writeText("");
        }

        for (PartitionTreePrior prior : options.getPartitionTreePriors()) {
            if (prior.getNodeHeightPrior() == TreePriorType.EXTENDED_SKYLINE) {
                writer.writeIDref(CoalescentLikelihoodParser.COALESCENT_LIKELIHOOD, prior.getPrefix() + COALESCENT); // only 1 coalescent
            } else if (prior.getNodeHeightPrior() == TreePriorType.SKYGRID) {
                writer.writeIDref(GMRFSkyrideLikelihoodParser.SKYGRID_LIKELIHOOD, prior.getPrefix() + "skygrid");
            }
        }

        writer.writeCloseTag(LoggerParser.LOG);

        generateInsertionPoint(ComponentGenerator.InsertionPoint.AFTER_FILE_LOG, writer);

    }

    public void writeDemographicLogToFile(XMLWriter writer,
                               TreePriorGenerator treePriorGenerator,
                               BranchRatesModelGenerator branchRatesModelGenerator,
                               SubstitutionModelGenerator substitutionModelGenerator,
                               TreeLikelihoodGenerator treeLikelihoodGenerator) {
        writer.writeComment("demographic log file");

        if (options.demographicLogFileName == null) {
            options.demographicLogFileName = options.fileNameStem + ".demo.log";
        }

        String header = "Demographic Model: " + options.demographicModelName;
        writer.writeOpenTag(LoggerParser.LOG,
                new Attribute[]{
                        new Attribute.Default<String>(XMLParser.ID, "fileLog"),
                        new Attribute.Default<String>(LoggerParser.HEADER, header + ""),
                        new Attribute.Default<String>(LoggerParser.LOG_EVERY, options.logEvery + ""),
                        new Attribute.Default<String>(LoggerParser.FILE_NAME, options.logFileName),
                        new Attribute.Default<Boolean>(LoggerParser.ALLOW_OVERWRITE_LOG, options.allowOverwriteLog)
                });

        if (options.hasData()) {
            writer.writeIDref(CompoundLikelihoodParser.POSTERIOR, "posterior");
        }
        writer.writeIDref(CompoundLikelihoodParser.PRIOR, "prior");

        for (PartitionTreeModel model : options.getPartitionTreeModels()) {
            writer.writeIDref(ParameterParser.PARAMETER, model.getPrefix() + TreeModel.TREE_MODEL + "." + TreeModelParser.ROOT_HEIGHT);
        }

        if (options.useStarBEAST) {
            for (Taxa taxa : options.speciesSets) {
                // make tmrca(tree.name) eay to read in log for Tracer
                writer.writeIDref(TMRCAStatisticParser.TMRCA_STATISTIC, "tmrca(" + taxa.getId() + ")");
            }
        } else {
            for (Taxa taxa : options.taxonSets) {
                // make tmrca(tree.name) eay to read in log for Tracer
                PartitionTreeModel treeModel = options.taxonSetsTreeModel.get(taxa);
                writer.writeIDref(TMRCAStatisticParser.TMRCA_STATISTIC, "tmrca(" + treeModel.getPrefix() + taxa.getId() + ")");
            }
        }

//        if ( options.shareSameTreePrior ) { // Share Same Tree Prior
//	        treePriorGenerator.setModelPrefix("");
//        	treePriorGenerator.writeParameterLog(options.activedSameTreePrior, writer);
//        } else { // no species
        for (PartitionTreePrior prior : options.getPartitionTreePriors()) {
//	        	treePriorGenerator.setModelPrefix(prior.getPrefix()); // priorName.treeModel
            treePriorGenerator.writeParameterLog(prior, writer);
        }
//	    }

        for (PartitionSubstitutionModel model : options.getPartitionSubstitutionModels()) {
            substitutionModelGenerator.writeLog(model, writer);
            if (model.hasCodon()) {
                writer.writeIDref(CompoundParameterParser.COMPOUND_PARAMETER, model.getPrefix() + "allMus");
            }
        }

        for (ClockModelGroup clockModelGroup : options.clockModelOptions.getClockModelGroups()) {
            if (clockModelGroup.getRateTypeOption() == FixRateType.FIX_MEAN) {
                writer.writeIDref(ParameterParser.PARAMETER, clockModelGroup.getName());
                for (PartitionClockModel model : options.getPartitionClockModels(clockModelGroup)) {
                    if (model.getClockType() == ClockType.UNCORRELATED) {
                        switch (model.getClockDistributionType()) {
                            case LOGNORMAL:
                                writer.writeIDref(ParameterParser.PARAMETER, model.getPrefix() + ClockType.UCLD_STDEV);
                                break;
                            case GAMMA:
                                throw new UnsupportedOperationException("Uncorrelated gamma model not implemented yet");
//                            break;
                            case CAUCHY:
                                throw new UnsupportedOperationException("Uncorrelated Cauchy model not implemented yet");
//                            break;
                            case EXPONENTIAL:
                                // nothing required
                                break;
                        }
                    }
                }
            } else {
                for (PartitionClockModel model : options.getPartitionClockModels(clockModelGroup)) {
                    branchRatesModelGenerator.writeLog(model, writer);
                }
            }
        }

        for (PartitionClockModel model : options.getPartitionClockModels()) {
            branchRatesModelGenerator.writeLogStatistic(model, writer);
        }

        generateInsertionPoint(ComponentGenerator.InsertionPoint.IN_FILE_LOG_PARAMETERS, writer);

        treeLikelihoodGenerator.writeTreeLikelihoodReferences(writer);
        branchRatesModelGenerator.writeClockLikelihoodReferences(writer);

        generateInsertionPoint(ComponentGenerator.InsertionPoint.IN_FILE_LOG_LIKELIHOODS, writer);

        // coalescentLikelihood
        for (PartitionTreeModel model : options.getPartitionTreeModels()) {
            PartitionTreePrior prior = model.getPartitionTreePrior();
            treePriorGenerator.writePriorLikelihoodReferenceLog(prior, model, writer);
            writer.writeText("");
        }

        for (PartitionTreePrior prior : options.getPartitionTreePriors()) {
            if (prior.getNodeHeightPrior() == TreePriorType.EXTENDED_SKYLINE) {
                writer.writeIDref(CoalescentLikelihoodParser.COALESCENT_LIKELIHOOD, prior.getPrefix() + COALESCENT); // only 1 coalescent
            } else if (prior.getNodeHeightPrior() == TreePriorType.SKYGRID) {
                writer.writeIDref(GMRFSkyrideLikelihoodParser.SKYGRID_LIKELIHOOD, prior.getPrefix() + "skygrid");
            }
        }

        writer.writeCloseTag(LoggerParser.LOG);

        generateInsertionPoint(ComponentGenerator.InsertionPoint.AFTER_FILE_LOG, writer);

    }

    /**
     * write tree log to file
     *
     * @param writer XMLWriter
     */
    public void writeTreeLogToFile(XMLWriter writer) {
        writer.writeComment("write tree log to file");

        if (options.useStarBEAST) { // species
            // species tree log
            writer.writeOpenTag(TreeLoggerParser.LOG_TREE,
                    new Attribute[]{
                            new Attribute.Default<String>(XMLParser.ID, TraitData.TRAIT_SPECIES + "." + TREE_FILE_LOG), // speciesTreeFileLog
                            new Attribute.Default<String>(TreeLoggerParser.LOG_EVERY, options.logEvery + ""),
                            new Attribute.Default<String>(TreeLoggerParser.NEXUS_FORMAT, "true"),
                            new Attribute.Default<String>(TreeLoggerParser.FILE_NAME, options.fileNameStem + "." + options.starBEASTOptions.SPECIES_TREE_FILE_NAME),
                            new Attribute.Default<String>(TreeLoggerParser.SORT_TRANSLATION_TABLE, "true")
                    });

            writer.writeIDref(SpeciesTreeModelParser.SPECIES_TREE, SP_TREE);

            if (options.hasData()) {
                // we have data...
                writer.writeIDref("posterior", "posterior");
            }
            writer.writeCloseTag(TreeLoggerParser.LOG_TREE);
        }

        // gene tree log
        //TODO make code consistent to MCMCPanel
        for (PartitionTreeModel tree : options.getPartitionTreeModels()) {
            String treeFileName;
            if (options.substTreeLog) {
                treeFileName = options.fileNameStem + "." + tree.getPrefix() + "(time)." + STARBEASTOptions.TREE_FILE_NAME;
            } else {
                treeFileName = options.fileNameStem + "." + tree.getPrefix() + STARBEASTOptions.TREE_FILE_NAME; // stem.partitionName.tree
            }

            if (options.treeFileName.get(0).endsWith(".txt")) {
                treeFileName += ".txt";
            }

            List<Attribute> attributes = new ArrayList<Attribute>();

            attributes.add(new Attribute.Default<String>(XMLParser.ID, tree.getPrefix() + TREE_FILE_LOG)); // partionName.treeFileLog
            attributes.add(new Attribute.Default<String>(TreeLoggerParser.LOG_EVERY, options.logEvery + ""));
            attributes.add(new Attribute.Default<String>(TreeLoggerParser.NEXUS_FORMAT, "true"));
            attributes.add(new Attribute.Default<String>(TreeLoggerParser.FILE_NAME, treeFileName));
            attributes.add(new Attribute.Default<String>(TreeLoggerParser.SORT_TRANSLATION_TABLE, "true"));

            //if (options.clockModelOptions.getRateOptionClockModel() == FixRateType.RElATIVE_TO && tree.containsUncorrelatedRelaxClock()) { //TODO: Sibon's discretized branch length stuff
            //    double aveFixedRate = options.clockModelOptions.getSelectedRate(options.getPartitionClockModels());
            //    attributes.add(new Attribute.Default<String>(TreeLoggerParser.NORMALISE_MEAN_RATE_TO, Double.toString(aveFixedRate)));
            //}

            // generate <logTree>
            writer.writeOpenTag(TreeLoggerParser.LOG_TREE, attributes);

//            writer.writeOpenTag(TreeLoggerParser.LOG_TREE,
//                    new Attribute[]{
//                            new Attribute.Default<String>(XMLParser.ID, tree.getPrefix() + TREE_FILE_LOG), // partionName.treeFileLog
//                            new Attribute.Default<String>(TreeLoggerParser.LOG_EVERY, options.logEvery + ""),
//                            new Attribute.Default<String>(TreeLoggerParser.NEXUS_FORMAT, "true"),
//                            new Attribute.Default<String>(TreeLoggerParser.FILE_NAME, treeFileName),
//                            new Attribute.Default<String>(TreeLoggerParser.SORT_TRANSLATION_TABLE, "true")
//                    });

            writer.writeIDref(TreeModel.TREE_MODEL, tree.getPrefix() + TreeModel.TREE_MODEL);

            writeTreeTraits(writer, tree);

            if (options.hasData()) {
                // we have data...
                writer.writeIDref("posterior", "posterior");
            }

            generateInsertionPoint(ComponentGenerator.InsertionPoint.IN_TREES_LOG, tree, writer);

            writer.writeCloseTag(TreeLoggerParser.LOG_TREE);
        } // end For loop

        if (options.substTreeLog) {
            if (options.useStarBEAST) { // species
                //TODO: species sub tree
            }

            // gene tree
            for (PartitionTreeModel tree : options.getPartitionTreeModels()) {
                // write tree log to file
                writer.writeOpenTag(TreeLoggerParser.LOG_TREE,
                        new Attribute[]{
                                new Attribute.Default<String>(XMLParser.ID, tree.getPrefix() + SUB_TREE_FILE_LOG),
                                new Attribute.Default<String>(TreeLoggerParser.LOG_EVERY, options.logEvery + ""),
                                new Attribute.Default<String>(TreeLoggerParser.NEXUS_FORMAT, "true"),
                                new Attribute.Default<String>(TreeLoggerParser.FILE_NAME, options.fileNameStem + "." + tree.getPrefix() +
                                        "(subst)." + STARBEASTOptions.TREE_FILE_NAME),
                                new Attribute.Default<String>(TreeLoggerParser.BRANCH_LENGTHS, TreeLoggerParser.SUBSTITUTIONS)
                        });
                writer.writeIDref(TreeModel.TREE_MODEL, tree.getPrefix() + TreeModel.TREE_MODEL);

                // assume the first clock model is the one... (not sure if this makes sense)
                PartitionClockModel model = options.getPartitionClockModels(options.getDataPartitions(tree)).get(0);
                String tag = "";
                String id = "";

                switch (model.getClockType()) {
                    case STRICT_CLOCK:
                        tag = StrictClockBranchRatesParser.STRICT_CLOCK_BRANCH_RATES;
                        id = model.getPrefix() + BranchRateModel.BRANCH_RATES;
                        break;

                    case UNCORRELATED:
                        tag = DiscretizedBranchRatesParser.DISCRETIZED_BRANCH_RATES;
                        id = model.getPrefix() + BranchRateModel.BRANCH_RATES;
                        break;

                    case RANDOM_LOCAL_CLOCK:
                        tag = RandomLocalClockModelParser.LOCAL_BRANCH_RATES;
                        id = model.getPrefix() + BranchRateModel.BRANCH_RATES;
                        break;

                    case FIXED_LOCAL_CLOCK:
                        tag = LocalClockModelParser.LOCAL_CLOCK_MODEL;
                        id = model.getPrefix() + BranchRateModel.BRANCH_RATES;
                        break;
                    case AUTOCORRELATED:
                        tag = ACLikelihoodParser.AC_LIKELIHOOD;
                        id =  options.noDuplicatedPrefix(model.getPrefix(), tree.getPrefix()) + BranchRateModel.BRANCH_RATES;
                        break;

                    default:
                        throw new IllegalArgumentException("Unknown clock model");
                }
                writer.writeIDref(tag, id);
                writeTreeTrait(writer, tag, id, BranchRateModel.RATE, model.getPrefix() + BranchRateModel.RATE);

                writer.writeCloseTag(TreeLoggerParser.LOG_TREE);
            }
        }

        generateInsertionPoint(ComponentGenerator.InsertionPoint.AFTER_TREES_LOG, writer);
    }

    private void writeTreeTraits(XMLWriter writer, PartitionTreeModel tree) {
        for (PartitionClockModel model : options.getPartitionClockModels(options.getDataPartitions(tree))) {

            switch (model.getClockType()) {
                case STRICT_CLOCK:
                    writeTreeTrait(writer, StrictClockBranchRatesParser.STRICT_CLOCK_BRANCH_RATES,
                            model.getPrefix() + BranchRateModel.BRANCH_RATES,
                            BranchRateModel.RATE, model.getPrefix() + BranchRateModel.RATE);
                    break;

                case UNCORRELATED:
                    writeTreeTrait(writer, DiscretizedBranchRatesParser.DISCRETIZED_BRANCH_RATES,
                            options.noDuplicatedPrefix(model.getPrefix(), tree.getPrefix()) + BranchRateModel.BRANCH_RATES,
                            BranchRateModel.RATE, model.getPrefix() + BranchRateModel.RATE);
                    break;

                case RANDOM_LOCAL_CLOCK:
                    writeTreeTrait(writer, RandomLocalClockModelParser.LOCAL_BRANCH_RATES,
                            model.getPrefix() + BranchRateModel.BRANCH_RATES,
                            BranchRateModel.RATE, model.getPrefix() + BranchRateModel.RATE);
                    break;

                case FIXED_LOCAL_CLOCK:
                    writeTreeTrait(writer, LocalClockModelParser.LOCAL_CLOCK_MODEL,
                            options.noDuplicatedPrefix(model.getPrefix(), tree.getPrefix()) + BranchRateModel.BRANCH_RATES,
                            BranchRateModel.RATE, model.getPrefix() + BranchRateModel.RATE);
                    break;

                case AUTOCORRELATED:
                    writer.writeIDref(ACLikelihoodParser.AC_LIKELIHOOD,
                            options.noDuplicatedPrefix(model.getPrefix(), tree.getPrefix()) + BranchRateModel.BRANCH_RATES);
                    writeTreeTrait(writer, ACLikelihoodParser.AC_LIKELIHOOD,
                            options.noDuplicatedPrefix(model.getPrefix(), tree.getPrefix()) + BranchRateModel.BRANCH_RATES,
                            BranchRateModel.RATE, model.getPrefix() + BranchRateModel.RATE);
                    break;

                default:
                    throw new IllegalArgumentException("Unknown clock model");
            }
        }
    }

    private void writeTreeTrait(XMLWriter writer, String treeTraitTag, String treeTraitID, String traitName, String traitTag) {

        writer.writeOpenTag(TreeLoggerParser.TREE_TRAIT,
                new Attribute[]{
                        new Attribute.Default<String>(TreeLoggerParser.NAME, traitName),
                        new Attribute.Default<String>(TreeLoggerParser.TAG, traitTag)
                });
        writer.writeIDref(treeTraitTag, treeTraitID);

        writer.writeCloseTag(TreeLoggerParser.TREE_TRAIT);

    }
}