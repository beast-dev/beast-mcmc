/*
 * BeastParser.java
 *
 * Copyright (C) 2002-2006 Alexei Drummond and Andrew Rambaut
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

package dr.app.beast;

import dr.xml.PropertyParser;
import dr.xml.UserInput;
import dr.xml.XMLObjectParser;
import dr.xml.XMLParser;

import java.util.*;
import java.io.*;
import java.lang.reflect.Field;

/**
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @version $Id: BeastParser.java,v 1.76 2006/08/30 16:01:59 rambaut Exp $
 */
public class BeastParser extends XMLParser {

    public BeastParser(String[] args, List<String> additionalParsers, boolean verbose, boolean strictXML) {
        super(strictXML);

        setup(args);

        if (verbose) {
            System.out.println("Built-in parsers:");
            Iterator iterator = getParsers();
            while (iterator.hasNext()) {
                XMLObjectParser parser = (XMLObjectParser) iterator.next();
                System.out.println(parser.getParserName());
            }

        }

        // Try to find and load the additional 'core' parsers
        try {
            Properties properties = new Properties();
            properties.load(this.getClass().getResourceAsStream("beast.properties"));

            // get the parsers file prefix from the beast.properties file
            String parsers = properties.getProperty("parsers");

            if (System.getProperty("parsers") != null) {
                // If a system property has been set then allow this to override the default
                // e.g. -Dparsers=development
                parsers = properties.getProperty("parsers");
            }

            if (parsers != null) {
                // load the parsers
                loadProperties(this.getClass(), parsers +  "_parsers.properties", verbose);
            }

            if (additionalParsers != null) {
                for (String addParsers : additionalParsers) {
                    loadProperties(this.getClass(), addParsers +  "_parsers.properties", verbose);
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
        }

        // Now search the package hierarchy for 'beast.properties' files.
//        try {
//            loadProperties(this.getClass(), verbose);
//        } catch (IOException e) {
//            e.printStackTrace();
//        }
    }

    private void loadProperties(Class c, String parsersFile, boolean verbose) throws IOException {

        if (verbose) {
            System.out.println();
            System.out.println("Loading additional parsers (" + parsersFile + "):");
        }
        final InputStream stream = c.getResourceAsStream(parsersFile);
        if (stream == null) {
            throw new RuntimeException("Parsers file not found: " + parsersFile);
        }
        BufferedReader reader = new BufferedReader(new InputStreamReader(stream));
        String line = reader.readLine();
        while (line != null) {

            if (line.trim().length() > 0 && !line.trim().startsWith("#")) {
                try {
                    if (line.contains("Vector")) {
                        System.out.println("");
                    }
                    Class parser = Class.forName(line);
                    if (XMLObjectParser.class.isAssignableFrom(parser)) {
                        // if this class is an XMLObjectParser then create an instance
                        boolean replaced = addXMLObjectParser((XMLObjectParser)parser.newInstance(), true);
                        if (verbose) {
                            System.out.println((replaced ? "Replaced" : "Loaded") + " parser: " + parser.getName());
                        }
                    } else {
                        boolean parserFound = false;
                        // otherwise look for a static member which is an instance of XMLObjectParser
                        Field[] fields = parser.getDeclaredFields();
                        for (Field field: fields) {
                            if (XMLObjectParser.class.isAssignableFrom(field.getType())) {
                                try {
                                    boolean replaced = addXMLObjectParser((XMLObjectParser)field.get(null), true);
                                    if (verbose) {
                                        System.out.println((replaced ? "Replaced" : "Loaded") + " parser: " + parser.getName() + "." + field.getName());
                                    }
                                } catch (IllegalArgumentException iae) {
                                    System.err.println("Failed to install parser: " + iae.getMessage());
                                }
                                parserFound = true;
                            }
                        }

                        if (!parserFound) {
                            throw new IllegalArgumentException(parser.getName() + " is not of type XMLObjectParser " +
                                    "and doesn't contain any static members of this type");
                        }
                    }

                } catch (Exception e) {
                    System.err.println("Failed to load parser: " + e.getMessage());
                }
            }
            line = reader.readLine();

        }
    }

    private void setup(String[] args) {

        for (int i = 0; i < args.length; i++) {
            storeObject(Integer.toString(i), args[i]);
        }

        // add all the XMLObject parsers you need

        addXMLObjectParser(new PropertyParser());
        addXMLObjectParser(UserInput.STRING_PARSER);
        addXMLObjectParser(UserInput.DOUBLE_PARSER);
        addXMLObjectParser(UserInput.INTEGER_PARSER);

        //addXMLObjectParser(ColouringTest.PARSER);

        addXMLObjectParser(new dr.evoxml.GeneralDataTypeParser());
        addXMLObjectParser(new dr.evoxml.AlignmentParser());
        addXMLObjectParser(new dr.evoxml.SitePatternsParser());
        addXMLObjectParser(new dr.evoxml.PatternSubSetParser());
        addXMLObjectParser(new dr.evoxml.AscertainedSitePatternsParser());
        addXMLObjectParser(new dr.evoxml.ConvertAlignmentParser());
        addXMLObjectParser(new dr.evoxml.MergePatternsParser());
        addXMLObjectParser(new dr.evoxml.AttributePatternsParser());
        addXMLObjectParser(new dr.evoxml.SequenceParser());
        addXMLObjectParser(new dr.evoxml.SimpleNodeParser());
        addXMLObjectParser(new dr.evoxml.SimpleTreeParser());
        addXMLObjectParser(new dr.evoxml.UPGMATreeParser());
        addXMLObjectParser(new dr.evoxml.NeighborJoiningParser());
        addXMLObjectParser(new dr.evoxml.NewickParser());
        addXMLObjectParser(new dr.evoxml.TaxonParser());
        addXMLObjectParser(new dr.evoxml.TaxaParser());
        addXMLObjectParser(new dr.evoxml.DateParser());
        addXMLObjectParser(new dr.evoxml.LocationParser());
        addXMLObjectParser(new dr.evoxml.DistanceMatrixParser());
        addXMLObjectParser(new dr.evoxml.HiddenNucleotideParser());
        addXMLObjectParser(dr.evolution.util.RandomTaxaSample.PARSER);

        addXMLObjectParser(new dr.xml.AttributeParser());
        addXMLObjectParser(new dr.xml.AttributesParser());

        // speciation
        addXMLObjectParser(new dr.evomodelxml.YuleModelParser());
        addXMLObjectParser(new dr.evomodelxml.BirthDeathModelParser());
        addXMLObjectParser(dr.evomodel.speciation.SpeciationLikelihood.PARSER);

        // coalescent
        addXMLObjectParser(dr.evomodel.coalescent.CoalescentSimulator.PARSER);
        addXMLObjectParser(dr.evomodel.coalescent.CoalescentLikelihood.PARSER);
        addXMLObjectParser(dr.evomodel.coalescent.BayesianSkylineLikelihood.PARSER);
        addXMLObjectParser(dr.evomodel.coalescent.ConstantPopulationModel.PARSER);
        addXMLObjectParser(dr.evomodel.coalescent.ExponentialGrowthModel.PARSER);
        addXMLObjectParser(dr.evomodel.coalescent.LogisticGrowthModel.PARSER);
        addXMLObjectParser(dr.evomodel.coalescent.ConstantExponentialModel.PARSER);
        addXMLObjectParser(dr.evomodel.coalescent.ConstantLogisticModel.PARSER);
        addXMLObjectParser(dr.evomodel.coalescent.ExpansionModel.PARSER);

        addXMLObjectParser(dr.evomodel.coalescent.DemographicLogger.PARSER);

        // substitution models
        addXMLObjectParser(dr.evomodel.substmodel.FrequencyModel.PARSER);
        addXMLObjectParser(dr.evomodel.substmodel.GeneralSubstitutionModel.PARSER);
        addXMLObjectParser(new dr.evomodelxml.HKYParser());
        addXMLObjectParser(dr.evomodel.substmodel.GTR.PARSER);
        addXMLObjectParser(dr.evomodel.substmodel.EmpiricalAminoAcidModel.PARSER);
        addXMLObjectParser(dr.evomodel.substmodel.YangCodonModel.PARSER);

        // tree likelihood
        addXMLObjectParser(dr.evomodel.treelikelihood.TreeLikelihood.PARSER);

        // site models
        addXMLObjectParser(dr.evomodel.sitemodel.GammaSiteModel.PARSER);
        addXMLObjectParser(dr.evomodel.sitemodel.CategorySiteModel.PARSER);

        // molecular clocks
        addXMLObjectParser(dr.evomodel.clock.ACLikelihood.PARSER);
        addXMLObjectParser(dr.evomodel.clock.UCLikelihood.PARSER);
        addXMLObjectParser(new dr.evomodelxml.DiscretizedBranchRatesParser());
        addXMLObjectParser(dr.evomodel.branchratemodel.StrictClockBranchRates.PARSER);
        addXMLObjectParser(dr.evomodel.branchratemodel.RateEpochBranchRateModel.PARSER);
        addXMLObjectParser(dr.evomodel.branchratemodel.RandomLocalClockModel.PARSER);

        // tree models
        addXMLObjectParser(new dr.evomodelxml.TreeModelParser());
        addXMLObjectParser(dr.evomodel.tree.TipHeightLikelihood.PARSER);
        addXMLObjectParser(dr.evomodel.tree.TreeMetricStatistic.PARSER);
        addXMLObjectParser(dr.evomodel.tree.TreelengthStatistic.PARSER);
        addXMLObjectParser(dr.evomodel.tree.TreeShapeStatistic.PARSER);
        addXMLObjectParser(dr.evomodel.tree.TMRCAStatistic.PARSER);
        addXMLObjectParser(dr.evomodel.tree.MRCATraitStatistic.PARSER);
        addXMLObjectParser(dr.evomodel.tree.ExternalLengthStatistic.PARSER);
        addXMLObjectParser(dr.evomodel.tree.RateCovarianceStatistic.PARSER);
        addXMLObjectParser(dr.evomodel.tree.RateStatistic.PARSER);
        addXMLObjectParser(dr.evomodel.tree.MonophylyStatistic.PARSER);
        addXMLObjectParser(dr.evomodel.tree.CompatibilityStatistic.PARSER);
        addXMLObjectParser(dr.evomodel.tree.ParsimonyStatistic.PARSER);
        addXMLObjectParser(dr.evomodel.tree.ParsimonyStateStatistic.PARSER);
        addXMLObjectParser(dr.evomodel.tree.SpeciesTreeStatistic.PARSER);
        addXMLObjectParser(dr.evomodel.tree.UniformRootPrior.PARSER);

        // tree operators
        addXMLObjectParser(dr.evomodel.operators.SubtreeSlideOperator.PARSER);
        addXMLObjectParser(dr.evomodel.operators.ExchangeOperator.NARROW_EXCHANGE_PARSER);
        addXMLObjectParser(dr.evomodel.operators.ExchangeOperator.WIDE_EXCHANGE_PARSER);
        addXMLObjectParser(dr.evomodel.operators.NNI.NNI_PARSER);
        addXMLObjectParser(dr.evomodel.operators.FNPR.FNPR_PARSER);
        addXMLObjectParser(dr.evomodel.operators.WilsonBalding.PARSER);
        addXMLObjectParser(dr.evomodel.operators.ImportancePruneAndRegraft.IMPORTANCE_PRUNE_AND_REGRAFT_PARSER);
        addXMLObjectParser(dr.evomodel.operators.ImportanceSubtreeSwap.IMPORTANCE_SUBTREE_SWAP_PARSER);
        addXMLObjectParser(dr.evomodel.operators.GibbsSubtreeSwap.GIBBS_SUBTREE_SWAP_PARSER);
        addXMLObjectParser(dr.evomodel.operators.GibbsPruneAndRegraft.GIBBS_PRUNE_AND_REGRAFT_PARSER);
        addXMLObjectParser(dr.evomodel.operators.RateExchangeOperator.PARSER);
        addXMLObjectParser(dr.evomodel.operators.TreeBitMoveOperator.PARSER);
        addXMLObjectParser(dr.evomodel.operators.TreeBitRandomWalkOperator.PARSER);
  //      addXMLObjectParser(dr.evomodel.operators.TreeUniform.PARSSR);
  //      addXMLObjectParser(dr.evomodel.operators.INS.INS_PARSER);

        // rate operators
        addXMLObjectParser(dr.evomodel.operators.RateScaleOperator.PARSER);
        addXMLObjectParser(dr.evomodel.operators.RateVarianceScaleOperator.PARSER);
        addXMLObjectParser(dr.evomodel.operators.RateSampleOperator.PARSER);

        // likelihoods, models and distributions
        addXMLObjectParser(dr.inference.model.CompoundParameter.PARSER);
        addXMLObjectParser(dr.inference.model.CompoundLikelihood.PARSER);
        addXMLObjectParser(dr.inference.model.BooleanLikelihood.PARSER);
        addXMLObjectParser(dr.inference.model.DummyLikelihood.PARSER);
        addXMLObjectParser(dr.inference.model.JeffreysPriorLikelihood.PARSER);

        addXMLObjectParser(dr.inference.model.OrnsteinUhlenbeckPriorLikelihood.PARSER);

        addXMLObjectParser(dr.inference.distribution.DistributionLikelihood.PARSER);
        addXMLObjectParser(dr.inference.distribution.MixedDistributionLikelihood.PARSER);
        addXMLObjectParser(dr.inference.distribution.UniformDistributionModel.PARSER);
        addXMLObjectParser(dr.inference.distribution.ExponentialDistributionModel.PARSER);
        addXMLObjectParser(dr.inference.distribution.GammaDistributionModel.PARSER);
        addXMLObjectParser(dr.inference.distribution.NormalDistributionModel.PARSER);
        addXMLObjectParser(dr.inference.distribution.LogNormalDistributionModel.PARSER);
        addXMLObjectParser(dr.inference.distribution.ExponentialMarkovModel.PARSER);

        addXMLObjectParser(dr.inference.distribution.DistributionLikelihood.UNIFORM_PRIOR_PARSER);
        addXMLObjectParser(dr.inference.distribution.DistributionLikelihood.EXPONENTIAL_PRIOR_PARSER);
        addXMLObjectParser(dr.inference.distribution.DistributionLikelihood.NORMAL_PRIOR_PARSER);
        addXMLObjectParser(dr.inference.distribution.DistributionLikelihood.POISSON_PRIOR_PARSER);
        addXMLObjectParser(dr.inference.distribution.DistributionLikelihood.LOG_NORMAL_PRIOR_PARSER);
        addXMLObjectParser(dr.inference.distribution.DistributionLikelihood.GAMMA_PRIOR_PARSER);
  //      addXMLObjectParser(dr.inference.distribution.DistributionLikelihood.INVGAMMA_PRIOR_PARSER);

        addXMLObjectParser(dr.inference.distribution.BinomialLikelihood.PARSER);

        addXMLObjectParser(new dr.inference.model.StatisticParser());
        addXMLObjectParser(new dr.inference.model.ParameterParser());

        addXMLObjectParser(dr.inference.model.MeanStatistic.PARSER);
        addXMLObjectParser(dr.inference.model.VarianceStatistic.PARSER);
        addXMLObjectParser(dr.inference.model.ProductStatistic.PARSER);
        addXMLObjectParser(dr.inference.model.SumStatistic.PARSER);
        addXMLObjectParser(dr.inference.model.DifferenceStatistic.PARSER);
        addXMLObjectParser(dr.inference.model.RatioStatistic.PARSER);
        addXMLObjectParser(dr.inference.model.ReciprocalStatistic.PARSER);
        addXMLObjectParser(dr.inference.model.NegativeStatistic.PARSER);
        addXMLObjectParser(dr.inference.model.ExponentialStatistic.PARSER);
        addXMLObjectParser(dr.inference.model.LogarithmStatistic.PARSER);
        addXMLObjectParser(dr.inference.model.ExpressionStatistic.PARSER);

        addXMLObjectParser(dr.inference.model.RPNcalculatorStatistic.PARSER);

        addXMLObjectParser(dr.inference.model.TestStatistic.PARSER);
        addXMLObjectParser(dr.inference.model.NotStatistic.PARSER);
        addXMLObjectParser(dr.inference.model.SubStatistic.PARSER);

        // Markov chains and loggers
        addXMLObjectParser(dr.inference.mcmc.MCMC.PARSER);
        addXMLObjectParser(dr.inference.ml.MLOptimizer.PARSER);

        addXMLObjectParser(new dr.evomodelxml.LoggerParser());
        addXMLObjectParser(dr.inference.loggers.MLLogger.ML_LOGGER_PARSER);
        addXMLObjectParser(new dr.evomodelxml.TreeLoggerParser());
        addXMLObjectParser(dr.inference.loggers.Columns.PARSER);
        addXMLObjectParser(dr.inference.operators.SimpleOperatorSchedule.PARSER);
//		addXMLObjectParser(new dr.inference.markovchain.ConvergenceListenerParser());

        // operators
        addXMLObjectParser(dr.inference.operators.RandomWalkIntegerOperator.PARSER);
        addXMLObjectParser(dr.inference.operators.RandomWalkOperator.PARSER);
        addXMLObjectParser(dr.inference.operators.ScaleOperator.PARSER);
        addXMLObjectParser(dr.inference.operators.LogRandomWalkOperator.PARSER);
        addXMLObjectParser(dr.inference.operators.UniformOperator.PARSER);
        addXMLObjectParser(dr.inference.operators.UniformIntegerOperator.PARSER);
        addXMLObjectParser(dr.inference.operators.UpDownOperator.PARSER);
        addXMLObjectParser(dr.inference.operators.SetOperator.PARSER);
        addXMLObjectParser(dr.inference.operators.SwapOperator.PARSER);
        addXMLObjectParser(dr.inference.operators.DeltaExchangeOperator.PARSER);
        addXMLObjectParser(dr.inference.operators.CenteredScaleOperator.PARSER);
        addXMLObjectParser(dr.inference.operators.BitFlipOperator.PARSER);
        addXMLObjectParser(dr.inference.operators.BitMoveOperator.PARSER);
        addXMLObjectParser(dr.inference.operators.BitSwapOperator.PARSER);
        addXMLObjectParser(dr.inference.operators.JointOperator.PARSER);
   //     addXMLObjectParser(dr.inference.operators.TeamOperator.PARSER);

        // trace analysis
        addXMLObjectParser(new dr.evomodel.tree.TreeTraceAnalysisParser());
        addXMLObjectParser(new dr.inference.trace.TraceAnalysisParser());
        addXMLObjectParser(dr.inference.trace.LogFileTraceExporter.PARSER);
        addXMLObjectParser(new dr.exporters.CSVExporter());

        addXMLObjectParser(dr.inference.trace.MarginalLikelihoodAnalysis.PARSER);


        addXMLObjectParser(dr.inference.model.ThreadedCompoundLikelihood.PARSER);

    }

}

