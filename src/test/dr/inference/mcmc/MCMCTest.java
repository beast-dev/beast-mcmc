/*
 * MCMCTest.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
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
 *
 */

package test.dr.inference.mcmc;

import dr.evolution.alignment.SitePatterns;
import dr.evolution.datatype.Nucleotides;
import dr.evomodel.branchmodel.HomogeneousBranchModel;
import dr.evomodel.branchratemodel.DefaultBranchRateModel;
import dr.evomodel.operators.ExchangeOperator;
import dr.evomodel.operators.SubtreeSlideOperator;
import dr.evomodel.operators.WilsonBalding;
import dr.evomodel.siteratemodel.GammaSiteRateModel;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.evomodel.substmodel.nucleotide.HKY;
import dr.evomodel.tree.DefaultTreeModel;
import dr.evomodel.treedatalikelihood.BeagleDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.DataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.PreOrderSettings;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treelikelihood.PartialsRescalingScheme;
import dr.evomodelxml.treedatalikelihood.TreeDataLikelihoodParser;
import dr.oldevomodel.sitemodel.GammaSiteModel;
import dr.oldevomodel.treelikelihood.TreeLikelihood;
import dr.oldevomodelxml.sitemodel.GammaSiteModelParser;
import dr.oldevomodelxml.substmodel.HKYParser;
import dr.oldevomodelxml.treelikelihood.TreeLikelihoodParser;
import dr.inference.loggers.ArrayLogFormatter;
import dr.inference.loggers.MCLogger;
import dr.inference.loggers.TabDelimitedFormatter;
import dr.inference.mcmc.MCMC;
import dr.inference.mcmc.MCMCOptions;
import dr.inference.model.Parameter;
import dr.inference.operators.*;
import dr.inference.trace.ArrayTraceList;
import dr.inference.trace.Trace;
import dr.inference.trace.TraceCorrelation;
import dr.math.MathUtils;
import junit.framework.Test;
import junit.framework.TestSuite;
import test.dr.inference.trace.TraceCorrelationAssert;

import java.util.List;

/**
 * @author Walter Xie
 * convert testMCMC.xml in the folder /example
 */

public class MCMCTest extends TraceCorrelationAssert {

    public MCMCTest(String name) {
        super(name);
    }

    public void setUp() throws Exception {
        super.setUp();

        MathUtils.setSeed(666);

        createAlignment(PRIMATES_TAXON_SEQUENCE, Nucleotides.INSTANCE);

        createRandomInitialTree(0.0001); // popSize

//        createSpecifiedTree("((((human:0.02124198428146588,(bonobo:0.010505698073024256,chimp:0.010505698073024256)" +
//                ":0.010736286208441624):0.011019735965429791,gorilla:0.03226172024689567):0.022501552046463147," +
//                "orangutan:0.05476327229335882):0.009440823865408586,siamang:0.0642040961587674);");
    }


    public void testMCMC() {
        Parameter kappa = new Parameter.Default(HKYParser.KAPPA, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);

        dr.evomodel.substmodel.FrequencyModel f = new FrequencyModel(Nucleotides.INSTANCE, new Parameter.Default(alignment.getStateFrequencies()));
        dr.evomodel.substmodel.nucleotide.HKY hky = new HKY(kappa, f);

        //siteModel
        GammaSiteRateModel siteRateModel = new GammaSiteRateModel("homogeneious");

        //treeLikelihood
        TreeDataLikelihood treeDataLikelihood = getTreeDataLikelihood(hky, siteRateModel);

        treeDataLikelihood.setId(TreeDataLikelihoodParser.TREE_DATA_LIKELIHOOD);

        // Operators
        OperatorSchedule schedule = new SimpleOperatorSchedule();

        MCMCOperator operator = new ScaleOperator(kappa, 0.5);
        operator.setWeight(1.0);
        schedule.addOperator(operator);

        Parameter rootHeight = ((DefaultTreeModel)treeModel).getRootHeightParameter();
        rootHeight.setId(TREE_HEIGHT);
        operator = new ScaleOperator(rootHeight, 0.5);
        operator.setWeight(1.0);
        schedule.addOperator(operator);

        Parameter internalHeights = ((DefaultTreeModel)treeModel).createNodeHeightsParameter(false, true, false);
        operator = new UniformOperator(internalHeights, 10.0);
        schedule.addOperator(operator);

        operator = new SubtreeSlideOperator(((DefaultTreeModel)treeModel), 1, 1, true, false, false, false, AdaptationMode.ADAPTATION_ON, AdaptableMCMCOperator.DEFAULT_ADAPTATION_TARGET);
        schedule.addOperator(operator);

        operator = new ExchangeOperator(ExchangeOperator.NARROW, treeModel, 1.0);
//        operator.doOperation();
        schedule.addOperator(operator);

        operator = new ExchangeOperator(ExchangeOperator.WIDE, treeModel, 1.0);
//        operator.doOperation();
        schedule.addOperator(operator);

        operator = new WilsonBalding(treeModel, 1.0);
//        operator.doOperation();
        schedule.addOperator(operator);

        // Log
        ArrayLogFormatter formatter = new ArrayLogFormatter(false);

        MCLogger[] loggers = new MCLogger[2];
        loggers[0] = new MCLogger(formatter, 1000, false);
        loggers[0].add(treeDataLikelihood);
        loggers[0].add(rootHeight);
        loggers[0].add(kappa);

        loggers[1] = new MCLogger(new TabDelimitedFormatter(System.out), 100000, false);
        loggers[1].add(treeDataLikelihood);
        loggers[1].add(rootHeight);
        loggers[1].add(kappa);

        // MCMC
        MCMC mcmc = new MCMC("mcmc1");
        MCMCOptions options = new MCMCOptions(10000000);

        mcmc.setShowOperatorAnalysis(true);
        mcmc.init(options, treeDataLikelihood, schedule, loggers);
        mcmc.run();

        // time
        System.out.println(mcmc.getTimer().toString());

        // Tracer
        List<Trace> traces = formatter.getTraces();
        ArrayTraceList traceList = new ArrayTraceList("MCMCTest", traces, 0);

        for (int i = 1; i < traces.size(); i++) {
            traceList.analyseTrace(i);
        }

//      <expectation name="likelihood" value="-1815.75"/>
//		<expectation name="treeModel.rootHeight" value="6.42048E-2"/>
//		<expectation name="hky.kappa" value="32.8941"/>

        TraceCorrelation likelihoodStats = traceList.getCorrelationStatistics(traceList.getTraceIndex(TreeDataLikelihoodParser.TREE_DATA_LIKELIHOOD));
        assertExpectation(TreeLikelihoodParser.TREE_LIKELIHOOD, likelihoodStats, -1815.75);

        TraceCorrelation treeHeightStats = traceList.getCorrelationStatistics(traceList.getTraceIndex(TREE_HEIGHT));
        assertExpectation(TREE_HEIGHT, treeHeightStats, 6.42048E-2);

        TraceCorrelation kappaStats = traceList.getCorrelationStatistics(traceList.getTraceIndex(HKYParser.KAPPA));
        assertExpectation(HKYParser.KAPPA, kappaStats, 32.8941);
    }

    private TreeDataLikelihood getTreeDataLikelihood(SubstitutionModel substitutionModel, GammaSiteRateModel siteRateModel) {
        SitePatterns patterns = new SitePatterns(alignment, null, 0, -1, 1, true);

        DataLikelihoodDelegate dataLikelihoodDelegate = new BeagleDataLikelihoodDelegate(
                treeModel,
                patterns,
                new HomogeneousBranchModel(substitutionModel),
                siteRateModel,
                false,
                false,
                PartialsRescalingScheme.DEFAULT,
                false,
                new PreOrderSettings(false, false, false, true)
        );

        TreeDataLikelihood treeDataLikelihood = new TreeDataLikelihood(
                dataLikelihoodDelegate,
                treeModel,
                new DefaultBranchRateModel());
        return treeDataLikelihood;
    }
    public static Test suite() {
        return new TestSuite(MCMCTest.class);
    }
}