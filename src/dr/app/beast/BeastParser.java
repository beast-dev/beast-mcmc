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

import dr.xml.*;

import java.util.Iterator;

/**
 * @author Alexei Drummond
 * @author Andrew Rambaut
 *
 * @version $Id: BeastParser.java,v 1.76 2006/08/30 16:01:59 rambaut Exp $
 */
public class BeastParser extends XMLParser {

	public BeastParser(String[] args)  throws XMLParseException {
		super();
		setup(args);
	}

	public BeastParser(String[] args, boolean verbose) throws XMLParseException {
		super(verbose);
		setup(args);

		if (verbose) {
			Iterator iterator = getParsers();
			while (iterator.hasNext()) {
				XMLObjectParser parser = (XMLObjectParser)iterator.next();
				System.out.println(parser.toString());
			}
		}
	}

	private final void setup(String[] args) throws XMLParseException {

		for (int i = 0; i < args.length; i++) {
			storeObject(Integer.toString(i), args[i]);
		}

			// add all the XMLObject parsers you need

		addXMLObjectParser(new PropertyParser());
		addXMLObjectParser(UserInput.STRING_PARSER);
		addXMLObjectParser(UserInput.DOUBLE_PARSER);
		addXMLObjectParser(UserInput.INTEGER_PARSER);

		addXMLObjectParser(new dr.evoxml.GeneralDataTypeParser());
		addXMLObjectParser(new dr.evoxml.AlignmentParser());
		addXMLObjectParser(new dr.evoxml.SitePatternsParser());
		addXMLObjectParser(new dr.evoxml.ConvertAlignmentParser());
		addXMLObjectParser(new dr.evoxml.MergePatternsParser());
		addXMLObjectParser(new dr.evoxml.SequenceParser());
		addXMLObjectParser(new dr.evoxml.SimpleNodeParser());
		addXMLObjectParser(new dr.evoxml.SimpleTreeParser());
		addXMLObjectParser(new dr.evoxml.UPGMATreeParser());
		addXMLObjectParser(new dr.evoxml.NeighborJoiningParser());
		addXMLObjectParser(new dr.evoxml.NewickParser());
		addXMLObjectParser(new dr.evoxml.TaxonParser());
		addXMLObjectParser(new dr.evoxml.TaxaParser());
		addXMLObjectParser(new dr.evoxml.DateParser());
        addXMLObjectParser(new dr.evoxml.DistanceMatrixParser());

		addXMLObjectParser(new dr.xml.AttributeParser());
		addXMLObjectParser(new dr.xml.AttributesParser());


        addXMLObjectParser(dr.evomodel.speciation.SpeciationLikelihood.PARSER);
        addXMLObjectParser(dr.evomodel.speciation.YuleModel.PARSER);
        addXMLObjectParser(dr.evomodel.speciation.BirthDeathModel.PARSER);
        addXMLObjectParser(dr.evomodel.speciation.BranchingLikelihood.PARSER);
        addXMLObjectParser(dr.evomodel.speciation.BetaSplittingModel.PARSER);

		addXMLObjectParser(dr.evomodel.coalescent.CoalescentSimulator.PARSER);
		addXMLObjectParser(dr.evomodel.coalescent.CoalescentLikelihood.PARSER);
		addXMLObjectParser(dr.evomodel.coalescent.CoalescentMRCALikelihood.PARSER);
		addXMLObjectParser(dr.evomodel.coalescent.SkylineLikelihood.PARSER);
		addXMLObjectParser(dr.evomodel.coalescent.BayesianSkylineLikelihood.PARSER);
		addXMLObjectParser(dr.evomodel.coalescent.ConstantPopulationModel.PARSER);
		addXMLObjectParser(dr.evomodel.coalescent.ExponentialGrowthModel.PARSER);
		addXMLObjectParser(dr.evomodel.coalescent.CataclysmicDemographicModel.PARSER);
        addXMLObjectParser(dr.evomodel.coalescent.ExpConstExpDemographicModel.PARSER);
		addXMLObjectParser(dr.evomodel.coalescent.LogisticGrowthModel.PARSER);
		addXMLObjectParser(dr.evomodel.coalescent.ConstantExponentialModel.PARSER);
		addXMLObjectParser(dr.evomodel.coalescent.ExpansionModel.PARSER);
		addXMLObjectParser(dr.evomodel.coalescent.ConstantLogisticModel.PARSER);
		addXMLObjectParser(dr.evomodel.coalescent.PiecewisePopulationModel.PARSER);
		addXMLObjectParser(dr.evomodel.coalescent.ScaledPiecewiseModel.PARSER);
		addXMLObjectParser(dr.evomodel.coalescent.TwoEpochDemographicModel.PARSER);
		addXMLObjectParser(dr.evomodel.coalescent.ExponentialSawtoothModel.PARSER);
		addXMLObjectParser(dr.evomodel.coalescent.EmpiricalPiecewiseModel.PARSER);
//        addXMLObjectParser(dr.evomodel.coalescent.PopSizeStatistic.PARSER);

		// Structured coalescent
		addXMLObjectParser(dr.evomodel.coalescent.structure.StructuredCoalescentLikelihood.PARSER);
		addXMLObjectParser(dr.evomodel.coalescent.structure.ConstantMigrationModel.PARSER);
		addXMLObjectParser(dr.evomodel.coalescent.structure.ColourSamplerModel.PARSER);
		addXMLObjectParser(dr.evomodel.coalescent.structure.TreeColouringOperator.PARSER);
		addXMLObjectParser(dr.evomodel.coalescent.structure.MetaPopulationModel.PARSER);
		addXMLObjectParser(dr.evomodel.operators.ColouredExchangeOperator.NARROW_EXCHANGE_PARSER);
		addXMLObjectParser(dr.evomodel.operators.ColouredExchangeOperator.WIDE_EXCHANGE_PARSER);
		addXMLObjectParser(dr.evomodel.operators.ColouredSubtreeSlideOperator.PARSER);
		addXMLObjectParser(dr.evomodel.operators.ColouredOperator.PARSER);
		addXMLObjectParser(dr.evomodel.operators.FixedColouredOperator.PARSER);

        addXMLObjectParser(dr.evomodel.substmodel.FrequencyModel.PARSER);
		addXMLObjectParser(dr.evomodel.substmodel.GeneralSubstitutionModel.PARSER);
		addXMLObjectParser(dr.evomodel.substmodel.HKY.PARSER);
		addXMLObjectParser(dr.evomodel.substmodel.GTR.PARSER);
		addXMLObjectParser(dr.evomodel.substmodel.EmpiricalAminoAcidModel.PARSER);
		addXMLObjectParser(dr.evomodel.substmodel.YangCodonModel.PARSER);

		addXMLObjectParser(dr.evomodel.treelikelihood.TreeLikelihood.PARSER);
		addXMLObjectParser(dr.evomodel.treelikelihood.AdvancedTreeLikelihood.PARSER);

		addXMLObjectParser(dr.evomodel.sitemodel.GammaSiteModel.PARSER);

		addXMLObjectParser(dr.evomodel.clock.EDLikelihood.PARSER);
		addXMLObjectParser(dr.evomodel.clock.NDLikelihood.PARSER);
		addXMLObjectParser(dr.evomodel.clock.GDLikelihood.PARSER);
        addXMLObjectParser(dr.evomodel.branchratemodel.DiscretizedBranchRates.PARSER);
        addXMLObjectParser(dr.evomodel.branchratemodel.StrictClockBranchRates.PARSER);

		addXMLObjectParser(dr.evomodel.tree.TreeModel.PARSER);
        addXMLObjectParser(dr.evomodel.tree.TipHeightLikelihood.PARSER);
		addXMLObjectParser(dr.evomodel.tree.TreeMetricStatistic.PARSER);
		addXMLObjectParser(dr.evomodel.tree.TreeShapeStatistic.PARSER);
		addXMLObjectParser(dr.evomodel.tree.TMRCAStatistic.PARSER);
        addXMLObjectParser(dr.evomodel.tree.RateCovarianceStatistic.PARSER);
        addXMLObjectParser(dr.evomodel.tree.RateStatistic.PARSER);
		addXMLObjectParser(dr.evomodel.tree.MonophylyStatistic.PARSER);
		addXMLObjectParser(dr.evomodel.tree.ParsimonyStatistic.PARSER);
		addXMLObjectParser(dr.evomodel.tree.ParsimonyStateStatistic.PARSER);
		addXMLObjectParser(dr.evomodel.tree.SpeciesTreeStatistic.PARSER);
        addXMLObjectParser(dr.evomodel.tree.UniformRootPrior.PARSER);

		addXMLObjectParser(dr.evomodel.operators.ExchangeOperator.NARROW_EXCHANGE_PARSER);
		addXMLObjectParser(dr.evomodel.operators.ExchangeOperator.WIDE_EXCHANGE_PARSER);
		addXMLObjectParser(dr.evomodel.operators.WilsonBalding.PARSER);
		addXMLObjectParser(dr.evomodel.operators.SubtreeSlideOperator.PARSER);
		addXMLObjectParser(dr.evomodel.operators.RateExchangeOperator.PARSER);

		addXMLObjectParser(dr.evomodel.indel.TKF91Likelihood.PARSER);
		addXMLObjectParser(dr.evomodel.indel.TKF91Model.PARSER);
		addXMLObjectParser(dr.evomodel.indel.IstvanOperator.PARSER);

		addXMLObjectParser(dr.inference.model.CompoundParameter.PARSER);
		addXMLObjectParser(dr.inference.model.CompoundLikelihood.PARSER);
		addXMLObjectParser(dr.inference.model.BooleanLikelihood.PARSER);
        addXMLObjectParser(dr.inference.model.DummyLikelihood.PARSER);
		addXMLObjectParser(dr.inference.model.JeffreysPriorLikelihood.PARSER);

		addXMLObjectParser(dr.inference.distribution.DistributionLikelihood.PARSER);
		addXMLObjectParser(dr.inference.distribution.EmpiricalDistributionLikelihood.PARSER);
        addXMLObjectParser(dr.inference.distribution.UniformDistributionModel.PARSER);
		addXMLObjectParser(dr.inference.distribution.ExponentialDistributionModel.PARSER);
		addXMLObjectParser(dr.inference.distribution.GammaDistributionModel.PARSER);
		addXMLObjectParser(dr.inference.distribution.NormalDistributionModel.PARSER);
        addXMLObjectParser(dr.inference.distribution.LogNormalDistributionModel.PARSER);
		addXMLObjectParser(dr.inference.distribution.ExponentialMarkovModel.PARSER);

        addXMLObjectParser(dr.inference.distribution.DistributionLikelihood.UNIFORM_PRIOR_PARSER);
        addXMLObjectParser(dr.inference.distribution.DistributionLikelihood.EXPONENTIAL_PRIOR_PARSER);
        addXMLObjectParser(dr.inference.distribution.DistributionLikelihood.NORMAL_PRIOR_PARSER);
        addXMLObjectParser(dr.inference.distribution.DistributionLikelihood.LOG_NORMAL_PRIOR_PARSER);
        addXMLObjectParser(dr.inference.distribution.DistributionLikelihood.GAMMA_PRIOR_PARSER);

		addXMLObjectParser(dr.inference.distribution.BinomialLikelihood.PARSER);

		addXMLObjectParser(new dr.inference.model.StatisticParser());
		addXMLObjectParser(new dr.inference.model.ParameterParser());
		addXMLObjectParser(dr.inference.model.TestStatistic.PARSER);
        addXMLObjectParser(dr.inference.model.VarianceStatistic.PARSER);

		addXMLObjectParser(dr.inference.model.MeanStatistic.PARSER);
        addXMLObjectParser(dr.inference.model.ProductStatistic.PARSER);
        addXMLObjectParser(dr.inference.model.SumStatistic.PARSER);
        addXMLObjectParser(dr.inference.model.ReciprocalStatistic.PARSER);
        addXMLObjectParser(dr.inference.model.NegateStatistic.PARSER);
        addXMLObjectParser(dr.inference.model.ExponentiationStatistic.PARSER);
        addXMLObjectParser(dr.inference.model.ExpressionStatistic.PARSER);

        addXMLObjectParser(dr.inference.ml.MLOptimizer.PARSER);
        addXMLObjectParser(dr.inference.mcmc.MCMC.PARSER);
		addXMLObjectParser(dr.inference.loggers.MCLogger.PARSER);
		addXMLObjectParser(dr.inference.loggers.MLLogger.ML_LOGGER_PARSER);
		addXMLObjectParser(dr.evomodel.tree.TreeLogger.PARSER);
		addXMLObjectParser(dr.inference.loggers.Columns.PARSER);
		addXMLObjectParser(dr.inference.operators.SimpleOperatorSchedule.PARSER);

		addXMLObjectParser(dr.inference.operators.RandomWalkOperator.PARSER);
		addXMLObjectParser(dr.inference.operators.ScaleOperator.PARSER);
		addXMLObjectParser(dr.inference.operators.UniformOperator.PARSER);
		addXMLObjectParser(dr.inference.operators.UpDownOperator.PARSER);
		addXMLObjectParser(dr.inference.operators.SetOperator.PARSER);
        addXMLObjectParser(dr.inference.operators.SwapOperator.PARSER);
		addXMLObjectParser(dr.inference.operators.DeltaExchangeOperator.PARSER);
		addXMLObjectParser(dr.inference.operators.CenteredScaleOperator.PARSER);

		addXMLObjectParser(new dr.evomodel.tree.TreeTraceAnalysisParser());
		addXMLObjectParser(new dr.inference.trace.TraceAnalysisParser());
	}

}

