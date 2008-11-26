package dr.app.beauti.generator;

import dr.app.beauti.XMLWriter;
import dr.app.beauti.options.BeautiOptions;
import dr.app.beauti.options.ModelOptions;
import dr.app.beauti.options.TreePrior;
import dr.app.beauti.options.StartingTreeType;
import dr.evolution.util.Units;
import dr.evomodel.coalescent.*;
import dr.evomodel.speciation.BirthDeathGernhard08Model;
import dr.evomodel.speciation.SpeciationLikelihood;
import dr.evomodel.speciation.YuleModel;
import dr.evomodel.tree.TreeModel;
import dr.evomodelxml.BirthDeathModelParser;
import dr.evomodelxml.YuleModelParser;
import dr.exporters.CSVExporter;
import dr.inference.distribution.ExponentialDistributionModel;
import dr.inference.distribution.ExponentialMarkovModel;
import dr.inference.model.BooleanLikelihood;
import dr.inference.model.ParameterParser;
import dr.inference.model.SumStatistic;
import dr.inference.model.TestStatistic;
import dr.inference.trace.EBSPAnalysis;
import dr.util.Attribute;

/**
 * @author Alexei Drummond
 */
public class TreePriorGenerator extends Generator {

	public TreePriorGenerator(BeautiOptions options) {
		super(options);
	}

	void writeTreePrior(XMLWriter writer) {

		writeNodeHeightPrior(writer);
		if (options.nodeHeightPrior == TreePrior.LOGISTIC) {
			writer.writeText("");
			writeBooleanLikelihood(writer);
		} else if (options.nodeHeightPrior == TreePrior.SKYLINE) {
			writer.writeText("");
			writeExponentialMarkovLikelihood(writer);
		}
	}

	/**
	 * Write a tree prior (coalescent or speciational) model
	 *
	 * @param writer the writer
	 */
	void writeTreePriorModel(XMLWriter writer) {

		String initialPopSize = null;

		TreePrior nodeHeightPrior = options.nodeHeightPrior;
		Units.Type units = options.units;
		int parameterization = options.parameterization;

		if (nodeHeightPrior == TreePrior.CONSTANT) {

			writer.writeComment("A prior assumption that the population size has remained constant");
			writer.writeComment("throughout the time spanned by the genealogy.");
			writer.writeOpenTag(
					ConstantPopulationModel.CONSTANT_POPULATION_MODEL,
					new Attribute[]{
							new Attribute.Default<String>("id", "constant"),
							new Attribute.Default<String>("units", Units.Utils.getDefaultUnitName(options.units))
					}
			);

			writer.writeOpenTag(ConstantPopulationModel.POPULATION_SIZE);
			writeParameter("constant.popSize", options, writer);
			writer.writeCloseTag(ConstantPopulationModel.POPULATION_SIZE);
			writer.writeCloseTag(ConstantPopulationModel.CONSTANT_POPULATION_MODEL);

		} else if (options.nodeHeightPrior == TreePrior.EXPONENTIAL) {
			// generate an exponential prior tree

			writer.writeComment("A prior assumption that the population size has grown exponentially");
			writer.writeComment("throughout the time spanned by the genealogy.");
			writer.writeOpenTag(
					ExponentialGrowthModel.EXPONENTIAL_GROWTH_MODEL,
					new Attribute[]{
							new Attribute.Default<String>("id", "exponential"),
							new Attribute.Default<String>("units", Units.Utils.getDefaultUnitName(options.units))
					}
			);

			// write pop size socket
			writer.writeOpenTag(ExponentialGrowthModel.POPULATION_SIZE);
			writeParameter("exponential.popSize", options, writer);
			writer.writeCloseTag(ExponentialGrowthModel.POPULATION_SIZE);

			if (parameterization == ModelOptions.GROWTH_RATE) {
				// write growth rate socket
				writer.writeOpenTag(ExponentialGrowthModel.GROWTH_RATE);
				writeParameter("exponential.growthRate", options, writer);
				writer.writeCloseTag(ExponentialGrowthModel.GROWTH_RATE);
			} else {
				// write doubling time socket
				writer.writeOpenTag(ExponentialGrowthModel.DOUBLING_TIME);
				writeParameter("exponential.doublingTime", options, writer);
				writer.writeCloseTag(ExponentialGrowthModel.DOUBLING_TIME);
			}

			writer.writeCloseTag(ExponentialGrowthModel.EXPONENTIAL_GROWTH_MODEL);
		} else if (options.nodeHeightPrior == TreePrior.LOGISTIC) {
			// generate an exponential prior tree

			writer.writeComment("A prior assumption that the population size has grown logistically");
			writer.writeComment("throughout the time spanned by the genealogy.");
			writer.writeOpenTag(
					LogisticGrowthModel.LOGISTIC_GROWTH_MODEL,
					new Attribute[]{
							new Attribute.Default<String>("id", "logistic"),
							new Attribute.Default<String>("units", Units.Utils.getDefaultUnitName(options.units))
					}
			);

			// write pop size socket
			writer.writeOpenTag(LogisticGrowthModel.POPULATION_SIZE);
			writeParameter("logistic.popSize", options, writer);
			writer.writeCloseTag(LogisticGrowthModel.POPULATION_SIZE);

			if (parameterization == ModelOptions.GROWTH_RATE) {
				// write growth rate socket
				writer.writeOpenTag(LogisticGrowthModel.GROWTH_RATE);
				writeParameter("logistic.growthRate", options, writer);
				writer.writeCloseTag(LogisticGrowthModel.GROWTH_RATE);
			} else {
				// write doubling time socket
				writer.writeOpenTag(LogisticGrowthModel.DOUBLING_TIME);
				writeParameter("logistic.doublingTime", options, writer);
				writer.writeCloseTag(LogisticGrowthModel.DOUBLING_TIME);
			}

			// write logistic t50 socket
			writer.writeOpenTag(LogisticGrowthModel.TIME_50);
			writeParameter("logistic.t50", options, writer);
			writer.writeCloseTag(LogisticGrowthModel.TIME_50);

			writer.writeCloseTag(LogisticGrowthModel.LOGISTIC_GROWTH_MODEL);

			initialPopSize = "logistic.popSize";

		} else if (options.nodeHeightPrior == TreePrior.EXPANSION) {
			// generate an exponential prior tree

			writer.writeComment("A prior assumption that the population size has grown exponentially");
			writer.writeComment("from some ancestral population size in the past.");
			writer.writeOpenTag(
					ExpansionModel.EXPANSION_MODEL,
					new Attribute[]{
							new Attribute.Default<String>("id", "expansion"),
							new Attribute.Default<String>("units", Units.Utils.getDefaultUnitName(options.units))
					}
			);

			// write pop size socket
			writeParameter(ExpansionModel.POPULATION_SIZE, "expansion.popSize", options, writer);

			if (options.parameterization == ModelOptions.GROWTH_RATE) {
				// write growth rate socket
				writeParameter(ExpansionModel.GROWTH_RATE, "expansion.growthRate", options, writer);
			} else {
				// write doubling time socket
				writeParameter(ExpansionModel.DOUBLING_TIME, "expansion.doublingTime", options, writer);
			}

			// write ancestral proportion socket
			writeParameter(ExpansionModel.ANCESTRAL_POPULATION_PROPORTION,
					"expansion.ancestralProportion", options, writer);

			writer.writeCloseTag(ExpansionModel.EXPANSION_MODEL);

			initialPopSize = "expansion.popSize";

		} else if (options.nodeHeightPrior == TreePrior.YULE) {
			writer.writeComment("A prior on the distribution node heights defined given");
			writer.writeComment("a Yule speciation process (a pure birth process).");
			writer.writeOpenTag(
					YuleModel.YULE_MODEL,
					new Attribute[]{
							new Attribute.Default<String>("id", "yule"),
							new Attribute.Default<String>("units", Units.Utils.getDefaultUnitName(units))
					}
			);

			writeParameter(YuleModelParser.BIRTH_RATE, "yule.birthRate", options, writer);
			writer.writeCloseTag(YuleModel.YULE_MODEL);
		} else if (nodeHeightPrior == TreePrior.BIRTH_DEATH) {
			writer.writeComment("A prior on the distribution node heights defined given");
			writer.writeComment("a Birth-Death speciation process (Gernhard 2008).");
			writer.writeOpenTag(
					BirthDeathGernhard08Model.BIRTH_DEATH_MODEL,
					new Attribute[]{
							new Attribute.Default<String>("id", "birthDeath"),
							new Attribute.Default<String>("units", Units.Utils.getDefaultUnitName(units))
					}
			);

			writeParameter(BirthDeathModelParser.BIRTHDIFF_RATE,
					BirthDeathModelParser.BIRTHDIFF_RATE_PARAM_NAME, options, writer);
			writeParameter(BirthDeathModelParser.RELATIVE_DEATH_RATE,
					BirthDeathModelParser.RELATIVE_DEATH_RATE_PARAM_NAME, options, writer);
			writer.writeCloseTag(BirthDeathGernhard08Model.BIRTH_DEATH_MODEL);
		}

		if (nodeHeightPrior != TreePrior.CONSTANT && nodeHeightPrior != TreePrior.EXPONENTIAL) {
			// If the node height prior is not one of these two then we need to simulate a
			// random starting tree under a constant size coalescent.

			writer.writeComment("This is a simple constant population size coalescent model");
			writer.writeComment("that is used to generate an initial tree for the chain.");
			writer.writeOpenTag(
					ConstantPopulationModel.CONSTANT_POPULATION_MODEL,
					new Attribute[]{
							new Attribute.Default<String>("id", "initialDemo"),
							new Attribute.Default<String>("units", Units.Utils.getDefaultUnitName(units))
					}
			);

			writer.writeOpenTag(ConstantPopulationModel.POPULATION_SIZE);
			if (initialPopSize != null) {
				writer.writeTag(ParameterParser.PARAMETER,
						new Attribute[]{
								new Attribute.Default<String>("idref", initialPopSize),
						}, true);
			} else {
				writeParameter("initialDemo.popSize", 1, 100.0, Double.NaN, Double.NaN, writer);
			}
			writer.writeCloseTag(ConstantPopulationModel.POPULATION_SIZE);
			writer.writeCloseTag(ConstantPopulationModel.CONSTANT_POPULATION_MODEL);
		}
	}

	/**
	 * Write the prior on node heights (coalescent or speciational models)
	 *
	 * @param writer the writer
	 */
	private void writeNodeHeightPrior(XMLWriter writer) {

		TreePrior treePrior = options.nodeHeightPrior;

		if (treePrior == TreePrior.YULE || treePrior == TreePrior.BIRTH_DEATH) {
			// generate a speciational process

			writer.writeOpenTag(
					SpeciationLikelihood.SPECIATION_LIKELIHOOD,
					new Attribute[]{
							new Attribute.Default<String>("id", "speciation")
					}
			);

			// write pop size socket
			writer.writeOpenTag(SpeciationLikelihood.MODEL);
			writeNodeHeightPriorModelRef(writer);
			writer.writeCloseTag(SpeciationLikelihood.MODEL);
			writer.writeOpenTag(SpeciationLikelihood.TREE);
			writer.writeTag(TreeModel.TREE_MODEL, new Attribute.Default<String>("idref", "treeModel"), true);
			writer.writeCloseTag(SpeciationLikelihood.TREE);

			writer.writeCloseTag(SpeciationLikelihood.SPECIATION_LIKELIHOOD);

		} else if (treePrior == TreePrior.SKYLINE) {
			// generate a Bayesian skyline plot

			writer.writeOpenTag(
					BayesianSkylineLikelihood.SKYLINE_LIKELIHOOD,
					new Attribute[]{
							new Attribute.Default<String>("id", "skyline"),
							new Attribute.Default<String>("linear",
									options.skylineModel == ModelOptions.LINEAR_SKYLINE ? "true" : "false")
					}
			);

			// write pop size socket
			writer.writeOpenTag(BayesianSkylineLikelihood.POPULATION_SIZES);
			if (options.skylineModel == ModelOptions.LINEAR_SKYLINE) {
				writeParameter("skyline.popSize", options.skylineGroupCount + 1, writer);
			} else {
				writeParameter("skyline.popSize", options.skylineGroupCount, writer);
			}
			writer.writeCloseTag(BayesianSkylineLikelihood.POPULATION_SIZES);

			// write group size socket
			writer.writeOpenTag(BayesianSkylineLikelihood.GROUP_SIZES);
			writeParameter("skyline.groupSize", options.skylineGroupCount, writer);
			writer.writeCloseTag(BayesianSkylineLikelihood.GROUP_SIZES);

			writer.writeOpenTag(CoalescentLikelihood.POPULATION_TREE);
			writer.writeTag(TreeModel.TREE_MODEL, new Attribute.Default<String>("idref", "treeModel"), true);
			writer.writeCloseTag(CoalescentLikelihood.POPULATION_TREE);

			writer.writeCloseTag(BayesianSkylineLikelihood.SKYLINE_LIKELIHOOD);
		} else if (treePrior == TreePrior.EXTENDED_SKYLINE) {
			final String tagName = VariableDemographicModel.PARSER.getParserName();

			writer.writeOpenTag(
					tagName,
					new Attribute[]{
							new Attribute.Default<String>("id", VariableDemographicModel.demoElementName),
							new Attribute.Default<String>(VariableDemographicModel.TYPE, options.extendedSkylineModel),
                            // use midpoint by default (todo) would be nice to have a user 'tickable' option
                            new Attribute.Default<String>(VariableDemographicModel.USE_MIDPOINTS, "true")
                    }
			);

			writer.writeOpenTag(VariableDemographicModel.POPULATION_SIZES);
			final int nTax = options.taxonList.getTaxonCount();
			final int nPops = nTax - (options.extendedSkylineModel.equals(VariableDemographicModel.STEPWISE) ? 1 : 0);
			writeParameter(VariableDemographicModel.demoElementName + ".popSize", nPops, writer);
			writer.writeCloseTag(VariableDemographicModel.POPULATION_SIZES);

			writer.writeOpenTag(VariableDemographicModel.INDICATOR_PARAMETER);
			writeParameter(VariableDemographicModel.demoElementName + ".indicators", nPops - 1, writer);
			writer.writeCloseTag(VariableDemographicModel.INDICATOR_PARAMETER);

			writer.writeOpenTag(VariableDemographicModel.POPULATION_TREES);

			writer.writeOpenTag(VariableDemographicModel.POP_TREE);
			writer.writeTag(TreeModel.TREE_MODEL, new Attribute.Default<String>("idref", "treeModel"), true);
			writer.writeCloseTag(VariableDemographicModel.POP_TREE);

			writer.writeCloseTag(VariableDemographicModel.POPULATION_TREES);

			writer.writeCloseTag(tagName);

			writer.writeOpenTag(CoalescentLikelihood.COALESCENT_LIKELIHOOD, new Attribute.Default<String>("id", "coalescent"));
			writer.writeOpenTag(CoalescentLikelihood.MODEL);
			writer.writeTag(tagName, new Attribute.Default<String>("idref", VariableDemographicModel.demoElementName), true);
			writer.writeCloseTag(CoalescentLikelihood.MODEL);
			writer.writeComment("Take population Tree from demographic");
			writer.writeCloseTag(CoalescentLikelihood.COALESCENT_LIKELIHOOD);

			writer.writeOpenTag(SumStatistic.SUM_STATISTIC,
					new Attribute[]{
							new Attribute.Default<String>("id", VariableDemographicModel.demoElementName + ".populationSizeChanges"),
							new Attribute.Default<String>("elementwise", "true")
					});
			writer.writeTag(ParameterParser.PARAMETER,
					new Attribute.Default<String>("idref", VariableDemographicModel.demoElementName + ".indicators"), true);
			writer.writeCloseTag(SumStatistic.SUM_STATISTIC);
			writer.writeOpenTag(ExponentialDistributionModel.EXPONENTIAL_DISTRIBUTION_MODEL,
					new Attribute[]{
							new Attribute.Default<String>("id", VariableDemographicModel.demoElementName + ".populationMeanDist")
							//,new Attribute.Default<String>("elementwise", "true")
					});
			writer.writeOpenTag(ExponentialDistributionModel.MEAN);
			writer.writeTag(ParameterParser.PARAMETER,
					new Attribute[]{
							new Attribute.Default<String>("id", VariableDemographicModel.demoElementName + ".populationMean"),
							new Attribute.Default<String>("value", "1")}, true);
			writer.writeCloseTag(ExponentialDistributionModel.MEAN);
			writer.writeCloseTag(ExponentialDistributionModel.EXPONENTIAL_DISTRIBUTION_MODEL);
		} else if (treePrior == TreePrior.GMRF_SKYRIDE) {
			writer.writeOpenTag(
					GMRFSkyrideLikelihood.SKYLINE_LIKELIHOOD,
					new Attribute[]{
							new Attribute.Default<String>("id", "skyride"),
							new Attribute.Default<String>(GMRFSkyrideLikelihood.TIME_AWARE_SMOOTHING,
									options.skyrideSmoothing == ModelOptions.SKYRIDE_TIME_AWARE_SMOOTHING ? "true" : "false"),
                            new Attribute.Default<String>(GMRFSkyrideLikelihood.RANDOMIZE_TREE,
                                    options.startingTreeType == StartingTreeType.UPGMA ? "true" : "false"),
					}
			);

			writer.writeOpenTag(GMRFSkyrideLikelihood.POPULATION_PARAMETER);
			writeParameter("skyride.popSize", options.skyrideIntervalCount, writer);
			writer.writeCloseTag(GMRFSkyrideLikelihood.POPULATION_PARAMETER);

			writer.writeOpenTag(GMRFSkyrideLikelihood.GROUP_SIZES);
			writeParameter("skyride.groupSize", options.skyrideIntervalCount, writer);
			writer.writeCloseTag(GMRFSkyrideLikelihood.GROUP_SIZES);

			writer.writeOpenTag(GMRFSkyrideLikelihood.PRECISION_PARAMETER);
			writeParameter("skyride.precision", 1, writer);
			writer.writeCloseTag(GMRFSkyrideLikelihood.PRECISION_PARAMETER);

			writer.writeOpenTag(GMRFSkyrideLikelihood.POPULATION_TREE);
			writer.writeTag(TreeModel.TREE_MODEL, new Attribute.Default<String>("idref", "treeModel"), true);
			writer.writeCloseTag(GMRFSkyrideLikelihood.POPULATION_TREE);

			writer.writeCloseTag(GMRFSkyrideLikelihood.SKYLINE_LIKELIHOOD);

		} else {
			// generate a coalescent process

			writer.writeOpenTag(
					CoalescentLikelihood.COALESCENT_LIKELIHOOD,
					new Attribute[]{new Attribute.Default<String>("id", "coalescent")}
			);
			writer.writeOpenTag(CoalescentLikelihood.MODEL);
			writeNodeHeightPriorModelRef(writer);
			writer.writeCloseTag(CoalescentLikelihood.MODEL);
			writer.writeOpenTag(CoalescentLikelihood.POPULATION_TREE);
			writer.writeTag(TreeModel.TREE_MODEL, new Attribute.Default<String>("idref", "treeModel"), true);
			writer.writeCloseTag(CoalescentLikelihood.POPULATION_TREE);
			writer.writeCloseTag(CoalescentLikelihood.COALESCENT_LIKELIHOOD);
		}
	}

	void writeNodeHeightPriorModelRef(XMLWriter writer) {

		TreePrior treePrior = options.nodeHeightPrior;

		switch (treePrior) {
			case CONSTANT:
				writer.writeTag(ConstantPopulationModel.CONSTANT_POPULATION_MODEL,
						new Attribute[]{new Attribute.Default<String>("idref", "constant")}, true);
				break;
			case EXPONENTIAL:
				writer.writeTag(ExponentialGrowthModel.EXPONENTIAL_GROWTH_MODEL,
						new Attribute[]{new Attribute.Default<String>("idref", "exponential")}, true);
				break;
			case LOGISTIC:
				writer.writeTag(LogisticGrowthModel.LOGISTIC_GROWTH_MODEL,
						new Attribute[]{new Attribute.Default<String>("idref", "logistic")}, true);
				break;
			case EXPANSION:
				writer.writeTag(ExpansionModel.EXPANSION_MODEL,
						new Attribute[]{new Attribute.Default<String>("idref", "expansion")}, true);
				break;
			case SKYLINE:
				writer.writeTag(BayesianSkylineLikelihood.SKYLINE_LIKELIHOOD,
						new Attribute[]{new Attribute.Default<String>("idref", "skyline")}, true);
				break;
			case GMRF_SKYRIDE:
				writer.writeTag(GMRFSkyrideLikelihood.SKYLINE_LIKELIHOOD,
						new Attribute[]{new Attribute.Default<String>("idref", "skyride")}, true);
				break;
			case YULE:
				writer.writeTag(YuleModel.YULE_MODEL,
						new Attribute[]{new Attribute.Default<String>("idref", "yule")}, true);
				break;
			case BIRTH_DEATH:
				writer.writeTag(BirthDeathGernhard08Model.BIRTH_DEATH_MODEL,
						new Attribute[]{new Attribute.Default<String>("idref", "birthDeath")}, true);
				break;
			default:
				throw new RuntimeException("No tree prior has been specified so cannot refer to it");
		}
	}

	void writeParameterLog(XMLWriter writer) {

		switch (options.nodeHeightPrior) {

			case CONSTANT:
				writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default<String>("idref", "constant.popSize"), true);
				break;
			case EXPONENTIAL:
				writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default<String>("idref", "exponential.popSize"), true);
				if (options.parameterization == ModelOptions.GROWTH_RATE) {
					writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default<String>("idref", "exponential.growthRate"), true);
				} else {
					writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default<String>("idref", "exponential.doublingTime"), true);
				}
				break;
			case LOGISTIC:
				writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default<String>("idref", "logistic.popSize"), true);
				if (options.parameterization == ModelOptions.GROWTH_RATE) {
					writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default<String>("idref", "logistic.growthRate"), true);
				} else {
					writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default<String>("idref", "logistic.doublingTime"), true);
				}
				writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default<String>("idref", "logistic.t50"), true);
				break;
			case EXPANSION:
				writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default<String>("idref", "expansion.popSize"), true);
				if (options.parameterization == ModelOptions.GROWTH_RATE) {
					writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default<String>("idref", "expansion.growthRate"), true);
				} else {
					writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default<String>("idref", "expansion.doublingTime"), true);
				}
				writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default<String>("idref", "expansion.ancestralProportion"), true);
				break;
			case SKYLINE:
				writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default<String>("idref", "skyline.popSize"), true);
				writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default<String>("idref", "skyline.groupSize"), true);
				break;
			case EXTENDED_SKYLINE:
				writeSumStatisticColumn(writer, "demographic.populationSizeChanges", "popSize_changes");
				writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default<String>("idref", "demographic.populationMean"), true);
				writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default<String>("idref", "demographic.popSize"), true);
				writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default<String>("idref", "demographic.indicators"), true);
				break;
			case GMRF_SKYRIDE:
				writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default<String>("idref", "skyride.precision"), true);
				writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default<String>("idref", "skyride.popSize"), true);
				writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default<String>("idref", "skyride.groupSize"), true);
				break;
			case YULE:
				writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default<String>("idref", "yule.birthRate"), true);
				break;
			case BIRTH_DEATH:
				writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default<String>("idref", BirthDeathModelParser.BIRTHDIFF_RATE_PARAM_NAME), true);
				writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default<String>("idref", BirthDeathModelParser.RELATIVE_DEATH_RATE_PARAM_NAME), true);
				break;
		}

	}

	void writeAnalysisToCSVfile(XMLWriter writer) {

		String logFileName = options.logFileName;

		if (options.nodeHeightPrior == TreePrior.EXTENDED_SKYLINE) {
			writer.writeOpenTag(EBSPAnalysis.VD_ANALYSIS, new Attribute[]{
					new Attribute.Default<String>("id", "demographic.analysis"),
					new Attribute.Default<Double>(EBSPAnalysis.BURN_IN, 0.1)}
			);

			writer.writeOpenTag(EBSPAnalysis.LOG_FILE_NAME);
			writer.writeText(logFileName);
			writer.writeCloseTag(EBSPAnalysis.LOG_FILE_NAME);

			writer.writeOpenTag(EBSPAnalysis.TREE_FILE_NAMES);
			writer.writeOpenTag(EBSPAnalysis.TREE_LOG);
			writer.writeText(options.treeFileName);
			writer.writeCloseTag(EBSPAnalysis.TREE_LOG);
			writer.writeCloseTag(EBSPAnalysis.TREE_FILE_NAMES);

			writer.writeOpenTag(EBSPAnalysis.MODEL_TYPE);
			writer.writeText(options.extendedSkylineModel);
			writer.writeCloseTag(EBSPAnalysis.MODEL_TYPE);

			writer.writeOpenTag(EBSPAnalysis.POPULATION_FIRST_COLUMN);
			writer.writeText(VariableDemographicModel.demoElementName + ".popSize" + 1);
			writer.writeCloseTag(EBSPAnalysis.POPULATION_FIRST_COLUMN);

			writer.writeOpenTag(EBSPAnalysis.INDICATORS_FIRST_COLUMN);
			writer.writeText(VariableDemographicModel.demoElementName + ".indicators" + 1);
			writer.writeCloseTag(EBSPAnalysis.INDICATORS_FIRST_COLUMN);

			writer.writeCloseTag(EBSPAnalysis.VD_ANALYSIS);

			writer.writeOpenTag(CSVExporter.CSV_EXPORT,
					new Attribute[]{
							new Attribute.Default<String>(CSVExporter.FILE_NAME,
									logFileName.subSequence(0, logFileName.length() - 4) + ".csv"),
							new Attribute.Default<String>(CSVExporter.SEPARATOR, ",")
					});
			writer.writeOpenTag(CSVExporter.COLUMNS);
			writer.writeTag(EBSPAnalysis.VD_ANALYSIS,
					new Attribute[]{new Attribute.Default<String>("idref", "demographic.analysis")}, true);
			writer.writeCloseTag(CSVExporter.COLUMNS);
			writer.writeCloseTag(CSVExporter.CSV_EXPORT);
		}
	}

	private void writeExponentialMarkovLikelihood(XMLWriter writer) {
		writer.writeOpenTag(
				ExponentialMarkovModel.EXPONENTIAL_MARKOV_MODEL,
				new Attribute[]{new Attribute.Default<String>("id", "eml1"),
						new Attribute.Default<String>("jeffreys", "true")}
		);
		writer.writeOpenTag(ExponentialMarkovModel.CHAIN_PARAMETER);
		writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default<String>("idref", "skyline.popSize"), true);
		writer.writeCloseTag(ExponentialMarkovModel.CHAIN_PARAMETER);
		writer.writeCloseTag(ExponentialMarkovModel.EXPONENTIAL_MARKOV_MODEL);
	}

	/**
	 * Write the boolean likelihood
	 *
	 * @param writer the writer
	 */
	private void writeBooleanLikelihood(XMLWriter writer) {
		writer.writeOpenTag(
				BooleanLikelihood.BOOLEAN_LIKELIHOOD,
				new Attribute[]{new Attribute.Default<String>("id", "booleanLikelihood1")}
		);
		writer.writeOpenTag(
				TestStatistic.TEST_STATISTIC,
				new Attribute[]{
						new Attribute.Default<String>("id", "test1"),
						new Attribute.Default<String>("name", "test1")
				}
		);
		writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default<String>("idref", "logistic.t50"), true);
		writer.writeOpenTag("lessThan");
		writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default<String>("idref", "treeModel.rootHeight"), true);
		writer.writeCloseTag("lessThan");
		writer.writeCloseTag(TestStatistic.TEST_STATISTIC);
		writer.writeCloseTag(BooleanLikelihood.BOOLEAN_LIKELIHOOD);
	}

	public void writeLikelihoodLog(XMLWriter writer) {
		if (options.nodeHeightPrior == TreePrior.YULE || options.nodeHeightPrior == TreePrior.BIRTH_DEATH) {
			writer.writeTag(SpeciationLikelihood.SPECIATION_LIKELIHOOD, new Attribute.Default<String>("idref", "speciation"), true);
		} else if (options.nodeHeightPrior == TreePrior.SKYLINE) {
			writer.writeTag(BayesianSkylineLikelihood.SKYLINE_LIKELIHOOD, new Attribute.Default<String>("idref", "skyline"), true);
		} else if (options.nodeHeightPrior == TreePrior.GMRF_SKYRIDE) {
//	        writer.writeTag(GMRFSkyrideLikelihood.SKYLINE_LIKELIHOOD, new Attribute.Default<String>("idref","skyride"), true);
			// Currently nothing additional needs logging
		} else {
			writer.writeTag(CoalescentLikelihood.COALESCENT_LIKELIHOOD, new Attribute.Default<String>("idref", "coalescent"), true);
		}

	}
}
