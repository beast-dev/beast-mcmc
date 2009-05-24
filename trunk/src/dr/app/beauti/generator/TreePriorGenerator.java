package dr.app.beauti.generator;

import dr.app.beauti.XMLWriter;
import dr.app.beauti.components.ComponentFactory;
import dr.app.beauti.options.BeautiOptions;
import dr.app.beauti.options.ModelOptions;
import dr.app.beauti.options.StartingTreeType;
import dr.app.beauti.options.TreePrior;
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
import dr.util.Attribute;
import dr.xml.XMLParser;

/**
 * @author Alexei Drummond
 */
public class TreePriorGenerator extends Generator {

	public TreePriorGenerator(BeautiOptions options, ComponentFactory[] components) {
		super(options, components);
	}

	void writeTreePrior(XMLWriter writer) {	// for species, partitionName.treeModel

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

		if (nodeHeightPrior == TreePrior.CONSTANT ) {

			writer.writeComment("A prior assumption that the population size has remained constant");
			writer.writeComment("throughout the time spanned by the genealogy.");
			writer.writeOpenTag(
					ConstantPopulationModel.CONSTANT_POPULATION_MODEL,
					new Attribute[]{
							new Attribute.Default<String>(XMLParser.ID, genePrefix + "constant"),
							new Attribute.Default<String>("units", Units.Utils.getDefaultUnitName(options.units))
					}
			);

			writer.writeOpenTag(ConstantPopulationModel.POPULATION_SIZE);
			writeParameter("constant.popSize", options, writer);
			writer.writeCloseTag(ConstantPopulationModel.POPULATION_SIZE);
			writer.writeCloseTag(ConstantPopulationModel.CONSTANT_POPULATION_MODEL);

		} else if (nodeHeightPrior == TreePrior.EXPONENTIAL) {
			// generate an exponential prior tree

			writer.writeComment("A prior assumption that the population size has grown exponentially");
			writer.writeComment("throughout the time spanned by the genealogy.");
			writer.writeOpenTag(
					ExponentialGrowthModel.EXPONENTIAL_GROWTH_MODEL,
					new Attribute[]{
							new Attribute.Default<String>(XMLParser.ID, genePrefix + "exponential"),
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
		} else if (nodeHeightPrior == TreePrior.LOGISTIC) {
			// generate an exponential prior tree

			writer.writeComment("A prior assumption that the population size has grown logistically");
			writer.writeComment("throughout the time spanned by the genealogy.");
			writer.writeOpenTag(
					LogisticGrowthModel.LOGISTIC_GROWTH_MODEL,
					new Attribute[]{
							new Attribute.Default<String>(XMLParser.ID, genePrefix + "logistic"),
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

		} else if (nodeHeightPrior == TreePrior.EXPANSION) {
			// generate an exponential prior tree

			writer.writeComment("A prior assumption that the population size has grown exponentially");
			writer.writeComment("from some ancestral population size in the past.");
			writer.writeOpenTag(
					ExpansionModel.EXPANSION_MODEL,
					new Attribute[]{
							new Attribute.Default<String>(XMLParser.ID, genePrefix + "expansion"),
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

		} else if (nodeHeightPrior == TreePrior.YULE) {
			writer.writeComment("A prior on the distribution node heights defined given");
			writer.writeComment("a Yule speciation process (a pure birth process).");
			writer.writeOpenTag(
					YuleModel.YULE_MODEL,
					new Attribute[]{
							new Attribute.Default<String>(XMLParser.ID, genePrefix + "yule"),
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
							new Attribute.Default<String>(XMLParser.ID, genePrefix + "birthDeath"),
							new Attribute.Default<String>("units", Units.Utils.getDefaultUnitName(units))
					}
			);

			writeParameter(BirthDeathModelParser.BIRTHDIFF_RATE, genePrefix + BirthDeathModelParser.BIRTHDIFF_RATE_PARAM_NAME, options, writer);
			writeParameter(BirthDeathModelParser.RELATIVE_DEATH_RATE, genePrefix + BirthDeathModelParser.RELATIVE_DEATH_RATE_PARAM_NAME, options, writer);
			writer.writeCloseTag(BirthDeathGernhard08Model.BIRTH_DEATH_MODEL);
		} else if (nodeHeightPrior == TreePrior.SPECIES_BIRTH_DEATH ){
//				|| nodeHeightPrior == TreePrior.SPECIES_YULE) {
			
				writer.writeComment("A prior assumption that the population size has remained constant");
				writer.writeComment("throughout the time spanned by the genealogy.");
				writer.writeOpenTag(
						ConstantPopulationModel.CONSTANT_POPULATION_MODEL,
						new Attribute[]{
								new Attribute.Default<String>(XMLParser.ID, genePrefix + "constant"),
								new Attribute.Default<String>("units", Units.Utils.getDefaultUnitName(options.units))
						}
				);

				writer.writeOpenTag(ConstantPopulationModel.POPULATION_SIZE);
				writeParameter("constant.popSize", options, writer);
				writer.writeCloseTag(ConstantPopulationModel.POPULATION_SIZE);
				writer.writeCloseTag(ConstantPopulationModel.CONSTANT_POPULATION_MODEL);
		}
		
		if ((!options.isSpeciesAnalysis()) && nodeHeightPrior != TreePrior.CONSTANT && nodeHeightPrior != TreePrior.EXPONENTIAL) {
			// If the node height prior is not one of these two then we need to simulate a
			// random starting tree under a constant size coalescent.

			writer.writeComment("This is a simple constant population size coalescent model");
			writer.writeComment("that is used to generate an initial tree for the chain.");
			writer.writeOpenTag(
					ConstantPopulationModel.CONSTANT_POPULATION_MODEL,
					new Attribute[]{
							new Attribute.Default<String>(XMLParser.ID, genePrefix + "initialDemo"),
							new Attribute.Default<String>("units", Units.Utils.getDefaultUnitName(units))
					}
			);

			writer.writeOpenTag(ConstantPopulationModel.POPULATION_SIZE);
			if (initialPopSize != null) {
				writer.writeTag(ParameterParser.PARAMETER,
						new Attribute[]{
								new Attribute.Default<String>(XMLParser.IDREF, genePrefix + initialPopSize),
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
			writer.writeComment("Generate a speciational process");
			writer.writeOpenTag(
					SpeciationLikelihood.SPECIATION_LIKELIHOOD,
					new Attribute[]{
							new Attribute.Default<String>(XMLParser.ID, genePrefix + "speciation")
					}
			);

			// write pop size socket
			writer.writeOpenTag(SpeciationLikelihood.MODEL);
			writeNodeHeightPriorModelRef(writer);
			writer.writeCloseTag(SpeciationLikelihood.MODEL);
			writer.writeOpenTag(SpeciationLikelihood.TREE);
			writer.writeTag(TreeModel.TREE_MODEL, new Attribute.Default<String>(XMLParser.IDREF, genePrefix + TreeModel.TREE_MODEL), true);
			writer.writeCloseTag(SpeciationLikelihood.TREE);

			writer.writeCloseTag(SpeciationLikelihood.SPECIATION_LIKELIHOOD);

		} else if (treePrior == TreePrior.SKYLINE) {
			// generate a Bayesian skyline plot
			writer.writeComment("Generate a Bayesian skyline process");
			writer.writeOpenTag(
					BayesianSkylineLikelihood.SKYLINE_LIKELIHOOD,
					new Attribute[]{
							new Attribute.Default<String>(XMLParser.ID, genePrefix + "skyline"),
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
			writer.writeTag(TreeModel.TREE_MODEL, new Attribute.Default<String>(XMLParser.IDREF, genePrefix + TreeModel.TREE_MODEL), true);
			writer.writeCloseTag(CoalescentLikelihood.POPULATION_TREE);

			writer.writeCloseTag(BayesianSkylineLikelihood.SKYLINE_LIKELIHOOD);
		} else if (treePrior == TreePrior.EXTENDED_SKYLINE) {
			final String tagName = VariableDemographicModel.PARSER.getParserName();
			writer.writeComment("Generate an extended Bayesian skyline process");
			writer.writeOpenTag(
					tagName,
					new Attribute[]{
							new Attribute.Default<String>(XMLParser.ID, genePrefix + VariableDemographicModel.demoElementName),
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
			writer.writeTag(TreeModel.TREE_MODEL, new Attribute.Default<String>(XMLParser.IDREF, genePrefix + TreeModel.TREE_MODEL), true);
			writer.writeCloseTag(VariableDemographicModel.POP_TREE);

			writer.writeCloseTag(VariableDemographicModel.POPULATION_TREES);

			writer.writeCloseTag(tagName);

			writer.writeOpenTag(CoalescentLikelihood.COALESCENT_LIKELIHOOD, new Attribute.Default<String>(XMLParser.ID, genePrefix + COALESCENT));
			writer.writeOpenTag(CoalescentLikelihood.MODEL);
			writer.writeTag(tagName, new Attribute.Default<String>(XMLParser.IDREF, genePrefix + VariableDemographicModel.demoElementName), true);
			writer.writeCloseTag(CoalescentLikelihood.MODEL);
			writer.writeComment("Take population Tree from demographic");
			writer.writeCloseTag(CoalescentLikelihood.COALESCENT_LIKELIHOOD);

			writer.writeOpenTag(SumStatistic.SUM_STATISTIC,
					new Attribute[]{
							new Attribute.Default<String>(XMLParser.ID, genePrefix + VariableDemographicModel.demoElementName + ".populationSizeChanges"),
							new Attribute.Default<String>("elementwise", "true")
					});
			writer.writeTag(ParameterParser.PARAMETER,
					new Attribute.Default<String>(XMLParser.IDREF, genePrefix + VariableDemographicModel.demoElementName + ".indicators"), true);
			writer.writeCloseTag(SumStatistic.SUM_STATISTIC);
			writer.writeOpenTag(ExponentialDistributionModel.EXPONENTIAL_DISTRIBUTION_MODEL,
					new Attribute[]{
							new Attribute.Default<String>(XMLParser.ID, genePrefix + VariableDemographicModel.demoElementName + ".populationMeanDist")
							//,new Attribute.Default<String>("elementwise", "true")
					});
			writer.writeOpenTag(ExponentialDistributionModel.MEAN);
			writer.writeTag(ParameterParser.PARAMETER,
					new Attribute[]{
							new Attribute.Default<String>(XMLParser.ID, genePrefix + VariableDemographicModel.demoElementName + ".populationMean"),
							new Attribute.Default<String>("value", "1")}, true);
			writer.writeCloseTag(ExponentialDistributionModel.MEAN);
			writer.writeCloseTag(ExponentialDistributionModel.EXPONENTIAL_DISTRIBUTION_MODEL);
		} else if (treePrior == TreePrior.GMRF_SKYRIDE) {
			writer.writeComment("Generate a GMRF Bayesian Skyride process");
			writer.writeOpenTag(
					GMRFSkyrideLikelihood.SKYLINE_LIKELIHOOD,
					new Attribute[]{
							new Attribute.Default<String>(XMLParser.ID, genePrefix + "skyride"),
							new Attribute.Default<String>(GMRFSkyrideLikelihood.TIME_AWARE_SMOOTHING,
									options.skyrideSmoothing == ModelOptions.SKYRIDE_TIME_AWARE_SMOOTHING ? "true" : "false"),
                            new Attribute.Default<String>(GMRFSkyrideLikelihood.RANDOMIZE_TREE,
                                    options.startingTreeType == StartingTreeType.UPGMA ? "true" : "false"),
					}
			);

            int skyrideIntervalCount = options.taxonList.getTaxonCount() - 1;
			writer.writeOpenTag(GMRFSkyrideLikelihood.POPULATION_PARAMETER);
			writeParameter("skyride.popSize", skyrideIntervalCount, writer);
			writer.writeCloseTag(GMRFSkyrideLikelihood.POPULATION_PARAMETER);

			writer.writeOpenTag(GMRFSkyrideLikelihood.GROUP_SIZES);
			writeParameter("skyride.groupSize", skyrideIntervalCount, writer);
			writer.writeCloseTag(GMRFSkyrideLikelihood.GROUP_SIZES);

			writer.writeOpenTag(GMRFSkyrideLikelihood.PRECISION_PARAMETER);
			writeParameter("skyride.precision", 1, writer);
			writer.writeCloseTag(GMRFSkyrideLikelihood.PRECISION_PARAMETER);

			writer.writeOpenTag(GMRFSkyrideLikelihood.POPULATION_TREE);
			writer.writeTag(TreeModel.TREE_MODEL, new Attribute.Default<String>(XMLParser.IDREF, genePrefix + TreeModel.TREE_MODEL), true);
			writer.writeCloseTag(GMRFSkyrideLikelihood.POPULATION_TREE);

			writer.writeCloseTag(GMRFSkyrideLikelihood.SKYLINE_LIKELIHOOD);
		
		} else if (options.isSpeciesAnalysis()) {
//			writer.writeComment("Gene tree prior uses speices tree prior.");
		} else {
			// generate a coalescent process
			writer.writeComment("Generate a coalescent process");
			writer.writeOpenTag(
					CoalescentLikelihood.COALESCENT_LIKELIHOOD,
					new Attribute[]{new Attribute.Default<String>(XMLParser.ID, genePrefix + COALESCENT)}
			);
			writer.writeOpenTag(CoalescentLikelihood.MODEL);
			writeNodeHeightPriorModelRef(writer);
			writer.writeCloseTag(CoalescentLikelihood.MODEL);
			writer.writeOpenTag(CoalescentLikelihood.POPULATION_TREE);
			writer.writeTag(TreeModel.TREE_MODEL, new Attribute.Default<String>(XMLParser.IDREF, genePrefix + TreeModel.TREE_MODEL), true);
			writer.writeCloseTag(CoalescentLikelihood.POPULATION_TREE);
			writer.writeCloseTag(CoalescentLikelihood.COALESCENT_LIKELIHOOD);
		}
	}

	void writeNodeHeightPriorModelRef(XMLWriter writer) {

		TreePrior treePrior = options.nodeHeightPrior;

		switch (treePrior) {
			case CONSTANT:
//			case SPECIES_YULE:
			case SPECIES_BIRTH_DEATH:
				writer.writeIDref(ConstantPopulationModel.CONSTANT_POPULATION_MODEL, genePrefix + "constant");
				break;
			case EXPONENTIAL:
				writer.writeIDref(ExponentialGrowthModel.EXPONENTIAL_GROWTH_MODEL, genePrefix + "exponential");
				break;
			case LOGISTIC:
				writer.writeIDref(LogisticGrowthModel.LOGISTIC_GROWTH_MODEL, genePrefix + "logistic");
				break;
			case EXPANSION:
				writer.writeIDref(ExpansionModel.EXPANSION_MODEL, genePrefix + "expansion");
				break;
			case SKYLINE:
				writer.writeIDref(BayesianSkylineLikelihood.SKYLINE_LIKELIHOOD, genePrefix + "skyline");
				break;
			case GMRF_SKYRIDE:
				writer.writeIDref(GMRFSkyrideLikelihood.SKYLINE_LIKELIHOOD, genePrefix + "skyride");
				break;
			case YULE:
				writer.writeIDref(YuleModel.YULE_MODEL, genePrefix + "yule");
				break;
			case BIRTH_DEATH:
				writer.writeIDref(BirthDeathGernhard08Model.BIRTH_DEATH_MODEL, genePrefix + "birthDeath");
				break;
			default:
				throw new RuntimeException("No tree prior has been specified so cannot refer to it");
		}
	}

	void writeParameterLog(XMLWriter writer) {

		switch (options.nodeHeightPrior) {

			case CONSTANT:
//			case SPECIES_YULE:
			case SPECIES_BIRTH_DEATH:
				writer.writeIDref(ParameterParser.PARAMETER, genePrefix + "constant.popSize");
				break;
			case EXPONENTIAL:
				writer.writeIDref(ParameterParser.PARAMETER, genePrefix + "exponential.popSize");
				if (options.parameterization == ModelOptions.GROWTH_RATE) {
					writer.writeIDref(ParameterParser.PARAMETER, genePrefix + "exponential.growthRate");
				} else {
					writer.writeIDref(ParameterParser.PARAMETER, genePrefix + "exponential.doublingTime");
				}
				break;
			case LOGISTIC:
				writer.writeIDref(ParameterParser.PARAMETER, genePrefix + "logistic.popSize");
				if (options.parameterization == ModelOptions.GROWTH_RATE) {
					writer.writeIDref(ParameterParser.PARAMETER, genePrefix + "logistic.growthRate");
				} else {
					writer.writeIDref(ParameterParser.PARAMETER, genePrefix + "logistic.doublingTime");
				}
				writer.writeIDref(ParameterParser.PARAMETER, genePrefix + "logistic.t50");
				break;
			case EXPANSION:
				writer.writeIDref(ParameterParser.PARAMETER, genePrefix + "expansion.popSize");
				if (options.parameterization == ModelOptions.GROWTH_RATE) {
					writer.writeIDref(ParameterParser.PARAMETER, genePrefix + "expansion.growthRate");
				} else {
					writer.writeIDref(ParameterParser.PARAMETER, genePrefix + "expansion.doublingTime");
				}
				writer.writeIDref(ParameterParser.PARAMETER, genePrefix + "expansion.ancestralProportion");
				break;
			case SKYLINE:
				writer.writeIDref(ParameterParser.PARAMETER, genePrefix + "skyline.popSize");
				writer.writeIDref(ParameterParser.PARAMETER, genePrefix + "skyline.groupSize");
				break;
			case EXTENDED_SKYLINE:
				writeSumStatisticColumn(writer, "demographic.populationSizeChanges", "popSize_changes");
				writer.writeIDref(ParameterParser.PARAMETER, genePrefix + "demographic.populationMean");
				writer.writeIDref(ParameterParser.PARAMETER, genePrefix + "demographic.popSize");
				writer.writeIDref(ParameterParser.PARAMETER, genePrefix + "demographic.indicators");
				break;
			case GMRF_SKYRIDE:
				writer.writeIDref(ParameterParser.PARAMETER, genePrefix + "skyride.precision");
				writer.writeIDref(ParameterParser.PARAMETER, genePrefix + "skyride.popSize");
				writer.writeIDref(ParameterParser.PARAMETER, genePrefix + "skyride.groupSize");
				break;
			case YULE:
				writer.writeIDref(ParameterParser.PARAMETER, genePrefix + "yule.birthRate");
				break;
			case BIRTH_DEATH:
				writer.writeIDref(ParameterParser.PARAMETER, genePrefix + BirthDeathModelParser.BIRTHDIFF_RATE_PARAM_NAME);
				writer.writeIDref(ParameterParser.PARAMETER, genePrefix + BirthDeathModelParser.RELATIVE_DEATH_RATE_PARAM_NAME);
				break;
		}

	}

	void writeAnalysisToCSVfile(XMLWriter writer) {

		String logFileName = options.logFileName;

		if (options.nodeHeightPrior == TreePrior.EXTENDED_SKYLINE) {
			writer.writeOpenTag(EBSPAnalysis.VD_ANALYSIS, new Attribute[]{
					new Attribute.Default<String>(XMLParser.ID, genePrefix + "demographic.analysis"),
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
			writer.writeIDref(EBSPAnalysis.VD_ANALYSIS, genePrefix + "demographic.analysis");
			writer.writeCloseTag(CSVExporter.COLUMNS);
			writer.writeCloseTag(CSVExporter.CSV_EXPORT);
		}
	}

	private void writeExponentialMarkovLikelihood(XMLWriter writer) {
		writer.writeOpenTag(
				ExponentialMarkovModel.EXPONENTIAL_MARKOV_MODEL,
				new Attribute[]{new Attribute.Default<String>(XMLParser.ID, genePrefix + "eml1"),
						new Attribute.Default<String>("jeffreys", "true")}
		);
		writer.writeOpenTag(ExponentialMarkovModel.CHAIN_PARAMETER);
		writer.writeIDref(ParameterParser.PARAMETER, genePrefix + "skyline.popSize");
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
				new Attribute[]{new Attribute.Default<String>(XMLParser.ID, genePrefix + "booleanLikelihood1")}
		);
		writer.writeOpenTag(
				TestStatistic.TEST_STATISTIC,
				new Attribute[]{
						new Attribute.Default<String>(XMLParser.ID, "test1"),
						new Attribute.Default<String>("name", "test1")
				}
		);
		writer.writeIDref(ParameterParser.PARAMETER, genePrefix + "logistic.t50");
		writer.writeOpenTag("lessThan");
		writer.writeIDref(ParameterParser.PARAMETER, genePrefix + "treeModel.rootHeight");
		writer.writeCloseTag("lessThan");
		writer.writeCloseTag(TestStatistic.TEST_STATISTIC);
		writer.writeCloseTag(BooleanLikelihood.BOOLEAN_LIKELIHOOD);
	}

	public void writeLikelihoodLog(XMLWriter writer) {
		if (options.nodeHeightPrior == TreePrior.YULE || options.nodeHeightPrior == TreePrior.BIRTH_DEATH) {
			writer.writeIDref(SpeciationLikelihood.SPECIATION_LIKELIHOOD, genePrefix + "speciation");
		} else if (options.nodeHeightPrior == TreePrior.SKYLINE) {
			writer.writeIDref(BayesianSkylineLikelihood.SKYLINE_LIKELIHOOD, genePrefix + "skyline");
		} else if (options.nodeHeightPrior == TreePrior.GMRF_SKYRIDE) {
//	        writer.writeIDref(GMRFSkyrideLikelihood.SKYLINE_LIKELIHOOD, "skyride");
			// Currently nothing additional needs logging
		} else if (options.isSpeciesAnalysis()) {
			// no
		} else {
			writer.writeIDref(CoalescentLikelihood.COALESCENT_LIKELIHOOD, genePrefix + COALESCENT);
		}

	}
}
