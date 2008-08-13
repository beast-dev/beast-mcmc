/*
 * BeastOptions.java
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

package dr.app.tools;

import dr.evolution.alignment.Alignment;
import dr.evolution.alignment.SimpleAlignment;
import dr.evolution.alignment.SitePatterns;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.util.Date;
import dr.evolution.util.Taxon;
import dr.evolution.util.TaxonList;
import dr.evolution.util.Units;
import dr.evomodel.arg.coalescent.VeryOldCoalescentLikelihood;
import dr.evomodel.coalescent.CoalescentSimulator;
import dr.evomodel.coalescent.ConstantPopulationModel;
import dr.evomodel.coalescent.ExponentialGrowthModel;
import dr.evomodel.coalescent.LogisticGrowthModel;
import dr.evomodel.operators.ExchangeOperator;
import dr.evomodel.operators.SubtreeSlideOperator;
import dr.evomodel.operators.WilsonBalding;
import dr.evomodel.sitemodel.GammaSiteModel;
import dr.evomodel.substmodel.EmpiricalAminoAcidModel;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.tree.TreeLogger;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.treelikelihood.TreeLikelihood;
import dr.evoxml.SitePatternsParser;
import dr.evoxml.TaxaParser;
import dr.inference.loggers.Columns;
import dr.inference.loggers.MCLogger;
import dr.inference.model.*;
import dr.inference.operators.*;
import dr.util.Attribute;

import java.io.Writer;
import java.util.ArrayList;

public class BeastOptions {

	public static final int JC = 0;
	public static final int HKY = 1;
	public static final int GTR = 2;

	public static final int BLOSUM_62 = 0;
	public static final int DAYHOFF = 1;
	public static final int JTT = 2;
	public static final int MT_REV_24 = 3;
	public static final int CP_REV_45 = 4;
	public static final int WAG = 5;

	public static final int NONE = -1;
	public static final int CONSTANT = 0;
	public static final int EXPONENTIAL = 1;
	public static final int LOGISTIC = 2;

	public BeastOptions() {
	}

	public String fileNameStem = null;
	public String logFileName = null;
	public String treeFileName = null;

	public TaxonList taxonList = null;
	public SimpleAlignment originalAlignment = null;
	public Alignment alignment = null;
	public Tree tree = null;

	// MCMC options
	public int chainLength = 100000;
	public int logEvery = 100;
	public int echoEvery = 100;
	public int burnIn = 1000;
	public String fileName = null;
	public boolean userTree = false;

	public int nucSubstitutionModel = HKY;
	public int aaSubstitutionModel = BLOSUM_62;
	public boolean gammaHetero = false;
	public int gammaCategories = 4;
	public boolean invarHetero = false;
	public boolean codonHetero = false;

	public double meanMutationRate = 1.0;

	public boolean unlinkedSubstitutionModel = false;
	public boolean unlinkedHeterogeneityModel = false;
	public boolean unlinkedFrequencyModel = false;

	public int coalescentModel = CONSTANT;
	public boolean fixedTree = false;

	//Operators (ends in O for the tweak value, W for the weight)
	//Numbered ones are for codon positions

	public double popSizeO = 0.5;
	public double popSizeW = 1.0;

	public double growthRateO = 0.1;
	public double growthRateW = 1.0;

	public double logisticShapeO = 0.1;
	public double logisticShapeW = 1.0;

	public double kappaO = 0.5;
	public double kappaW = 1.0;
	public double kappa1O = 0.5;
	public double kappa1W = 1.0;
	public double kappa2O = 0.5;
	public double kappa2W = 1.0;
	public double kappa3O = 0.5;
	public double kappa3W = 1.0;

	public double gtrACO = 0.5;
	public double gtrACW = 1.0;
	public double gtrAC1O = 0.5;
	public double gtrAC1W = 1.0;
	public double gtrAC2O = 0.5;
	public double gtrAC2W = 1.0;
	public double gtrAC3O = 0.5;
	public double gtrAC3W = 1.0;

	public double gtrAGO = 0.5;
	public double gtrAGW = 1.0;
	public double gtrAG1O = 0.5;
	public double gtrAG1W = 1.0;
	public double gtrAG2O = 0.5;
	public double gtrAG2W = 1.0;
	public double gtrAG3O = 0.5;
	public double gtrAG3W = 1.0;

	public double gtrATO = 0.5;
	public double gtrATW = 1.0;
	public double gtrAT1O = 0.5;
	public double gtrAT1W = 1.0;
	public double gtrAT2O = 0.5;
	public double gtrAT2W = 1.0;
	public double gtrAT3O = 0.5;
	public double gtrAT3W = 1.0;

	public double gtrCGO = 0.5;
	public double gtrCGW = 1.0;
	public double gtrCG1O = 0.5;
	public double gtrCG1W = 1.0;
	public double gtrCG2O = 0.5;
	public double gtrCG2W = 1.0;
	public double gtrCG3O = 0.5;
	public double gtrCG3W = 1.0;

	public double gtrGTO = 0.5;
	public double gtrGTW = 1.0;
	public double gtrGT1O = 0.5;
	public double gtrGT1W = 1.0;
	public double gtrGT2O = 0.5;
	public double gtrGT2W = 1.0;
	public double gtrGT3O = 0.5;
	public double gtrGT3W = 1.0;

	public double alphaO = 0.5;
	public double alphaW = 1.0;
	public double alpha1O = 0.5;
	public double alpha1W = 1.0;
	public double alpha2O = 0.5;
	public double alpha2W = 1.0;
	public double alpha3O = 0.5;
	public double alpha3W = 1.0;

	public double pinvO = 0.5;
	public double pinvW = 1.0;
	public double pinv1O = 0.5;
	public double pinv1W = 1.0;
	public double pinv2O = 0.5;
	public double pinv2W = 1.0;
	public double pinv3O = 0.5;
	public double pinv3W = 1.0;

	public double muO = 0.5;
	public double muW = 1.0;
	public double mu1O = 0.5;
	public double mu1W = 1.0;
	public double mu2O = 0.5;
	public double mu2W = 1.0;
	public double mu3O = 0.5;
	public double mu3W = 1.0;

	public double muHeightsO = 0.1;
	public double muHeightsW = 2.0;

	public double musCenteredO = 0.5;
	public double musCenteredW = 1.0;
	public double musDeltaO = 0.5;
	public double musDeltaW = 1.0;

	public double rootO = 0.5;
	public double rootW = 1.0;

	public double heightsW = 2.0;

	public double wideExchangeW = 2.0;
	public double narrowExchangeW = 2.0;
	public double wilsonBaldingW = 2.0;

	public double subtreeSlideO = 1.0;
	public double subtreeSlideW = 2.0;

	//Miscellaneous

	public int units = Units.SUBSTITUTIONS;
	public boolean usingDates = false;

	/**
	 * Generate a beast xml file from these beast options
	 */
	public void generateXML(Writer w) {

		XMLWriter writer = new XMLWriter(w);

		writer.writeText("<?xml version=\"1.0\" standalone=\"yes\"?>");
		writer.writeComment("Generated by BEAUTi v1.1");
		writer.writeComment("Copyright 2002-2004 Alexei Drummond and Andrew Rambaut");
		writer.writeComment("      Department of Zoology, University of Oxford");
		writer.writeComment("      http://evolve.zoo.ox.ac.uk/beast/");
		writer.writeOpenTag("beast");
		writer.writeText("");
		writeTaxa(writer);
		if (alignment != null) {
			writer.writeText("");
			writeAlignment(writer);
			writer.writeText("");
			writePatternList(writer);

		}
		writer.writeText("");
		writeDemographicModel(writer);
		writer.writeText("");
		writeStartingTree(writer);
		writer.writeText("");
		writeTreeModel(writer);
		if (coalescentModel != NONE) {
			writer.writeText("");
			writeCoalescentLikelihood(writer);
		}
		if (coalescentModel == LOGISTIC) {
			writer.writeText("");
			writeBooleanLikelihood(writer);
		}
		if (alignment != null) {
			writer.writeText("");
			writeSubstitutionModel(writer);
			writer.writeText("");
			writeSiteModel(writer);
			writer.writeText("");
			writeTreeLikelihood(writer);
		}
		writer.writeText("");
		writeOperatorSchedule(writer);
		writer.writeText("");
		writeMCMC(writer);
		writer.writeText("");
		writeTimerReport(writer);
		writer.writeText("");
		writeTraceAnalysis(writer);
		writer.writeCloseTag("beast");
		writer.flush();
	}

	/**
	 * Generate a taxa block from these beast options
	 */
	public void writeTaxa(XMLWriter writer) {

		writer.writeComment("The list of taxa analyse (can also include dates/ages).");
		writer.writeComment("ntax=" + taxonList.getTaxonCount());
		writer.writeOpenTag("taxa", new Attribute[]{new Attribute.Default("id", "taxa1")});

		boolean firstDate = true;
		for (int i = 0; i < taxonList.getTaxonCount(); i++) {
			Taxon taxon = taxonList.getTaxon(i);

			boolean hasDate = false;

			if (usingDates) {
				hasDate = TaxonList.Utils.hasAttribute(taxonList, i, Date.DATE);
			}

			writer.writeTag("taxon", new Attribute[]{new Attribute.Default("id", taxon.getId())}, !hasDate);

			if (hasDate) {
				Date date = (Date) taxon.getAttribute(Date.DATE);

				if (firstDate) {
					units = date.getUnits();
					firstDate = false;
				} else {
					if (units != date.getUnits()) {
						System.err.println("Error: Units in dates do not match.");
					}
				}

				Attribute[] attributes = new Attribute[]{
						new Attribute.Default("value", date.getTimeValue() + ""),
						new Attribute.Default("direction", date.isBackwards() ? "backwards" : "forwards"),
						new Attribute.Default("units", Units.Utils.getDefaultUnitName(units))
						/*,
											new Attribute.Default("origin", date.getOrigin()+"")*/
				};

				writer.writeTag(Date.DATE, attributes, true);
				writer.writeCloseTag("taxon");
			}
		}

		writer.writeCloseTag("taxa");
	}

	/**
	 * Generate an alignment block from these beast options
	 */
	public void writeAlignment(XMLWriter writer) {

		writer.writeComment("The sequence alignment (each sequence refers to a taxon above).");
		writer.writeComment("ntax=" + alignment.getTaxonCount() + " nchar=" + alignment.getSiteCount());
		writer.writeOpenTag(
				"alignment",
				new Attribute[]{
						new Attribute.Default("id", "alignment1"),
						new Attribute.Default("dataType", alignment.getDataType().getDescription())
				}
		);

		for (int i = 0; i < alignment.getTaxonCount(); i++) {
			Taxon taxon = alignment.getTaxon(i);

			writer.writeOpenTag("sequence");
			writer.writeTag("taxon", new Attribute[]{new Attribute.Default("idref", taxon.getId())}, true);
			writer.writeText(alignment.getAlignedSequenceString(i));
			writer.writeCloseTag("sequence");
		}
		writer.writeCloseTag("alignment");
	}

	/**
	 * Write a demographic model
	 */
	public void writeDemographicModel(XMLWriter writer) {
		if (coalescentModel == CONSTANT) {

			writer.writeComment("A prior assumption that the population size has remained constant");
			writer.writeComment("throughout the time spanned by the genealogy.");
			writer.writeOpenTag(
					ConstantPopulationModel.CONSTANT_POPULATION_MODEL,
					new Attribute[]{
							new Attribute.Default("id", "demo1"),
							new Attribute.Default("units", Units.Utils.getDefaultUnitName(units))
					}
			);

			writer.writeOpenTag(ConstantPopulationModel.POPULATION_SIZE);
			writeParameter("demo1.popSize", 1.0, 0.0, Double.MAX_VALUE, writer);
			writer.writeCloseTag(ConstantPopulationModel.POPULATION_SIZE);
			writer.writeCloseTag(ConstantPopulationModel.CONSTANT_POPULATION_MODEL);

		} else if (coalescentModel == EXPONENTIAL) {
			// generate an exponential prior tree

			writer.writeComment("A prior assumption that the population size has grown exponentially");
			writer.writeComment("throughout the time spanned by the genealogy.");
			writer.writeOpenTag(
					ExponentialGrowthModel.EXPONENTIAL_GROWTH_MODEL,
					new Attribute[]{
							new Attribute.Default("id", "demo1"),
							new Attribute.Default("units", Units.Utils.getDefaultUnitName(units))
					}
			);

			// write pop size socket
			writer.writeOpenTag(ExponentialGrowthModel.POPULATION_SIZE);
			writeParameter("demo1.popSize", 1.0, 0.0, Double.MAX_VALUE, writer);
			writer.writeCloseTag(ExponentialGrowthModel.POPULATION_SIZE);

			// write growth rate socket
			writer.writeOpenTag(ExponentialGrowthModel.GROWTH_RATE);
			writeParameter("demo1.growthRate", 0.0, -Double.MAX_VALUE, Double.MAX_VALUE, writer);
			writer.writeCloseTag(ExponentialGrowthModel.GROWTH_RATE);

			writer.writeCloseTag(ExponentialGrowthModel.EXPONENTIAL_GROWTH_MODEL);
		} else if (coalescentModel == LOGISTIC) {
			// generate an exponential prior tree

			writer.writeComment("A prior assumption that the population size has grown logistically");
			writer.writeComment("throughout the time spanned by the genealogy.");
			writer.writeOpenTag(
					LogisticGrowthModel.LOGISTIC_GROWTH_MODEL,
					new Attribute[]{
							new Attribute.Default("id", "demo1"),
							new Attribute.Default("units", Units.Utils.getDefaultUnitName(units))
					}
			);

			// write pop size socket
			writer.writeOpenTag(LogisticGrowthModel.POPULATION_SIZE);
			writeParameter("demo1.popSize", 1.0, 0.0, Double.MAX_VALUE, writer);
			writer.writeCloseTag(LogisticGrowthModel.POPULATION_SIZE);

			// write growth rate socket
			writer.writeOpenTag(LogisticGrowthModel.DOUBLING_TIME);
			writeParameter("demo1.doublingTime", 1.0, 0.0, Double.MAX_VALUE, writer);
			writer.writeCloseTag(LogisticGrowthModel.DOUBLING_TIME);

			// write logistic t50 socket
			writer.writeOpenTag(LogisticGrowthModel.TIME_50);
			writeParameter("demo1.t50", 0.001, 0.0, Double.MAX_VALUE, writer);
			writer.writeCloseTag(LogisticGrowthModel.TIME_50);

			writer.writeCloseTag(LogisticGrowthModel.LOGISTIC_GROWTH_MODEL);

			writer.writeText("");
			writer.writeComment("This is a simple constant population size coalescent model");
			writer.writeComment("that is used to generate an initial tree for the chain.");
			writer.writeOpenTag(
					ConstantPopulationModel.CONSTANT_POPULATION_MODEL,
					new Attribute[]{
							new Attribute.Default("id", "initialDemo"),
							new Attribute.Default("units", Units.Utils.getDefaultUnitName(units))
					}
			);

			writer.writeOpenTag(ConstantPopulationModel.POPULATION_SIZE);
			writer.writeTag(ParameterParser.PARAMETER,
					new Attribute[]{
							new Attribute.Default("idref", "demo1.popSize"),
					}, true);
			writer.writeCloseTag(ConstantPopulationModel.POPULATION_SIZE);
			writer.writeCloseTag(ConstantPopulationModel.CONSTANT_POPULATION_MODEL);

		} else if (coalescentModel == NONE) {

			writer.writeComment("This is a simple constant population size coalescent model");
			writer.writeComment("that is used to generate an initial tree for the chain.");
			writer.writeOpenTag(
					ConstantPopulationModel.CONSTANT_POPULATION_MODEL,
					new Attribute[]{
							new Attribute.Default("id", "initialDemo"),
							new Attribute.Default("units", Units.Utils.getDefaultUnitName(units))
					}
			);

			writer.writeOpenTag(ConstantPopulationModel.POPULATION_SIZE);
			writeParameter("initialDemo.popSize", 1.0, 0.0, Double.MAX_VALUE, writer);
			writer.writeCloseTag(ConstantPopulationModel.POPULATION_SIZE);
			writer.writeCloseTag(ConstantPopulationModel.CONSTANT_POPULATION_MODEL);

		}
	}

	/**
	 * Writes the pattern lists
	 */
	public void writePatternList(XMLWriter writer) {

		if (alignment.getDataType() == Nucleotides.INSTANCE && codonHetero) {
			for (int i = 1; i <= 3; i++) {
				writePatternList(i, 3, writer);
			}
		} else {
			writePatternList(1, 0, writer);
		}
	}

	/**
	 * Write a single pattern list
	 */
	private void writePatternList(int from, int every, XMLWriter writer) {

		SitePatterns patterns = new SitePatterns(alignment, from, 0, every);
		writer.writeComment("The unique patterns " + ((every == 3) ? " in codon position " + from : ""));
		writer.writeComment("npatterns=" + patterns.getPatternCount());
		if (every != 0) {
			writer.writeOpenTag(SitePatternsParser.PATTERNS,
					new Attribute[]{
							new Attribute.Default("id", "patterns" + from),
							new Attribute.Default("from", "" + from),
							new Attribute.Default("every", "" + every)
					}
			);
		} else {
			writer.writeOpenTag(SitePatternsParser.PATTERNS,
					new Attribute[]{
							new Attribute.Default("id", "patterns" + from),
							new Attribute.Default("from", "" + from)
					}
			);
		}

		writer.writeTag("alignment", new Attribute.Default("idref", "alignment1"), true);
		writer.writeCloseTag(SitePatternsParser.PATTERNS);
	}

	/**
	 * Write tree model XML block.
	 */
	private void writeTreeModel(XMLWriter writer) {

		writer.writeTag(TreeModel.TREE_MODEL, new Attribute.Default("id", "treeModel1"), false);

		if (userTree) {
			writer.writeTag("tree", new Attribute.Default("idref", "startingTree"), true);
		} else {
			writer.writeTag(CoalescentSimulator.COALESCENT_TREE, new Attribute.Default("idref", "startingTree"), true);
		}

		writer.writeOpenTag(TreeModel.ROOT_HEIGHT);
		writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default("id", "treeModel1.rootHeight"), true);
		writer.writeCloseTag(TreeModel.ROOT_HEIGHT);


		writer.writeOpenTag(TreeModel.NODE_HEIGHTS, new Attribute.Default(TreeModel.INTERNAL_NODES, "true"));
		writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default("id", "treeModel1.internalNodeHeights"), true);
		writer.writeCloseTag(TreeModel.NODE_HEIGHTS);

		writer.writeOpenTag(TreeModel.NODE_HEIGHTS,
				new Attribute[]{
						new Attribute.Default(TreeModel.INTERNAL_NODES, "true"),
						new Attribute.Default(TreeModel.ROOT_NODE, "true")
				});
		writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default("id", "treeModel1.allInternalNodeHeights"), true);
		writer.writeCloseTag(TreeModel.NODE_HEIGHTS);

		writer.writeCloseTag(TreeModel.TREE_MODEL);
	}

	/**
	 * Writes the substitution model to XML.
	 */
	public void writeSubstitutionModel(XMLWriter writer) {

		if (alignment.getDataType() == Nucleotides.INSTANCE) {

			// Jukes-Cantor model
			if (nucSubstitutionModel == JC) {
				writer.writeComment("The JC substitution model (Jukes & Cantor, 1969)");
				writer.writeOpenTag(
						dr.evomodel.substmodel.HKY.HKY_MODEL,
						new Attribute[]{new Attribute.Default("id", "jc1")}
				);
				writer.writeOpenTag(dr.evomodel.substmodel.HKY.FREQUENCIES);
				writer.writeOpenTag(
						FrequencyModel.FREQUENCY_MODEL,
						new Attribute[]{
								new Attribute.Default("dataType", alignment.getDataType().getDescription())
						}
				);
				writer.writeOpenTag(FrequencyModel.FREQUENCIES);
				writer.writeTag(
						ParameterParser.PARAMETER,
						new Attribute[]{
								new Attribute.Default("id", "jc1.frequencies"),
								new Attribute.Default("value", "0.25 0.25 0.25 0.25")
						},
						true
				);
				writer.writeCloseTag(FrequencyModel.FREQUENCIES);

				writer.writeCloseTag(FrequencyModel.FREQUENCY_MODEL);
				writer.writeCloseTag(dr.evomodel.substmodel.HKY.FREQUENCIES);

				writer.writeOpenTag(dr.evomodel.substmodel.HKY.KAPPA);
				writeParameter("hky1.kappa", 1.0, 1.0, 1.0, writer);
				writer.writeCloseTag(dr.evomodel.substmodel.HKY.KAPPA);
				writer.writeCloseTag(dr.evomodel.substmodel.HKY.HKY_MODEL);

			} else
				// Hasegawa Kishino and Yano 85 model
				if (nucSubstitutionModel == HKY) {

					for (int i = 1; i <= (unlinkedSubstitutionModel ? 3 : 1); i++) {
						writeHKYModel(i, writer);
					}
				} else
					// General time reversible model
					if (nucSubstitutionModel == GTR) {
						for (int i = 1; i <= (unlinkedSubstitutionModel ? 3 : 1); i++) {
							writeGTRModel(i, writer);
						}
					}
		} else {
			// Amino Acid model
			String aaModel = "";

			switch (aaSubstitutionModel) {
				case 0:
					aaModel = EmpiricalAminoAcidModel.BLOSUM_62;
					break;
				case 1:
					aaModel = EmpiricalAminoAcidModel.DAYHOFF;
					break;
				case 2:
					aaModel = EmpiricalAminoAcidModel.JTT;
					break;
				case 3:
					aaModel = EmpiricalAminoAcidModel.MT_REV_24;
					break;
				case 4:
					aaModel = EmpiricalAminoAcidModel.CP_REV_45;
					break;
				case 5:
					aaModel = EmpiricalAminoAcidModel.WAG;
					break;
			}

			writer.writeComment("The " + aaModel + " substitution model");
			writer.writeTag(
					EmpiricalAminoAcidModel.EMPIRICAL_AMINO_ACID_MODEL,
					new Attribute[]{new Attribute.Default("id", "aa1"),
							new Attribute.Default("type", aaModel)}, true
			);

		}
	}

	/**
	 * Write the HKY model XML block.
	 */
	public void writeHKYModel(int num, XMLWriter writer) {
		// Hasegawa Kishino and Yano 85 model
		writer.writeComment("The HKY substitution model (Hasegawa, Kishino & Yano, 1985)");
		writer.writeOpenTag(
				dr.evomodel.substmodel.HKY.HKY_MODEL,
				new Attribute[]{new Attribute.Default("id", "hky" + num)}
		);
		writer.writeOpenTag(dr.evomodel.substmodel.HKY.FREQUENCIES);
		writer.writeOpenTag(
				FrequencyModel.FREQUENCY_MODEL,
				new Attribute[]{
						new Attribute.Default("dataType", alignment.getDataType().getDescription())
				}
		);
		writer.writeTag("alignment", new Attribute[]{new Attribute.Default("idref", "alignment1")}, true);
		writer.writeOpenTag(FrequencyModel.FREQUENCIES);
		writer.writeTag(
				ParameterParser.PARAMETER,
				new Attribute[]{
						new Attribute.Default("id", "hky" + num + ".frequencies"),
						new Attribute.Default("dimension", "4")
				},
				true
		);
		writer.writeCloseTag(FrequencyModel.FREQUENCIES);
		writer.writeCloseTag(FrequencyModel.FREQUENCY_MODEL);
		writer.writeCloseTag(dr.evomodel.substmodel.HKY.FREQUENCIES);

		writer.writeOpenTag(dr.evomodel.substmodel.HKY.KAPPA);
		writeParameter("hky" + num + ".kappa", 1.0, 0.0, 100.0, writer);
		writer.writeCloseTag(dr.evomodel.substmodel.HKY.KAPPA);
		writer.writeCloseTag(dr.evomodel.substmodel.HKY.HKY_MODEL);
	}

	/**
	 * Write the GTR model XML block.
	 */
	public void writeGTRModel(int num, XMLWriter writer) {
		writer.writeComment("The general time reversible (GTR) substitution model");
		writer.writeOpenTag(
				dr.evomodel.substmodel.GTR.GTR_MODEL,
				new Attribute[]{new Attribute.Default("id", "gtr" + num)}
		);
		writer.writeOpenTag(dr.evomodel.substmodel.GTR.FREQUENCIES);
		writer.writeOpenTag(
				FrequencyModel.FREQUENCY_MODEL,
				new Attribute[]{
						new Attribute.Default("dataType", alignment.getDataType().getDescription())
				}
		);
		writer.writeTag("alignment", new Attribute[]{new Attribute.Default("idref", "alignment1")}, true);
		writer.writeOpenTag(FrequencyModel.FREQUENCIES);
		writer.writeTag(
				ParameterParser.PARAMETER,
				new Attribute[]{
						new Attribute.Default("id", "gtr" + num + ".frequencies"),
						new Attribute.Default("dimension", "4")
				},
				true
		);
		writer.writeCloseTag(FrequencyModel.FREQUENCIES);
		writer.writeCloseTag(FrequencyModel.FREQUENCY_MODEL);
		writer.writeCloseTag(dr.evomodel.substmodel.GTR.FREQUENCIES);

		writer.writeOpenTag(dr.evomodel.substmodel.GTR.A_TO_C);
		writeParameter("gtr" + num + ".ac", 1.0, 0.0, 100.0, writer);
		writer.writeCloseTag(dr.evomodel.substmodel.GTR.A_TO_C);

		writer.writeOpenTag(dr.evomodel.substmodel.GTR.A_TO_G);
		writeParameter("gtr" + num + ".ag", 1.0, 0.0, 100.0, writer);
		writer.writeCloseTag(dr.evomodel.substmodel.GTR.A_TO_G);

		writer.writeOpenTag(dr.evomodel.substmodel.GTR.A_TO_T);
		writeParameter("gtr" + num + ".at", 1.0, 0.0, 100.0, writer);
		writer.writeCloseTag(dr.evomodel.substmodel.GTR.A_TO_T);

		writer.writeOpenTag(dr.evomodel.substmodel.GTR.C_TO_G);
		writeParameter("gtr" + num + ".cg", 1.0, 0.0, 100.0, writer);
		writer.writeCloseTag(dr.evomodel.substmodel.GTR.C_TO_G);

		writer.writeOpenTag(dr.evomodel.substmodel.GTR.G_TO_T);
		writeParameter("gtr" + num + ".gt", 1.0, 0.0, 100.0, writer);
		writer.writeCloseTag(dr.evomodel.substmodel.GTR.G_TO_T);
		writer.writeCloseTag(dr.evomodel.substmodel.GTR.GTR_MODEL);
	}

	/**
	 * Write the site model XML block.
	 */
	public void writeSiteModel(XMLWriter writer) {
		if (alignment.getDataType() == Nucleotides.INSTANCE) {
			if (codonHetero) {
				for (int i = 1; i <= 3; i++) {
					writeNucSiteModel(i, writer);
				}
				writer.println();
				writer.writeOpenTag(CompoundParameter.COMPOUND_PARAMETER, new Attribute[]{new Attribute.Default("id", "allMus")});
				for (int i = 1; i <= 3; i++) {
					writer.writeTag(ParameterParser.PARAMETER,
							new Attribute[]{new Attribute.Default("idref", "siteModel" + i + ".mu")}, true);
				}
				writer.writeCloseTag(CompoundParameter.COMPOUND_PARAMETER);
			} else {
				writeNucSiteModel(1, writer);
			}
		} else {
			writeAASiteModel(writer);
		}
	}

	/**
	 * Write the site model XML block.
	 */
	public void writeNucSiteModel(int num, XMLWriter writer) {

		writer.writeComment("site model");
		writer.writeOpenTag(GammaSiteModel.SITE_MODEL, new Attribute[]{new Attribute.Default("id", "siteModel" + num)});


		writer.writeOpenTag(GammaSiteModel.SUBSTITUTION_MODEL);

		int n = 1;
		if (unlinkedSubstitutionModel) {
			n = num;
		}
		switch (nucSubstitutionModel) {
			// JC cannot be unlinked because it has no parameters 
			case JC:
				writer.writeTag(dr.evomodel.substmodel.HKY.HKY_MODEL, new Attribute.Default("idref", "jc1"), true);
				break;
			case HKY:
				writer.writeTag(dr.evomodel.substmodel.HKY.HKY_MODEL, new Attribute.Default("idref", "hky" + n), true);
				break;
			case GTR:
				writer.writeTag(dr.evomodel.substmodel.GTR.GTR_MODEL, new Attribute.Default("idref", "gtr" + n), true);
				break;
			default:
				throw new IllegalArgumentException("Unknown substitution model.");
		}
		writer.writeCloseTag(GammaSiteModel.SUBSTITUTION_MODEL);

		writer.writeOpenTag(GammaSiteModel.MUTATION_RATE);
		writeParameter("siteModel" + num + ".mu", meanMutationRate, 0.0, Double.MAX_VALUE, writer);
		writer.writeCloseTag(GammaSiteModel.MUTATION_RATE);

		if (gammaHetero) {
			writer.writeOpenTag(GammaSiteModel.GAMMA_SHAPE, new Attribute.Default(GammaSiteModel.GAMMA_CATEGORIES, "" + gammaCategories));
			if (num == 1 || unlinkedHeterogeneityModel) {
				writeParameter("siteModel" + num + ".shape", 1.0, 0.0, 1000.0, writer);
			} else {
				writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default("idref", "siteModel1.shape"), true);
			}
			writer.writeCloseTag(GammaSiteModel.GAMMA_SHAPE);
		}

		if (invarHetero) {
			writer.writeOpenTag(GammaSiteModel.PROPORTION_INVARIANT);
			if (num == 1 || unlinkedHeterogeneityModel) {
				writeParameter("siteModel" + num + ".pinv", 0.1, 0.0, 1.0, writer);
			} else {
				writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default("idref", "siteModel1.pinv"), true);
			}
			writer.writeCloseTag(GammaSiteModel.PROPORTION_INVARIANT);
		}

		writer.writeCloseTag(GammaSiteModel.SITE_MODEL);
	}

	/**
	 * Write the site model XML block.
	 */
	public void writeAASiteModel(XMLWriter writer) {

		writer.writeComment("site model");
		writer.writeOpenTag(GammaSiteModel.SITE_MODEL, new Attribute[]{
				new Attribute.Default("id", "siteModel1")});


		writer.writeOpenTag(GammaSiteModel.SUBSTITUTION_MODEL);
		writer.writeTag(EmpiricalAminoAcidModel.EMPIRICAL_AMINO_ACID_MODEL,
				new Attribute.Default("idref", "aa1"), true);
		writer.writeCloseTag(GammaSiteModel.SUBSTITUTION_MODEL);

		writer.writeOpenTag(GammaSiteModel.MUTATION_RATE);
		writeParameter("siteModel1.mu", meanMutationRate, 0.0, Double.MAX_VALUE, writer);
		writer.writeCloseTag(GammaSiteModel.MUTATION_RATE);

		if (gammaHetero) {
			writer.writeOpenTag(GammaSiteModel.GAMMA_SHAPE, new Attribute.Default(GammaSiteModel.GAMMA_CATEGORIES, "" + gammaCategories));
			writeParameter("siteModel1.shape", 1.0, 0.0, 1000.0, writer);
			writer.writeCloseTag(GammaSiteModel.GAMMA_SHAPE);
		}

		if (invarHetero) {
			writer.writeOpenTag(GammaSiteModel.PROPORTION_INVARIANT);
			writeParameter("siteModel1.pinv", 0.1, 0.0, 1.0, writer);
			writer.writeCloseTag(GammaSiteModel.PROPORTION_INVARIANT);
		}

		writer.writeCloseTag(GammaSiteModel.SITE_MODEL);
	}

	/**
	 * Write the coalescent likelihood
	 */
	public void writeCoalescentLikelihood(XMLWriter writer) {
		writer.writeOpenTag(
				VeryOldCoalescentLikelihood.COALESCENT_LIKELIHOOD,
				new Attribute[]{new Attribute.Default("id", "coalescent1")}
		);
		writer.writeOpenTag(VeryOldCoalescentLikelihood.MODEL);
		writeDemoModelRef(writer);
		writer.writeCloseTag(VeryOldCoalescentLikelihood.MODEL);
		writer.writeOpenTag(VeryOldCoalescentLikelihood.POPULATION_TREE);
		writer.writeTag(TreeModel.TREE_MODEL, new Attribute.Default("idref", "treeModel1"), true);
		writer.writeCloseTag(VeryOldCoalescentLikelihood.POPULATION_TREE);
		writer.writeCloseTag(VeryOldCoalescentLikelihood.COALESCENT_LIKELIHOOD);
	}

	/**
	 * Write the boolean likelihood
	 */
	public void writeBooleanLikelihood(XMLWriter writer) {
		writer.writeOpenTag(
				BooleanLikelihood.BOOLEAN_LIKELIHOOD,
				new Attribute[]{new Attribute.Default("id", "booleanLikelihood1")}
		);
		writer.writeOpenTag(
				TestStatistic.TEST_STATISTIC,
				new Attribute[]{
						new Attribute.Default("id", "test1"),
						new Attribute.Default("name", "test1")
				}
		);
		writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default("idref", "demo1.t50"), true);
		writer.writeOpenTag("lessThan");
		writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default("idref", "treeModel1.rootHeight"), true);
		writer.writeCloseTag("lessThan");
		writer.writeCloseTag(TestStatistic.TEST_STATISTIC);
		writer.writeCloseTag(BooleanLikelihood.BOOLEAN_LIKELIHOOD);
	}

	/**
	 * Write the tree likelihood XML block.
	 */
	public void writeTreeLikelihood(XMLWriter writer) {
		boolean nucs = alignment.getDataType() == Nucleotides.INSTANCE;
		if (nucs && codonHetero) {
			for (int i = 1; i <= 3; i++) {
				writeTreeLikelihood(i, writer);
			}
		} else {
			writeTreeLikelihood(1, writer);
		}
	}

	/**
	 * Write the tree likelihood XML block.
	 */
	public void writeTreeLikelihood(int num, XMLWriter writer) {

		writer.writeOpenTag(
				TreeLikelihood.TREE_LIKELIHOOD,
				new Attribute[]{new Attribute.Default("id", "treeLikelihood" + num)}
		);
		writer.writeTag(SitePatternsParser.PATTERNS, new Attribute[]{new Attribute.Default("idref", "patterns" + num)}, true);
		writer.writeTag(TreeModel.TREE_MODEL, new Attribute[]{new Attribute.Default("idref", "treeModel1")}, true);
		writer.writeTag(GammaSiteModel.SITE_MODEL, new Attribute[]{new Attribute.Default("idref", "siteModel" + num)}, true);

		writer.writeCloseTag(TreeLikelihood.TREE_LIKELIHOOD);
	}

	/**
	 * Write the operator schedule XML block.
	 */
	public void writeOperatorSchedule(XMLWriter writer) {

		writer.writeOpenTag(
				SimpleOperatorSchedule.OPERATOR_SCHEDULE,
				new Attribute[]{new Attribute.Default("id", "operators1")}
		);

		if (alignment != null) {

			if (coalescentModel != NONE) {
				writeScaleOperator(popSizeO, popSizeW, "demo1.popSize", writer);
				if (coalescentModel == EXPONENTIAL) {
					writeRandomWalkOperator(growthRateO, growthRateW, "demo1.growthRate", writer);
				} else if (coalescentModel == LOGISTIC) {
					writeScaleOperator(growthRateO, growthRateW, "demo1.doublingTime", writer);
					writeScaleOperator(logisticShapeO, logisticShapeW, "demo1.t50", writer);
				}
			}

			boolean nucs = alignment.getDataType() == Nucleotides.INSTANCE;

			if (nucs) {
				if (nucSubstitutionModel == HKY) {
					if (unlinkedSubstitutionModel) {
						writeScaleOperator(kappa1O, kappa1W, "hky1.kappa", writer);
						writeScaleOperator(kappa2O, kappa2W, "hky2.kappa", writer);
						writeScaleOperator(kappa3O, kappa3W, "hky3.kappa", writer);
					} else {
						writeScaleOperator(kappaO, kappaW, "hky1.kappa", writer);
					}
				} else if (nucSubstitutionModel == GTR) {
					if (unlinkedSubstitutionModel) {
						writeScaleOperator(gtrAC1O, gtrAC1W, "gtr1.ac", writer);
						writeScaleOperator(gtrAG1O, gtrAG1W, "gtr1.ag", writer);
						writeScaleOperator(gtrAT1O, gtrAT1W, "gtr1.at", writer);
						writeScaleOperator(gtrCG1O, gtrCG1W, "gtr1.cg", writer);
						writeScaleOperator(gtrGT1O, gtrGT1W, "gtr1.gt", writer);

						writeScaleOperator(gtrAC2O, gtrAC2W, "gtr2.ac", writer);
						writeScaleOperator(gtrAG2O, gtrAG2W, "gtr2.ag", writer);
						writeScaleOperator(gtrAT2O, gtrAT2W, "gtr2.at", writer);
						writeScaleOperator(gtrCG2O, gtrCG2W, "gtr2.cg", writer);
						writeScaleOperator(gtrGT2O, gtrGT2W, "gtr2.gt", writer);

						writeScaleOperator(gtrAC3O, gtrAC3W, "gtr3.ac", writer);
						writeScaleOperator(gtrAG3O, gtrAG3W, "gtr3.ag", writer);
						writeScaleOperator(gtrAT3O, gtrAT3W, "gtr3.at", writer);
						writeScaleOperator(gtrCG3O, gtrCG3W, "gtr3.cg", writer);
						writeScaleOperator(gtrGT3O, gtrGT3W, "gtr3.gt", writer);
					} else {
						writeScaleOperator(gtrACO, gtrACW, "gtr1.ac", writer);
						writeScaleOperator(gtrAGO, gtrAGW, "gtr1.ag", writer);
						writeScaleOperator(gtrATO, gtrATW, "gtr1.at", writer);
						writeScaleOperator(gtrCGO, gtrCGW, "gtr1.cg", writer);
						writeScaleOperator(gtrGTO, gtrGTW, "gtr1.gt", writer);
					}
				}
			}

			// if gamma do shape move
			if (gammaHetero) {
				if (nucs && unlinkedHeterogeneityModel) {
					writeScaleOperator(alpha1O, alpha1W, "siteModel1.shape", writer);
					writeScaleOperator(alpha2O, alpha2W, "siteModel2.shape", writer);
					writeScaleOperator(alpha3O, alpha3W, "siteModel3.shape", writer);
				} else {
					writeScaleOperator(alphaO, alphaW, "siteModel1.shape", writer);
				}
			}
			// if pinv do pinv move
			if (invarHetero) {
				if (nucs && unlinkedHeterogeneityModel) {
					writeScaleOperator(pinv1O, pinv1W, "siteModel1.pinv", writer);
					writeScaleOperator(pinv2O, pinv2W, "siteModel2.pinv", writer);
					writeScaleOperator(pinv3O, pinv3W, "siteModel3.pinv", writer);
				} else {
					writeScaleOperator(pinvO, pinvW, "siteModel1.pinv", writer);
				}
			}

			// if using dates do mutation rate move and up/down move
			if (usingDates) {

				if (nucs && codonHetero) {
					writeScaleOperator(mu1O, mu1W, "siteModel1.mu", writer);
					writeScaleOperator(mu2O, mu2W, "siteModel2.mu", writer);
					writeScaleOperator(mu3O, mu3W, "siteModel3.mu", writer);
				} else {
					writeScaleOperator(muO, muW, "siteModel1.mu", writer);
				}

				writer.writeOpenTag(UpDownOperator.UP_DOWN_OPERATOR,
						new Attribute[]{
								new Attribute.Default("weight", Integer.toString((int) muHeightsW)),
								new Attribute.Default("scaleFactor", Double.toString(muHeightsO)),
								new Attribute.Default("adapt", "false")}
				);
				writer.writeOpenTag(UpDownOperator.UP);
/*				writer.writeOpenTag(LogScaleUpDownOperator.UP_DOWN_OPERATOR, 
					new Attribute[] { 
						new Attribute.Default("weight", Integer.toString((int)muHeightsW)),
						new Attribute.Default("ordersOfMagnitude", Double.toString(muHeightsO)),
						new Attribute.Default("adapt", "false")}
				);
				writer.writeOpenTag(LogScaleUpDownOperator.UP);
*/
				if (codonHetero) {
					writer.writeTag(ParameterParser.PARAMETER, new Attribute[]{new Attribute.Default("idref", "allMus")}, true);
				} else {
					writer.writeTag(ParameterParser.PARAMETER, new Attribute[]{new Attribute.Default("idref", "siteModel1.mu")}, true);
				}
				writer.writeCloseTag(UpDownOperator.UP);

				writer.writeOpenTag(UpDownOperator.DOWN);
				writer.writeTag(ParameterParser.PARAMETER, new Attribute[]{new Attribute.Default("idref", "treeModel1.allInternalNodeHeights")}, true);
				writer.writeCloseTag(UpDownOperator.DOWN);

				writer.writeCloseTag(UpDownOperator.UP_DOWN_OPERATOR);
			} else if (nucs && codonHetero) {

				writer.writeOpenTag(CenteredScaleOperator.CENTERED_SCALE,
						new Attribute[]{
								new Attribute.Default("weight", Integer.toString((int) musCenteredW)),
								new Attribute.Default(CenteredScaleOperator.SCALE_FACTOR, Double.toString(musCenteredO)),
								new Attribute.Default("adapt", "false")}
				);
				writer.writeTag(ParameterParser.PARAMETER, new Attribute[]{new Attribute.Default("idref", "allMus")}, true);
				writer.writeCloseTag(CenteredScaleOperator.CENTERED_SCALE);

				writer.writeOpenTag(DeltaExchangeOperator.DELTA_EXCHANGE,
						new Attribute[]{
								new Attribute.Default("weight", Integer.toString((int) musDeltaW)),
								new Attribute.Default(DeltaExchangeOperator.DELTA, Double.toString(musDeltaO)),
								new Attribute.Default("adapt", "false")}
				);
				writer.writeTag(ParameterParser.PARAMETER, new Attribute[]{new Attribute.Default("idref", "allMus")}, true);
				writer.writeCloseTag(DeltaExchangeOperator.DELTA_EXCHANGE);

			}
		}

		writeScaleOperator(rootO, rootW, "treeModel1.rootHeight", writer);

		writer.writeOpenTag("uniformOperator", new Attribute.Default("weight", Integer.toString((int) heightsW)));
		writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default("idref", "treeModel1.internalNodeHeights"), true);
		writer.writeCloseTag("uniformOperator");

		// if not a fixed tree then sample tree space
		if (!fixedTree) {
			writer.writeOpenTag(ExchangeOperator.NARROW_EXCHANGE, new Attribute[]{
					new Attribute.Default("weight", Integer.toString((int) wideExchangeW))});
			writer.writeTag(TreeModel.TREE_MODEL, new Attribute[]{
					new Attribute.Default("idref", "treeModel1")}, true);
			writer.writeCloseTag(ExchangeOperator.NARROW_EXCHANGE);

			writer.writeOpenTag(ExchangeOperator.WIDE_EXCHANGE, new Attribute[]{
					new Attribute.Default("weight", Integer.toString((int) narrowExchangeW))});
			writer.writeTag(TreeModel.TREE_MODEL, new Attribute[]{
					new Attribute.Default("idref", "treeModel1")}, true);
			writer.writeCloseTag(ExchangeOperator.WIDE_EXCHANGE);

			writer.writeOpenTag(WilsonBalding.WILSON_BALDING, new Attribute[]{
					new Attribute.Default("weight", Integer.toString((int) wilsonBaldingW))});
			if (coalescentModel != NONE) {
				writeDemoModelRef(writer);
			}
			writer.writeTag(TreeModel.TREE_MODEL, new Attribute[]{
					new Attribute.Default("idref", "treeModel1")}, true);
			writer.writeCloseTag(WilsonBalding.WILSON_BALDING);

			writer.writeOpenTag(SubtreeSlideOperator.SUBTREE_SLIDE,
					new Attribute[]{
							new Attribute.Default("weight", Integer.toString((int) subtreeSlideW)),
							new Attribute.Default("gaussian", "true"),
							new Attribute.Default("size", Double.toString(subtreeSlideO))
					}
			);
			writer.writeTag(TreeModel.TREE_MODEL, new Attribute[]{new Attribute.Default("idref", "treeModel1")}, true);
			writer.writeCloseTag(SubtreeSlideOperator.SUBTREE_SLIDE);
		}

		writer.writeCloseTag(SimpleOperatorSchedule.OPERATOR_SCHEDULE);
	}

	private void writeScaleOperator(double scaleFactor, double weight, String parameterName, XMLWriter writer) {
		writer.writeOpenTag(
				ScaleOperator.SCALE_OPERATOR,
				new Attribute[]{
						new Attribute.Default("scaleFactor", scaleFactor + ""),
						new Attribute.Default("weight", Integer.toString((int) weight)),
						new Attribute.Default("adapt", "false"),
				});
		writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default("idref", parameterName), true);
		writer.writeCloseTag(ScaleOperator.SCALE_OPERATOR);
	}

	private void writeRandomWalkOperator(double windowSize, double weight, String parameterName, XMLWriter writer) {
		writer.writeOpenTag(
				"randomWalkOperator",
				new Attribute[]{
						new Attribute.Default("windowSize", windowSize + ""),
						new Attribute.Default("weight", Integer.toString((int) weight)),
						new Attribute.Default("adapt", "false"),
				});
		writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default("idref", parameterName), true);
		writer.writeCloseTag("randomWalkOperator");
	}

	/**
	 * Write the timer report block.
	 */
	public void writeTimerReport(XMLWriter writer) {
		writer.writeOpenTag("report");
		writer.writeOpenTag("property", new Attribute.Default("name", "timer"));
		writer.writeTag("object", new Attribute.Default("idref", "mcmc1"), true);
		writer.writeCloseTag("property");
		writer.writeCloseTag("report");
	}

	/**
	 * Write the trace analysis block.
	 */
	public void writeTraceAnalysis(XMLWriter writer) {
		writer.writeTag(
				"traceAnalysis",
				new Attribute[]{
						new Attribute.Default("fileName", fileNameStem + ".log")
				},
				true
		);
	}

	/**
	 * Write the MCMC block.
	 */
	public void writeMCMC(XMLWriter writer) {
		writer.writeOpenTag(
				"mcmc",
				new Attribute[]{
						new Attribute.Default("id", "mcmc1"),
						new Attribute.Default("chainLength", chainLength + "")
				});

		// write likelihood block	
		writer.writeOpenTag(CompoundLikelihood.COMPOUND_LIKELIHOOD, new Attribute.Default("id", "likelihood1"));

		if (coalescentModel != NONE) {
			writer.writeTag(VeryOldCoalescentLikelihood.COALESCENT_LIKELIHOOD, new Attribute.Default("idref", "coalescent1"), true);
		}
		boolean nucs = alignment.getDataType() == Nucleotides.INSTANCE;
		if (nucs && codonHetero) {
			for (int i = 1; i <= 3; i++) {
				writer.writeTag(TreeLikelihood.TREE_LIKELIHOOD, new Attribute.Default("idref", "treeLikelihood" + i), true);
			}
		} else writer.writeTag(TreeLikelihood.TREE_LIKELIHOOD, new Attribute.Default("idref", "treeLikelihood1"), true);

		if (coalescentModel == LOGISTIC) {
			writer.writeTag(BooleanLikelihood.BOOLEAN_LIKELIHOOD, new Attribute.Default("idref", "booleanLikelihood1"), true);
		}

		writer.writeCloseTag(CompoundLikelihood.COMPOUND_LIKELIHOOD);

		writer.writeTag(SimpleOperatorSchedule.OPERATOR_SCHEDULE, new Attribute.Default("idref", "operators1"), true);

		// write log to screen
		writer.writeOpenTag(MCLogger.LOG,
				new Attribute[]{
						new Attribute.Default("id", "screenLog"),
						new Attribute.Default(MCLogger.LOG_EVERY, echoEvery + "")
				});
		writer.writeOpenTag(Columns.COLUMN,
				new Attribute[]{
						new Attribute.Default(Columns.DECIMAL_PLACES, "4")
				}
		);
		writeLog(false, writer);
		writer.writeCloseTag(Columns.COLUMN);
		writer.writeCloseTag(MCLogger.LOG);

		// write log to file
		if (logFileName == null) {
			logFileName = fileNameStem + ".log";
		}
		writer.writeOpenTag(MCLogger.LOG,
				new Attribute[]{
						new Attribute.Default("id", "fileLog"),
						new Attribute.Default(MCLogger.LOG_EVERY, logEvery + ""),
						new Attribute.Default(MCLogger.FILE_NAME, logFileName)
				});
		writeLog(true, writer);
		writer.writeCloseTag(MCLogger.LOG);

		// write tree log to file
		if (treeFileName == null) {
			treeFileName = fileNameStem + ".trees";
		}
		writer.writeOpenTag(TreeLogger.LOG_TREE,
				new Attribute[]{
						new Attribute.Default("id", "treeFileLog"),
						new Attribute.Default(TreeLogger.LOG_EVERY, logEvery + ""),
						new Attribute.Default(TreeLogger.NEXUS_FORMAT, "true"),
						new Attribute.Default(TreeLogger.FILE_NAME, treeFileName)
				});
		writer.writeTag(TreeModel.TREE_MODEL, new Attribute.Default("idref", "treeModel1"), true);
		writer.writeCloseTag(TreeLogger.LOG_TREE);

		writer.writeCloseTag("mcmc");
	}

	/**
	 * Write the log
	 */
	private void writeLog(boolean verbose, XMLWriter writer) {
		if (alignment != null) {
			writer.writeTag(CompoundLikelihood.COMPOUND_LIKELIHOOD, new Attribute.Default("idref", "likelihood1"), true);
		}
		if (verbose && alignment != null) {
			boolean nucs = alignment.getDataType() == Nucleotides.INSTANCE;
			if (nucs && codonHetero) {
				for (int i = 1; i <= 3; i++) {
					writer.writeTag(ParameterParser.PARAMETER,
							new Attribute.Default("idref", "siteModel" + i + ".mu"), true);
				}
			} else if (usingDates) {
				writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default("idref", "siteModel1.mu"), true);
			}
			if (gammaHetero) {
				if (nucs && codonHetero && unlinkedHeterogeneityModel) {
					for (int i = 1; i <= 3; i++) {
						writer.writeTag(ParameterParser.PARAMETER,
								new Attribute.Default("idref", "siteModel" + i + ".shape"), true);
					}
				} else {
					writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default("idref", "siteModel1.shape"), true);
				}
			}
			if (invarHetero) {
				if (nucs && codonHetero && unlinkedHeterogeneityModel) {
					for (int i = 1; i <= 3; i++) {
						writer.writeTag(ParameterParser.PARAMETER,
								new Attribute.Default("idref", "siteModel" + i + ".pinv"), true);
					}
				} else {
					writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default("idref", "siteModel1.pinv"), true);
				}
			}

			if (coalescentModel != NONE) {

				writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default("idref", "demo1.popSize"), true);

				if (coalescentModel == EXPONENTIAL) {
					writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default("idref", "demo1.growthRate"), true);
				} else if (coalescentModel == LOGISTIC) {
					writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default("idref", "demo1.doublingTime"), true);
					writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default("idref", "demo1.t50"), true);
				}
			}

			if (nucs) {
				if (nucSubstitutionModel == HKY) {
					if (unlinkedSubstitutionModel) {
						writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default("idref", "hky1.kappa"), true);
						writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default("idref", "hky2.kappa"), true);
						writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default("idref", "hky3.kappa"), true);
					} else {
						writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default("idref", "hky1.kappa"), true);
					}
				} else if (nucSubstitutionModel == GTR) {
					if (unlinkedSubstitutionModel) {
						writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default("idref", "gtr1.ac"), true);
						writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default("idref", "gtr1.ag"), true);
						writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default("idref", "gtr1.at"), true);
						writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default("idref", "gtr1.cg"), true);
						writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default("idref", "gtr1.gt"), true);

						writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default("idref", "gtr2.ac"), true);
						writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default("idref", "gtr2.ag"), true);
						writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default("idref", "gtr2.at"), true);
						writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default("idref", "gtr2.cg"), true);
						writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default("idref", "gtr2.gt"), true);

						writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default("idref", "gtr3.ac"), true);
						writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default("idref", "gtr3.ag"), true);
						writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default("idref", "gtr3.at"), true);
						writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default("idref", "gtr3.cg"), true);
						writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default("idref", "gtr3.gt"), true);
					} else {
						writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default("idref", "gtr1.ac"), true);
						writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default("idref", "gtr1.ag"), true);
						writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default("idref", "gtr1.at"), true);
						writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default("idref", "gtr1.cg"), true);
						writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default("idref", "gtr1.gt"), true);
					}
				}
			}
		}

		writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default("idref", "treeModel1.rootHeight"), true);

		if (alignment != null) {
			boolean nucs = alignment.getDataType() == Nucleotides.INSTANCE;
			if (nucs && codonHetero) {
				for (int i = 1; i <= 3; i++) {
					writer.writeTag(TreeLikelihood.TREE_LIKELIHOOD, new Attribute.Default("idref", "treeLikelihood" + i), true);
				}
			} else
				writer.writeTag(TreeLikelihood.TREE_LIKELIHOOD, new Attribute.Default("idref", "treeLikelihood1"), true);
		}
		if (coalescentModel != NONE) {
			writer.writeTag(VeryOldCoalescentLikelihood.COALESCENT_LIKELIHOOD, new Attribute.Default("idref", "coalescent1"), true);
		}
	}

	/**
	 * write a parameter
	 */
	public void writeParameter(String id, double value, double lower, double upper, XMLWriter writer) {
		ArrayList attributes = new ArrayList();
		attributes.add(new Attribute.Default("id", id));
		attributes.add(new Attribute.Default("value", value + ""));
		if (lower > -Double.MAX_VALUE) {
			attributes.add(new Attribute.Default("lower", lower + ""));
		}
		if (upper < Double.MAX_VALUE) {
			attributes.add(new Attribute.Default("upper", upper + ""));
		}

		Attribute[] attrArray = new Attribute[attributes.size()];
		for (int i = 0; i < attrArray.length; i++) {
			attrArray[i] = (Attribute) attributes.get(i);
		}

		writer.writeTag(ParameterParser.PARAMETER, attrArray, true);
	}

	/**
	 * Generate XML for the starting tree
	 */
	public void writeStartingTree(XMLWriter writer) {
		if (userTree) {
			writeUserTree(tree, writer);
		} else {
			// generate a coalescent tree
			writer.writeComment("Generate a random starting tree under the coalescent process");
			writer.writeOpenTag(
					CoalescentSimulator.COALESCENT_TREE,
					new Attribute[]{new Attribute.Default("id", "startingTree")}
			);
			writer.writeTag(TaxaParser.TAXA, new Attribute[]{new Attribute.Default("idref", "taxa1")}, true);
			writeInitialDemoModelRef(writer);
			writer.writeCloseTag(CoalescentSimulator.COALESCENT_TREE);
		}
	}

	public void writeInitialDemoModelRef(XMLWriter writer) {
		if (coalescentModel == CONSTANT) {
			writer.writeTag(ConstantPopulationModel.CONSTANT_POPULATION_MODEL, new Attribute[]{new Attribute.Default("idref", "demo1")}, true);
		} else if (coalescentModel == EXPONENTIAL) {
			writer.writeTag(ExponentialGrowthModel.EXPONENTIAL_GROWTH_MODEL, new Attribute[]{new Attribute.Default("idref", "demo1")}, true);
		} else if (coalescentModel == LOGISTIC || coalescentModel == NONE) {
			writer.writeTag(ConstantPopulationModel.CONSTANT_POPULATION_MODEL, new Attribute[]{new Attribute.Default("idref", "initialDemo")}, true);
		}
	}

	public void writeDemoModelRef(XMLWriter writer) {
		writer.writeTag(ConstantPopulationModel.CONSTANT_POPULATION_MODEL, new Attribute[]{new Attribute.Default("idref", "demo1")}, true);
	}

	/**
	 * Generate XML for the user tree
	 */
	private void writeUserTree(Tree tree, XMLWriter writer) {

		writer.writeComment("The starting tree.");
		writer.writeOpenTag(
				"tree",
				new Attribute[]{
						new Attribute.Default("height", "startingTree"),
						new Attribute.Default("usingDates", "" + usingDates)
				}
		);
		writeNode(tree, tree.getRoot(), writer);
		writer.writeCloseTag("tree");
	}

	/**
	 * Generate XML for the node of a user tree.
	 */
	private void writeNode(Tree tree, NodeRef node, XMLWriter writer) {

		writer.writeOpenTag(
				"node",
				new Attribute[]{new Attribute.Default("height", "" + tree.getNodeHeight(node))}
		);

		if (tree.getChildCount(node) == 0) {
			writer.writeTag("taxon", new Attribute[]{new Attribute.Default("idref", tree.getNodeTaxon(node).getId())}, true);
		}
		for (int i = 0; i < tree.getChildCount(node); i++) {
			writeNode(tree, tree.getChild(node, i), writer);
		}
		writer.writeCloseTag("node");
	}
}

