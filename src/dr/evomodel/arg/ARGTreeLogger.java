/*
 * ARGTreeLogger.java
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

package dr.evomodel.arg;

import dr.evolution.tree.Tree;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.coalescent.structure.ColourSamplerModel;
import dr.inference.loggers.LogFormatter;
import dr.inference.loggers.MLLogger;
import dr.inference.loggers.TabDelimitedFormatter;
import dr.inference.model.Likelihood;
import dr.xml.*;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintWriter;

/**
 * A logger that logs tree and clade frequencies from a given partition in an ARG
 *
 * @author Marc Suchard
 * @version $Id: ARGTreeLogger.java,v 1.25 2006/09/05 13:29:34 rambaut Exp $
 */
public class ARGTreeLogger extends OldTreeLogger {

	public static final String LOG_ARG = "logArgTree";
	public static final String PARTITION = "partition";


	private int partition;

	/**
	 * Constructor
	 */
	public ARGTreeLogger(Tree tree, int partition, BranchRateModel branchRateModel, String rateLabel,
	                     ColourSamplerModel colourSamplerModel, String colouringLabel,
	                     Likelihood likelihood, String likelihoodLabel,
	                     LogFormatter formatter, int logEvery, boolean nexusFormat, boolean substitutions) {

		super(tree, branchRateModel, rateLabel, colourSamplerModel, colouringLabel, likelihood, likelihoodLabel, formatter, logEvery, nexusFormat, substitutions);
		this.partition = partition;

	}


	@Override
	protected String additionalInfo() {
		return " [&PARTITION=" + partition + "]"
				+ " [&YULE=" + getLogYuleProbabilityString() + "]"
				+ " [&NUM_REC=" + getNumberOfReassortments() + "]"
				;
	}

	private int getNumberOfReassortments() {
		ARGModel arg = (ARGModel) getTree();
		return arg.getReassortmentNodeCount();
	}

	private String getLogYuleProbabilityString() {
		ARGTree tree = new ARGTree((ARGModel) getTree(), partition);
//		BetaSplittingModel betaModel = new BetaSplittingModel(
//				new Parameter.Default(1.0), tree);
//		betaModel.setBeta(0.0);
//		double otherLP = 0;
		double logProbability = 0;
		for (int i = 0, n = tree.getNodeCount(); i < n; i++) {
//			System.err.println(n);
			ARGModel.Node node = (ARGModel.Node) tree.getNode(i);
			int count = node.getDescendentTipCount();
//			System.err.println(count);
			if (count > 2)
				logProbability -= 2 * Math.log(count - 1);
//			otherLP += betaModel.logNodeProbability(tree,node);
		}
//		System.err.println("me : "+logProbability);
//		System.err.println("old: "+otherLP);
		return String.format("%5.4f", logProbability);
	}

	@Override
	protected Tree getPrintTree() {
		return new ARGTree((ARGModel) getTree(), partition);
	}

	@Override
	protected boolean useTaxonLabels() {
		return true;
	}

	public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

		public String getParserName() {
			return LOG_ARG;
		}

		/**
		 * @return an object based on the XML element it was passed.
		 */
		public Object parseXMLObject(XMLObject xo) throws XMLParseException {

			ARGModel tree = (ARGModel) xo.getChild(ARGModel.class);

			String fileName = null;
			String title = null;
			boolean nexusFormat = false;

			String colouringLabel = "demes";
			String rateLabel = "rate";
			String likelihoodLabel = "lnP";

			if (xo.hasAttribute(TITLE)) {
				title = xo.getStringAttribute(TITLE);
			}

			if (xo.hasAttribute(FILE_NAME)) {
				fileName = xo.getStringAttribute(FILE_NAME);
			}

			if (xo.hasAttribute(NEXUS_FORMAT)) {
				nexusFormat = xo.getBooleanAttribute(NEXUS_FORMAT);
			}

			boolean substitutions = false;
			if (xo.hasAttribute(BRANCH_LENGTHS)) {
				substitutions = xo.getStringAttribute(BRANCH_LENGTHS).equals(SUBSTITUTIONS);
			}

			BranchRateModel branchRateModel = (BranchRateModel) xo.getChild(BranchRateModel.class);

			ColourSamplerModel colourSamplerModel = (ColourSamplerModel) xo.getChild(ColourSamplerModel.class);

			Likelihood likelihood = (Likelihood) xo.getChild(Likelihood.class);

			// logEvery of zero only displays at the end
			int logEvery = 1;

			if (xo.hasAttribute(LOG_EVERY)) {
				logEvery = xo.getIntegerAttribute(LOG_EVERY);
			}

			int partition = 0;

			if (xo.hasAttribute(PARTITION)) {
				partition = xo.getIntegerAttribute(PARTITION);
			}

			if (partition > tree.getMaxPartitionNumber())
				throw new XMLParseException("ARGModel does not contain a partition #" + partition);

			PrintWriter pw = null;

			if (fileName != null) {

				try {
					File file = new File(fileName);
					String name = file.getName();
					String parent = file.getParent();

					if (!file.isAbsolute()) {
						parent = System.getProperty("user.dir");
					}

//					System.out.println("Writing log file to "+parent+System.getProperty("path.separator")+name);
					pw = new PrintWriter(new FileOutputStream(new File(parent, name)));
				} catch (FileNotFoundException fnfe) {
					throw new XMLParseException("File '" + fileName + "' can not be opened for " + getParserName() + " element.");
				}
			} else {
				pw = new PrintWriter(System.out);
			}

			LogFormatter formatter = new TabDelimitedFormatter(pw);

			ARGTreeLogger logger = new ARGTreeLogger(tree, partition, branchRateModel, rateLabel,
					colourSamplerModel, colouringLabel, likelihood, likelihoodLabel,
					formatter, logEvery, nexusFormat, substitutions);

			if (title != null) {
				logger.setTitle(title);
			}

			return logger;
		}

		//************************************************************************
		// AbstractXMLObjectParser implementation
		//************************************************************************
		public XMLSyntaxRule[] getSyntaxRules() {
			return rules;
		}

		private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
				AttributeRule.newIntegerRule(LOG_EVERY),
				AttributeRule.newIntegerRule(PARTITION),
				new StringAttributeRule(FILE_NAME,
						"The name of the file to send log output to. " +
								"If no file name is specified then log is sent to standard output", true),
				new StringAttributeRule(TITLE, "The title of the log", true),
				AttributeRule.newBooleanRule(NEXUS_FORMAT, true,
						"Whether to use the NEXUS format for the tree log"),
				new StringAttributeRule(BRANCH_LENGTHS, "What units should the branch lengths be in", new String[]{TIME, SUBSTITUTIONS}, true),
				new ElementRule(ARGModel.class, "The ARG which is to be logged"),
				new ElementRule(BranchRateModel.class, true),
				new ElementRule(ColourSamplerModel.class, true),
				new ElementRule(Likelihood.class, true)
		};

		public String getParserDescription() {
			return "Logs a tree to a file";
		}

		public String getExample() {
			return
					"<!-- The " + getParserName() + " element takes an argTreeModel to be logged -->\n" +
							"<" + getParserName() + " " + LOG_EVERY + "=\"100\" " + FILE_NAME + "=\"log.trees\" " + NEXUS_FORMAT + "=\"true\">\n" +
							"	<treeModel idref=\"treeModel1\"/>\n" +
							"</" + getParserName() + ">\n";
		}

		public Class getReturnType() {
			return MLLogger.class;
		}
	};

}