/*
 * TreeLogger.java
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

package dr.evomodel.tree;

import dr.evolution.colouring.*;
import dr.evolution.tree.*;
import dr.inference.loggers.*;
import dr.inference.model.Likelihood;
import dr.xml.*;

import java.io.*;
import java.util.ArrayList;
import java.util.List;

/**
 * A logger that logs tree and clade frequencies.
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 *
 * @version $Id: TreeLogger.java,v 1.25 2006/09/05 13:29:34 rambaut Exp $
 */
public class TreeLogger extends MCLogger {

	public static final String LOG_TREE = "logTree";
	public static final String NEXUS_FORMAT = "nexusFormat";
	public static final String USING_RATES = "usingRates";
	public static final String BRANCH_LENGTHS = "branchLengths";
	public static final String TIME = "time";
	public static final String SUBSTITUTIONS = "substitutions";

	private Tree tree;
	private BranchRateController branchRateProvider = null;

	private TreeAttributeProvider[] treeAttributeProviders;
	private NodeAttributeProvider[] nodeAttributeProviders;
	private BranchAttributeProvider[] branchAttributeProviders;

	private boolean nexusFormat = false;
	public boolean usingRates = false;
	public boolean substitutions = false;

	/**
	 * Constructor
	 */
	public TreeLogger(Tree tree, BranchRateController branchRateProvider,
	                  TreeAttributeProvider[] treeAttributeProviders,
	                  NodeAttributeProvider[] nodeAttributeProviders,
	                  BranchAttributeProvider[] branchAttributeProviders,
	                  LogFormatter formatter, int logEvery, boolean nexusFormat) {

		super(formatter, logEvery);

		this.nexusFormat = nexusFormat;

		this.branchRateProvider = branchRateProvider;

		this.treeAttributeProviders = treeAttributeProviders;
		this.nodeAttributeProviders = nodeAttributeProviders;
		this.branchAttributeProviders = branchAttributeProviders;

		if (this.branchRateProvider != null) {
			this.substitutions = true;
		}
		this.tree = tree;
	}

	public void startLogging() {

		if (nexusFormat) {
			int taxonCount = tree.getTaxonCount();
			logLine("#NEXUS");
			logLine("");
			logLine("Begin taxa;");
			logLine("\tDimensions ntax=" + taxonCount + ";");
			logLine("\tTaxlabels");
			for (int i = 0; i < taxonCount; i++) {
				logLine("\t\t" + tree.getTaxon(i).getId());
			}
			logLine("\t\t;");
			logLine("End;");
			logLine("");
			logLine("Begin trees;");

			// This is needed if the trees use numerical taxon labels
			logLine("\tTranslate");
			for (int i = 0; i < taxonCount; i++) {
				int k = i + 1;
				if (k < taxonCount) {
					logLine("\t\t" + k + " " + tree.getTaxonId(i) + ",");
				} else {
					logLine("\t\t" + k + " " + tree.getTaxonId(i));
				}
			}
			logLine("\t\t;");
		}
	}

	public void log(int state) {

		if (logEvery <= 0 || ((state % logEvery) == 0)) {
			StringBuffer buffer = new StringBuffer("tree STATE_");
			buffer.append(state);
			if (treeAttributeProviders != null) {
				boolean hasAttribute = false;
				for (TreeAttributeProvider tap : treeAttributeProviders) {
					if (!hasAttribute) {
						buffer.append(" [&");
						hasAttribute = true;
					} else {
						buffer.append(",");
					}
					buffer.append(tap.getTreeAttributeLabel());
					buffer.append("=");
					buffer.append(tap.getAttributeForTree(tree));

				}
				if (hasAttribute) {
					buffer.append("]");
				}
			}

			buffer.append(" = [&R] ");

			if (substitutions) {
				Tree.Utils.newick(tree, tree.getRoot(), false, Tree.Utils.LENGTHS_AS_SUBSTITUTIONS,
						branchRateProvider, null, null, buffer);
			} else {
//				if (treeColouringProvider != null) {
//					TreeColouring colouring = treeColouringProvider.getTreeColouring(tree);
//
//					Tree.Utils.newick(tree, tree.getRoot(), false, Tree.Utils.LENGTHS_AS_TIME,
//							branchRateProvider, rateLabel, colouring, colouringLabel, buffer);
//				} else {
//					Tree.Utils.newick(tree, tree.getRoot(), false, Tree.Utils.LENGTHS_AS_TIME,
//							branchRateProvider, rateLabel, null, null, buffer);
//				}
				Tree.Utils.newick(tree, tree.getRoot(), false, Tree.Utils.LENGTHS_AS_TIME,
						null, nodeAttributeProviders, branchAttributeProviders, buffer);
			}

			buffer.append(";");
			logLine(buffer.toString());
		}
	}

	public void stopLogging() {

		logLine("End;");
		super.stopLogging();
	}

	public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

		public String getParserName() { return LOG_TREE; }

		/**
		 * @return an object based on the XML element it was passed.
		 */
		public Object parseXMLObject(XMLObject xo) throws XMLParseException {

			Tree tree = (Tree)xo.getChild(Tree.class);

			String fileName = null;
			String title = null;
			boolean nexusFormat = false;

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

			List<TreeAttributeProvider> taps = new ArrayList<TreeAttributeProvider>();
			List<NodeAttributeProvider> naps = new ArrayList<NodeAttributeProvider>();
			List<BranchAttributeProvider> baps = new ArrayList<BranchAttributeProvider>();

			for (int i = 0; i < xo.getChildCount(); i++) {
				Object cxo = xo.getChild(i);
				if (cxo instanceof TreeColouringProvider) {
					final TreeColouringProvider colouringProvider = (TreeColouringProvider)cxo;
					baps.add(new BranchAttributeProvider() {

						public String getBranchAttributeLabel() {
							return "deme";
						}

						public String getAttributeForBranch(Tree tree, NodeRef node) {
							TreeColouring colouring = colouringProvider.getTreeColouring(tree);
							BranchColouring bcol = colouring.getBranchColouring(node);
							StringBuilder buffer = new StringBuilder();
							if (bcol != null) {
								buffer.append("{");
								buffer.append(bcol.getChildColour());
								for (int i = 1; i <= bcol.getNumEvents(); i++) {
									buffer.append(",");
									buffer.append(bcol.getBackwardTime(i));
									buffer.append(",");
									buffer.append(bcol.getBackwardColourAbove(i));
								}
								buffer.append("}");
							}
							return buffer.toString();
						}
					});

				} else if (cxo instanceof Likelihood) {
						final Likelihood likelihood = (Likelihood)cxo;
						taps.add(new TreeAttributeProvider() {

							public String getTreeAttributeLabel() {
								return "lnP";
							}

							public String getAttributeForTree(Tree tree) {
								return Double.toString(likelihood.getLogLikelihood());
							}
						});

				} else {
					if (cxo instanceof TreeAttributeProvider) {
						taps.add((TreeAttributeProvider)cxo);
					}
					if (cxo instanceof NodeAttributeProvider) {
						naps.add((NodeAttributeProvider)cxo);
					}
					if (cxo instanceof BranchAttributeProvider) {
						baps.add((BranchAttributeProvider)cxo);
					}
				}
			}

			BranchRateController branchRateProvider = (BranchRateController)xo.getChild(BranchRateController.class);
			if (substitutions && branchRateProvider == null) {
				throw new XMLParseException("To log trees in units of substitutions a BranchRateModel must be provided");
			}

			// logEvery of zero only displays at the end
			int logEvery = 1;

			if (xo.hasAttribute(LOG_EVERY)) {
				logEvery = xo.getIntegerAttribute(LOG_EVERY);
			}

			PrintWriter pw;

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

			TreeAttributeProvider[] treeAttributeProviders = new TreeAttributeProvider[taps.size()];
			taps.toArray(treeAttributeProviders);
			NodeAttributeProvider[] nodeAttributeProviders = new NodeAttributeProvider[naps.size()];
			naps.toArray(nodeAttributeProviders);
			BranchAttributeProvider[] branchAttributeProviders = new BranchAttributeProvider[baps.size()];
			baps.toArray(branchAttributeProviders);

			TreeLogger logger = new TreeLogger(tree,
					branchRateProvider,
					treeAttributeProviders, nodeAttributeProviders, branchAttributeProviders,
					formatter, logEvery, nexusFormat);

			if (title != null) {
				logger.setTitle(title);
			}

			return logger;
		}

		//************************************************************************
		// AbstractXMLObjectParser implementation
		//************************************************************************
		public XMLSyntaxRule[] getSyntaxRules() { return rules; }

		private XMLSyntaxRule[] rules = new XMLSyntaxRule[] {
				AttributeRule.newIntegerRule(LOG_EVERY),
				new StringAttributeRule(FILE_NAME,
						"The name of the file to send log output to. " +
								"If no file name is specified then log is sent to standard output", true),
				new StringAttributeRule(TITLE, "The title of the log", true),
				AttributeRule.newBooleanRule(NEXUS_FORMAT, true,
						"Whether to use the NEXUS format for the tree log"),
				new StringAttributeRule(BRANCH_LENGTHS, "What units should the branch lengths be in", new String[] { TIME, SUBSTITUTIONS }, true),
				new ElementRule(Tree.class, "The tree which is to be logged"),
				new ElementRule(BranchRateController.class, true),
				new ElementRule(TreeColouringProvider.class, true),
				new ElementRule(Likelihood.class, true),
				new ElementRule(TreeAttributeProvider.class, 1, Integer.MAX_VALUE),
				new ElementRule(NodeAttributeProvider.class, 1, Integer.MAX_VALUE),
				new ElementRule(BranchAttributeProvider.class, 1, Integer.MAX_VALUE)
		};

		public String getParserDescription() {
			return "Logs a tree to a file";
		}

		public String getExample() {
			return
					"<!-- The " + getParserName() + " element takes a treeModel to be logged -->\n" +
							"<" + getParserName() + " " + LOG_EVERY + "=\"100\" " + FILE_NAME + "=\"log.trees\" " + NEXUS_FORMAT + "=\"true\">\n" +
							"	<treeModel idref=\"treeModel1\"/>\n" +
							"</" + getParserName() + ">\n";
		}

		public Class getReturnType() { return MLLogger.class; }
	};

}