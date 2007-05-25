/*
 * LocalClockModel.java
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

package dr.evomodel.branchratemodel;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.util.TaxonList;
import dr.evolution.util.Taxa;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.xml.*;

import java.util.*;

/**
 * Implements a local clock model for rates of evolution on a tree. One or
 * more sets of taxa can be given additional rates. The rate will apply from
 * the TMRCA down (including the branch above if includeStem=true). Unless
 * these clades are constrained to be monophyletic, then the behaviour may
 * be unpredictable.
 * 
 * @author Andrew Rambaut
 *
 * @version $Id: LocalClockModel.java,v 1.1 2005/04/05 09:27:48 rambaut Exp $
 */
public class LocalClockModel extends AbstractModel implements BranchRateModel  {

	public static final String LOCAL_CLOCK_MODEL = "localClockModel";
	public static final String CLADE = "clade";
	public static final String INCLUDE_STEM = "includeStem";
	public static final String EXTERNAL_BRANCHES = "externalBranches";
	public static final String RATE = "rate";

	private final Parameter rateParameter;
	private TreeModel treeModel;
	private ArrayList localClocks = new ArrayList();
	private final Map localRateParameters = new HashMap();


	public LocalClockModel(TreeModel treeModel, Parameter rateParameter) {

		super(LOCAL_CLOCK_MODEL);
		this.treeModel = treeModel;

		addModel(treeModel);

		this.rateParameter = rateParameter;

		addParameter(rateParameter);
	}

	private void addExternalBranchClock(Parameter rateParameter, TaxonList taxonList) throws Tree.MissingTaxonException {
		Set leafSet = Tree.Utils.getLeavesForTaxa(treeModel, taxonList);
		LocalClock clock = new LocalClock(rateParameter, leafSet);
		localClocks.add(clock);

	}

	private void addCladeClock(Parameter rateParameter, TaxonList taxonList, boolean includeStem) throws Tree.MissingTaxonException {
		Set leafSet = Tree.Utils.getLeavesForTaxa(treeModel, taxonList);
		LocalClock clock = new LocalClock(rateParameter, leafSet, includeStem);
		localClocks.add(clock);
	}

	public void handleModelChangedEvent(Model model, Object object, int index) {
		if (model == treeModel) {
			localRateParameters.clear();
			findRateParameters(treeModel, treeModel.getRoot(), null);
		}
		fireModelChanged();
	}

	protected void handleParameterChangedEvent(Parameter parameter, int index) {
		fireModelChanged();
	}

	protected void storeState() {
	}

	protected void restoreState() {
		localRateParameters.clear();
		findRateParameters(treeModel, treeModel.getRoot(), null);
	}

	protected void acceptState() {}

	public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

		public String getParserName() { return LOCAL_CLOCK_MODEL; }

		public Object parseXMLObject(XMLObject xo) throws XMLParseException {

			TreeModel tree = (TreeModel)xo.getChild(TreeModel.class);

			Parameter rateParameter = (Parameter)xo.getSocketChild(RATE);

			LocalClockModel localClockModel =  new LocalClockModel(tree, rateParameter);

			for (int i = 0; i < xo.getChildCount(); i++) {
				if (xo.getChild(i) instanceof XMLObject) {

					XMLObject xoc = (XMLObject)xo.getChild(i);
					if (xoc.getName().equals(CLADE)) {

						Parameter localRateParameter = (Parameter)xoc.getSocketChild(RATE);
						TaxonList taxonList = (TaxonList)xoc.getChild(TaxonList.class);

						if (taxonList.getTaxonCount()==1) {
							throw new XMLParseException("A local clock for a clade must be defined by at least two taxa");
						}

						boolean includeStem = false;

						if (xoc.hasAttribute(INCLUDE_STEM)) {
							includeStem = xoc.getBooleanAttribute(INCLUDE_STEM);
						}

						try {
							localClockModel.addCladeClock(localRateParameter, taxonList, includeStem);

						} catch (Tree.MissingTaxonException mte) {
							throw new XMLParseException("Taxon, " + mte + ", in " + getParserName() + " was not found in the tree.");
						}
					} else if (xoc.getName().equals(EXTERNAL_BRANCHES)) {

						Parameter localRateParameter = (Parameter)xoc.getSocketChild(RATE);
						TaxonList taxonList = (TaxonList)xoc.getChild(TaxonList.class);

						try {
							localClockModel.addExternalBranchClock(localRateParameter, taxonList);

						} catch (Tree.MissingTaxonException mte) {
							throw new XMLParseException("Taxon, " + mte + ", in " + getParserName() + " was not found in the tree.");
						}
					}

				}
			}

			System.out.println("Using local clock branch rate model.");

			return localClockModel;
		}

		//************************************************************************
		// AbstractXMLObjectParser implementation
		//************************************************************************

		public String getParserDescription() {
			return
					"This element returns a branch rate model that adds a delta to each terminal branch length.";
		}

		public Class getReturnType() { return LocalClockModel.class; }

		public XMLSyntaxRule[] getSyntaxRules() { return rules; }

		private XMLSyntaxRule[] rules = new XMLSyntaxRule[] {
				new ElementRule(Parameter.class, "The global rate parameter"),
				new ElementRule(TreeModel.class, "The tree model"),
				new ElementRule(EXTERNAL_BRANCHES,
						new XMLSyntaxRule[] {
								new ElementRule(Taxa.class, "A local clock that will be applied only to the external branches for these taxa"),
								new ElementRule(RATE, Parameter.class, "The local rate parameter", false),
						}, 0, Integer.MAX_VALUE),
				new ElementRule(CLADE,
						new XMLSyntaxRule[] {
								AttributeRule.newBooleanRule(INCLUDE_STEM, true, "determines whether or not the stem branch above this clade is included in the siteModel."),
								new ElementRule(Taxa.class, "A set of taxa which defines a clade to apply a different site model to"),
								new ElementRule(RATE, Parameter.class, "The local rate parameter", false),
						}, 0, Integer.MAX_VALUE)
		};
	};

	public double getBranchRate(Tree tree, NodeRef node) {

		if (tree.isRoot(node)) {
			throw new IllegalArgumentException("root node doesn't have a rate!");
		}

		Parameter localRateParameter = (Parameter) localRateParameters.get(node);

		if (localRateParameter != null) {
			return localRateParameter.getParameterValue(0);
		}

		return rateParameter.getParameterValue(0);
	}

	private void findRateParameters(Tree tree, NodeRef node, Parameter currentParameter) {

		Parameter newParameter = currentParameter;

		for (int j = 0; j < localClocks.size(); j++) {
			LocalClock local = (LocalClock)localClocks.get(j);

			if (local.findMRCA() == node) {
				newParameter = local.getRateParameter();
				if (local.includeStem()) {
					currentParameter = newParameter;
				}
			}
		}

		if (!tree.isExternal(node)) {
			for (int i = 0; i < tree.getChildCount(node); i++) {
				findRateParameters(tree, tree.getChild(node, i), newParameter);
			}
		}

		if (currentParameter != null) {
			localRateParameters.put(node, currentParameter);
		}
	}

	private class LocalClock {
		LocalClock(Parameter rateParameter, Set leafSet) {
			this.rateParameter = rateParameter;
			this.leafSet = leafSet;
			this.isClade = false;
			this.includeStem = true;
		}

		LocalClock(Parameter rateParameter, Set leafSet, boolean includeStem) {
			this.rateParameter = rateParameter;
			this.leafSet = leafSet;
			this.isClade = true;
			this.includeStem = includeStem;
		}

		NodeRef findMRCA() {
			return Tree.Utils.getCommonAncestorNode(treeModel, leafSet);
		}

		boolean includeStem() { return this.includeStem; }
		boolean isClade() { return this.isClade; }
		Parameter getRateParameter() { return this.rateParameter; }

		Parameter rateParameter;
		Set leafSet;
		boolean isClade;
		boolean includeStem;
	};

}