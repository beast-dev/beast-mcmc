/*
 * TKF91Likelihood.java
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

package dr.evomodel.indel;

import dr.evolution.alignment.Alignment;
import dr.evomodel.sitemodel.GammaSiteModel;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.*;
import dr.xml.*;


/**
 * Calculates the likelihood of a set of continuous attributes on a tree.
 *
 * @version $Id: TKF91Likelihood.java,v 1.17 2005/06/20 07:24:25 rambaut Exp $
 *
 * @author Alexei Drummond
 */
public class TKF91Likelihood extends AbstractModel implements Likelihood {

	public static final String TKF91_LIKELIHOOD = "tkf91Likelihood";
	public static final String TKF91_LENGTH_DIST = "lengthDistribution";
	public static final String TKF91_DEATH = "deathRate";
	public static final String MU = "mutationRate";

	public TKF91Likelihood(TreeModel treeModel, Alignment alignment, GammaSiteModel siteModel, TKF91Model tkfModel) {

		super(TKF91_LIKELIHOOD);

		if (siteModel.getAlphaParameter() != null) throw new IllegalArgumentException("TKF91 model cannot handle gamma-distributed rates");
		if (siteModel.getPInvParameter() != null) throw new IllegalArgumentException("TKF91 model cannot handle invariant sites");

		addModel(siteModel);
		addModel(tkfModel);
		addModel(treeModel);

		this.treeModel = treeModel;
		this.alignment = alignment;
		this.siteModel = siteModel;
		this.tkfModel = tkfModel;

		recursion = new dr.evomodel.indel.HomologyRecursion();
		recursion.init(
			treeModel,
			alignment,
			siteModel.getSubstitutionModel(),
			siteModel.getMutationRateParameter().getParameterValue(0),
			tkfModel.getLengthDistributionValue(),
			tkfModel.getDeathRate(1));

		addStatistic(new AlignmentLengthStatistic());
	}

	public class AlignmentLengthStatistic extends Statistic.Abstract {

		public double getStatisticValue(int dim) { return alignment.getSiteCount(); }

		public int getDimension() { return 1; }

		public String getStatisticName() { return "alignmentLength"; }
	}

	public void adoptState(Model model) {
		throw new RuntimeException("Not implemented!");
	}

	public void acceptState() {
		//throw new RuntimeException("Not implemented!");
	}

	public void handleParameterChangedEvent(Parameter p, int i) {
		//throw new RuntimeException("Not implemented!");
	}

	public void handleModelChangedEvent(Model m, Object o, int i) {
		//throw new RuntimeException("Not implemented!");
	}

	/**
	 * Sets the alignment of this tkf91likelihood
	 */
	public void setAlignment(Alignment a) {
		alignment = a;
		//System.out.println("Set new alignment");
	}

	protected void storeState() {
		//super.storeState();
		//storedAlignment = alignment;
		//System.out.println("Stored alignment");

	}

	protected void restoreState() {
		//super.restoreState();
		//alignment = storedAlignment;
		//System.out.println("restored alignment");
	}

	public Alignment getAlignment() {
		return alignment;
	}

	public GammaSiteModel getSiteModel() {
		return siteModel;
	}

	public TreeModel getTreeModel() {
		return treeModel;
	}

	public void makeDirty() { // this is always dirty
	}

	public Model getModel() { return this; }

	/**
	 * Get the log likelihood.
	 * @return the log likelihood.
	 */
	public double getLogLikelihood() {

		recursion.init(
			treeModel,
			alignment,
			siteModel.getSubstitutionModel(),
			siteModel.getMutationRateParameter().getParameterValue(0),
			tkfModel.getLengthDistributionValue(),
			tkfModel.getDeathRate(1));

		double logL = recursion.recursion();

		//System.out.println("logL = " + logL);
		return logL;
	}

	// **************************************************************
    // Loggable IMPLEMENTATION
    // **************************************************************

	/**
	 * @return the log columns.
	 */
	public dr.inference.loggers.LogColumn[] getColumns() {
		return new dr.inference.loggers.LogColumn[] {
			new LikelihoodColumn(getId())
			//, new AlignmentColumn(getId()+".alignment")
		};
	}

	private class LikelihoodColumn extends dr.inference.loggers.NumberColumn {
		public LikelihoodColumn(String label) { super(label); }
		public double getDoubleValue() { return getLogLikelihood(); }
	}

	private class AlignmentColumn extends dr.inference.loggers.LogColumn.Abstract {
		public AlignmentColumn(String label) { super(label); }
		public String getFormattedValue() { return alignment.toString(); }
	}

	// **************************************************************
    // Identifiable IMPLEMENTATION
    // **************************************************************

	private String id = null;

	public void setId(String id) { this.id = id; }

	public String getId() { return id; }


	public String toString() {
		if (id != null) {
			return id;
		}
		return super.toString();
	}

	protected boolean getLikelihoodKnown() {
		return false;
	}

	public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

		public String getParserName() { return TKF91_LIKELIHOOD; }

		public Object parseXMLObject(XMLObject xo) throws XMLParseException {

			TreeModel tree = (TreeModel)xo.getChild(TreeModel.class);
			Alignment alignment = (Alignment)xo.getChild(Alignment.class);
			GammaSiteModel siteModel = (GammaSiteModel)xo.getChild(GammaSiteModel.class);
			TKF91Model tkfModel = (TKF91Model)xo.getChild(TKF91Model.class);
			return new TKF91Likelihood(tree, alignment, siteModel, tkfModel);
		}

		//************************************************************************
		// AbstractXMLObjectParser implementation
		//************************************************************************

		public String getParserDescription() {
			return "Returns the total likelihood of a single alignment under the TKF91 model, for a given tree. " +
				"In particular all possible ancestral histories of insertions and deletions leading to the " +
				"alignment of sequences at the tips are taken into account.";
		}

		public Class getReturnType() { return TKF91Likelihood.class; }

		public XMLSyntaxRule[] getSyntaxRules() { return rules; }

		private XMLSyntaxRule[] rules = new XMLSyntaxRule[] {
			new ElementRule(TreeModel.class),
			new ElementRule(Alignment.class),
			new ElementRule(GammaSiteModel.class),
			new ElementRule(TKF91Model.class)
		};
	};

	private TreeModel treeModel;
	private Alignment alignment;
	//private Alignment storedAlignment;
	private GammaSiteModel siteModel;
	private TKF91Model tkfModel;
	private dr.evomodel.indel.HomologyRecursion recursion = null;
}