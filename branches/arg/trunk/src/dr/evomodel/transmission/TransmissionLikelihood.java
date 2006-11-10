/*
 * TransmissionLikelihood.java
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

package dr.evomodel.transmission;

import dr.evolution.coalescent.Coalescent;
import dr.evolution.coalescent.DemographicFunction;
import dr.evolution.coalescent.Intervals;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxon;
import dr.evolution.util.TaxonList;
import dr.evolution.util.Units;
import dr.evomodel.coalescent.DemographicModel;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.AbstractModel;
import dr.inference.model.Likelihood;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * A likelihood function for a complete transmission history. Takes a viruses tree
 * and a demographic model. The transmission history consists of a number of
 * hosts with known history of transmission. The viruses tree should have tip
 * attributes specifying which host they are from (host="").
 *
 * @version $Id: TransmissionLikelihood.java,v 1.13 2005/06/15 17:20:54 rambaut Exp $
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 */
public class TransmissionLikelihood extends AbstractModel implements Likelihood, Units {

	// PUBLIC STUFF

	public static final String TRANSMISSION_LIKELIHOOD = "transmissionLikelihood";
	public static final String SOURCE_PATIENT = "sourcePatient";

	public TransmissionLikelihood(Tree hostTree, Tree virusTree,
									DemographicModel sourceDemographic,
									TransmissionDemographicModel transmissionModel)
				throws TaxonList.MissingTaxonException {
		this(TRANSMISSION_LIKELIHOOD, hostTree, virusTree, sourceDemographic, transmissionModel);
	}

	public TransmissionLikelihood(String name, Tree hostTree, Tree virusTree,
									DemographicModel sourceDemographic,
									TransmissionDemographicModel transmissionModel)
				throws TaxonList.MissingTaxonException {

		super(name);

		this.hostTree = hostTree;
		if (hostTree instanceof TreeModel) {
			addModel((TreeModel)hostTree);
		}

		this.virusTree = virusTree;
		if (virusTree instanceof TreeModel) {
			addModel((TreeModel)virusTree);
		}

		this.sourceDemographic = sourceDemographic;
		addModel(sourceDemographic);

		this.transmissionModel = transmissionModel;
		addModel(transmissionModel);

		for (int i = 0; i < virusTree.getExternalNodeCount(); i++) {
			Taxon hostTaxon = (Taxon)virusTree.getTaxonAttribute(i, "host");
			if (hostTaxon == null) throw new TaxonList.MissingTaxonException("One or more of the viruses tree's taxa are missing the 'host' attribute");

			int host = hostTree.getTaxonIndex(hostTaxon);
			if (host == -1) throw new TaxonList.MissingTaxonException("One of the viruses tree's host attribute, " +
					hostTaxon.getId() + ", was not found as a taxon in the host tree");
		}

		setupHosts();
	}

	private void setupHosts() {

		hostCount = hostTree.getTaxonCount();

		intervals = new Intervals[hostCount];
		for (int i = 0; i < hostCount; i++) {
			// 3 times virusTree tip count will be enough events...
			intervals[i] = new Intervals(virusTree.getExternalNodeCount() * 3);
		}
		donorHost = new int[hostCount];
		donorHost[0] = -1;
		transmissionTime = new double[hostCount];
		transmissionTime[0] = Double.POSITIVE_INFINITY;
		donorSize = new double[hostCount];

		setupHosts(hostTree.getRoot());
	}

	private int setupHosts(NodeRef node) {

		int host;

		if (hostTree.isExternal(node)) {
			host = node.getNumber();
		} else {

		// This traversal assumes that the first child is the donor
		// and the second is the recipient

			int host1 = setupHosts(hostTree.getChild(node, 0));
			int host2 = setupHosts(hostTree.getChild(node, 1));

			donorHost[host2] = host1;
			transmissionTime[host2] = hostTree.getNodeHeight(node);

			host = host1;
		}

		return host;
	}

	// **************************************************************
    // ModelListener IMPLEMENTATION
    // **************************************************************

	protected final void handleModelChangedEvent(Model model, Object object, int index) {
		if (model == virusTree) {
			// treeModel has changed so recalculate the intervals
		} else if (model == hostTree) {
			// hosts treeModel has changed so recalculate the hosts and intervals
		} else {
			// demographicModel has changed so we don't need to recalculate the intervals
		}

		likelihoodKnown = false;
	}

	// **************************************************************
    // ParameterListener IMPLEMENTATION
    // **************************************************************

	protected final void handleParameterChangedEvent(Parameter parameter, int index) { } // No parameters to respond to

	// **************************************************************
    // Model IMPLEMENTATION
    // **************************************************************

	/**
	 * Stores the precalculated state: in this case the intervals
	 */
	protected final void storeState() {
	}

	/**
	 * Restores the precalculated state: that is the intervals of the tree.
	 */
	protected final void restoreState() {
		likelihoodKnown = false;
	}

	protected final void acceptState() { } // nothing to do
	protected final void adoptState(Model source) { }

	// **************************************************************
    // Likelihood IMPLEMENTATION
    // **************************************************************

	public final Model getModel() { return this; }

	public final double getLogLikelihood() {
		if (!likelihoodKnown) {
			logLikelihood = calculateLogLikelihood();
			likelihoodKnown = true;
		}
		return logLikelihood;
	}

	public final void makeDirty() {
		likelihoodKnown = false;
	}

	/**
	 * Calculates the log likelihood of this set of coalescent intervals,
	 * given a demographic model.
	 */
	public double calculateLogLikelihood() {

		makeDirty();
		setupHosts();
		for (int i = 0; i < hostCount; i++) {
			intervals[i].resetEvents();
			donorSize[i] = -1;
		}
		try {
			setupIntervals(virusTree.getRoot());
		} catch (IncompatibleException re) {
			// register the compatibility failure
			return Double.NEGATIVE_INFINITY;
		}

		/*

		if (intervalsKnown == false) {

			if (!hostsKnown) {
				setupHosts();
				hostsKnown = true;
			}

			for (int i = 0; i < hostCount; i++) {
				intervals[i].resetEvents();
				donorSize[i] = -1;
			}
			try {
				setupIntervals(virusTree.getRoot());

				if (isCompatible(virusTree.getRoot()) == -1) {
					System.out.println("compatibility failed!");
					makeDirty();
					// register the compatibility failure
					return Double.NEGATIVE_INFINITY;
				}

				//System.out.println("intervals set up successfully!");
				intervalsKnown = true;
				savedHostTree = new FlexibleTree(hostTree);
			} catch (IncompatibleException re) {
				System.out.println("intervals setup failed!");
				if (savedHostTree == null) {
					throw new RuntimeException(re.getMessage());
				} else {
					makeDirty();
					// register the compatibility failure
					return Double.NEGATIVE_INFINITY;
				}
			}
		}*/

		for (int i = 0; i < hostCount; i++) {
			donorSize[i] = -1;
		}

		DemographicFunction demoFunction = sourceDemographic.getDemographicFunction();
		double logL = Coalescent.calculateLogLikelihood(intervals[0], demoFunction);

		for (int i = 1; i < hostCount; i++) {
			double ds = getDonorSize(i);
			demoFunction = transmissionModel.getDemographicFunction(transmissionTime[i], ds);
			logL += Coalescent.calculateLogLikelihood(intervals[i], demoFunction);

		}

		return logL;
	}

	private double getDonorSize(int host) {
		if (donorSize[host] > 0.0) {
			return donorSize[host];
		}

		DemographicFunction demoFunction;

		if (donorHost[host] == 0) {
			demoFunction = sourceDemographic.getDemographicFunction();
		} else {
			double ds = getDonorSize(donorHost[host]);
			demoFunction = transmissionModel.getDemographicFunction(transmissionTime[host], ds);
		}

		donorSize[host] = demoFunction.getDemographic(transmissionTime[host]);
		return donorSize[host];
	}

	private int setupIntervals(NodeRef node) throws IncompatibleException {

		double height = virusTree.getNodeHeight(node);
		int host;

		if (virusTree.isExternal(node)) {
			Taxon hostTaxon = (Taxon)virusTree.getTaxonAttribute(node.getNumber(), "host");
			host = hostTree.getTaxonIndex(hostTaxon);

			intervals[host].addSampleEvent(height);
		} else {

			// Tree should be bifurcating...
			int host1 = setupIntervals(virusTree.getChild(node, 0));
			int host2 = setupIntervals(virusTree.getChild(node, 1));

			while (height > transmissionTime[host1]) {

				double time = transmissionTime[host1];

				intervals[host1].addNothingEvent(time);

				host1 = donorHost[host1];

				intervals[host1].addSampleEvent(time);
			}

			while (height > transmissionTime[host2]) {

				double time = transmissionTime[host2];

				intervals[host2].addNothingEvent(time);

				host2 = donorHost[host2];

				intervals[host2].addSampleEvent(time);
			}

			if (host1 != host2) {
				throw new IncompatibleException("Virus tree is not compatible with transmission history");
			}

			host = host1;
			intervals[host].addCoalescentEvent(height);

		}

		return host;
	}

    // **************************************************************
    // Loggable IMPLEMENTATION
    // **************************************************************

	/**
	 * @return the log columns.
	 */
	public final dr.inference.loggers.LogColumn[] getColumns() {
		return new dr.inference.loggers.LogColumn[] {
			new LikelihoodColumn(getId())
		};
	}

	private final class LikelihoodColumn extends dr.inference.loggers.NumberColumn {
		public LikelihoodColumn(String label) { super(label); }
		public double getDoubleValue() { return getLogLikelihood(); }
	}

    // **************************************************************
    // XMLElement IMPLEMENTATION
    // **************************************************************

	public org.w3c.dom.Element createElement(org.w3c.dom.Document d) {
		throw new RuntimeException("createElement not implemented");
	}

    // **************************************************************
    // Units IMPLEMENTATION
    // **************************************************************

	/**
	 * Sets the units these coalescent intervals are
	 * measured in.
	 */
	public final void setUnits(int u)
	{
		transmissionModel.setUnits(u);
	}

	/**
	 * Returns the units these coalescent intervals are
	 * measured in.
	 */
	public final int getUnits()
	{
		return transmissionModel.getUnits();
	}

	// ****************************************************************
	// Private and protected stuff
	// ****************************************************************

	public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

		public String getParserName() { return TRANSMISSION_LIKELIHOOD; }

		public Object parseXMLObject(XMLObject xo) throws XMLParseException {

			DemographicModel demoModel0 = (DemographicModel)xo.getSocketChild(SOURCE_PATIENT);
			TransmissionDemographicModel demoModel1 = (TransmissionDemographicModel)xo.getChild(TransmissionDemographicModel.class);
			Tree hostTree = (Tree)xo.getSocketChild("hostTree");
			Tree virusTree = (Tree)xo.getSocketChild("parasiteTree");

			TransmissionLikelihood likelihood = null;
			try {
				likelihood = new TransmissionLikelihood(hostTree, virusTree, demoModel0, demoModel1);
			} catch (TaxonList.MissingTaxonException e) {
				throw new XMLParseException(e.toString());
			}
			return likelihood;
		}

		//************************************************************************
		// AbstractXMLObjectParser implementation
		//************************************************************************

		public String getParserDescription() {
			return "This element represents a likelihood function for transmission.";
		}

		public Class getReturnType() { return TransmissionLikelihood.class; }

		public XMLSyntaxRule[] getSyntaxRules() { return rules; }

		private XMLSyntaxRule[] rules = new XMLSyntaxRule[] {
			new ElementRule(SOURCE_PATIENT, DemographicModel.class,
				"This describes the demographic process for the source donor patient."),
			new ElementRule(TransmissionDemographicModel.class,
				"This describes the demographic process for the recipient patients."),
			new ElementRule("hostTree",
				new XMLSyntaxRule[] { new ElementRule(Tree.class) }),
			new ElementRule("parasiteTree",
				new XMLSyntaxRule[] { new ElementRule(Tree.class) })
		};
	};

	class IncompatibleException extends Exception {

		/**
		 * 
		 */
		private static final long serialVersionUID = 8439923064799668934L;

		public IncompatibleException(String name) {
			super(name);
		}
	}



	/** The demographic models. */
	private DemographicModel sourceDemographic = null;
	private TransmissionDemographicModel transmissionModel = null;

	/** The host tree. */
	private Tree hostTree = null;

	/** The viruses tree. */
	private Tree virusTree = null;

	/** The number of hosts. */
	private int hostCount;

	/** The intervals for each host. */
	private Intervals[] intervals;

	/** The donor host for each recipient host (-1 for initial host). */
	private int[] donorHost;

	/** The time of transmission into this host (POSITIVE_INFINITY for initial host). */
	private double[] transmissionTime;

	/** The size of the donor population at time of transmission into recipient host. */
	private double[] donorSize;

	private boolean likelihoodKnown = false;
	private double logLikelihood;
}