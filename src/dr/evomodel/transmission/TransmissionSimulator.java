/*
 * TransmissionSimulator.java
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

package dr.evomodel.transmission;

import dr.evolution.coalescent.Intervals;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxa;
import dr.evolution.util.TaxonList;
import dr.evomodel.coalescent.DemographicModel;
import dr.xml.*;


/**
 * A likelihood function for a complete transmission history. Takes a viruses tree
 * and a demographic model. The transmission history consists of a number of
 * hosts with known history of transmission. The viruses tree should have tip
 * attributes specifying which host they are from (host="").
 *
 * @version $Id: TransmissionSimulator.java,v 1.5 2005/06/15 17:20:55 rambaut Exp $
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 */
public class TransmissionSimulator {

	// PUBLIC STUFF

	public static final String TRANSMISSION_SIMULATOR = "transmissionSimulator";
	public static final String SOURCE_PATIENT = "sourcePatient";

	public TransmissionSimulator(Taxa taxa, Tree hostTree, DemographicModel demographic)
				throws TaxonList.MissingTaxonException {
		this(TRANSMISSION_SIMULATOR, taxa, hostTree, demographic);
	}

	public TransmissionSimulator(String name, Taxa taxa, Tree hostTree, DemographicModel demographic)
				throws TaxonList.MissingTaxonException {

		this.hostTree = hostTree;

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

	// ****************************************************************
	// Private and protected stuff
	// ****************************************************************

	public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

		public String getParserName() { return TRANSMISSION_SIMULATOR; }

		public Object parseXMLObject(XMLObject xo) throws XMLParseException {

			DemographicModel demoModel0 = (DemographicModel)xo.getElementFirstChild(SOURCE_PATIENT);
			TransmissionDemographicModel demoModel1 = (TransmissionDemographicModel)xo.getChild(TransmissionDemographicModel.class);
			Tree hostTree = (Tree)xo.getElementFirstChild("hostTree");
			Tree virusTree = (Tree)xo.getElementFirstChild("parasiteTree");

			TransmissionSimulator simulator = null;
//			try {
//				simulator = new TransmissionSimulator(hostTree, virusTree, demoModel0, demoModel1);
//			} catch (TaxonList.MissingTaxonException e) {
//				throw new XMLParseException(e.toString());
//			}
			return simulator;
		}

		//************************************************************************
		// AbstractXMLObjectParser implementation
		//************************************************************************

		public String getParserDescription() {
			return "This element represents a simulator for a transmission history.";
		}

		public Class getReturnType() { return TransmissionSimulator.class; }

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

}