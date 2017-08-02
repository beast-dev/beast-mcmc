/*
 * TransmissionStatistic.java
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

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxon;
import dr.evomodel.tree.TreeStatistic;
import dr.inference.model.BooleanStatistic;
import dr.inference.model.Statistic;
import dr.xml.*;

import java.util.HashSet;
import java.util.Set;


/**
 * A statistic for the compatibility of a viruses tree with a transmission
 * history. The transmission history consists of a number of
 * hosts with known history of transmission. The viruses tree should have tip
 * attributes specifying which host they are from (host="").
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: TransmissionStatistic.java,v 1.11 2005/06/27 21:19:15 rambaut Exp $
 */
public class TransmissionStatistic extends TreeStatistic implements BooleanStatistic {

    // PUBLIC STUFF

    public static final String TRANSMISSION_STATISTIC = "transmissionStatistic";

    public TransmissionStatistic(String name, TransmissionHistoryModel transmissionHistoryModel, Tree virusTree) {

        super(name);

        this.transmissionHistoryModel = transmissionHistoryModel;
        this.virusTree = virusTree;

        setupHosts();
    }

    public TransmissionStatistic(String name, Tree hostTree, Tree virusTree) {

        super(name);

        this.hostTree = hostTree;
        this.virusTree = virusTree;

        setupHosts();
    }

    private void setupHosts() {

        if (transmissionHistoryModel != null) {
            hostCount = transmissionHistoryModel.getHostCount();

        } else {
            hostCount = hostTree.getTaxonCount();
        }

        donorHost = new int[hostCount];
        donorHost[0] = -1;
        transmissionTime = new double[hostCount];
        transmissionTime[0] = java.lang.Double.POSITIVE_INFINITY;

        if (transmissionHistoryModel != null) {
            for (int i = 0; i < transmissionHistoryModel.getTransmissionEventCount(); i++)
            {
                TransmissionHistoryModel.TransmissionEvent event = transmissionHistoryModel.getTransmissionEvent(i);

                int host1 = transmissionHistoryModel.getHostIndex(event.getDonor());
                int host2 = transmissionHistoryModel.getHostIndex(event.getRecipient());

                donorHost[host2] = host1;
                transmissionTime[host2] = event.getTransmissionTime();
            }
        } else {
            setupHostsTree(hostTree.getRoot());
        }
    }

    private int setupHostsTree(NodeRef node) {

        int host;

        if (hostTree.isExternal(node)) {
            host = node.getNumber();
        } else {

            // This traversal assumes that the first child is the donor
            // and the second is the recipient

            int host1 = setupHostsTree(hostTree.getChild(node, 0));
            int host2 = setupHostsTree(hostTree.getChild(node, 1));

            donorHost[host2] = host1;
            transmissionTime[host2] = hostTree.getNodeHeight(node);

            host = host1;
        }

        return host;
    }

    public void setTree(Tree tree) {
        this.virusTree = tree;
    }

    public Tree getTree() {
        return virusTree;
    }

    public String getDimensionName(int dim) {
        String recipient = transmissionHistoryModel.getHost(dim).getId();

        String donor = (donorHost[dim] == -1 ? "" : transmissionHistoryModel.getHost(donorHost[dim]).getId() + "->");
        return "transmission(" + donor + recipient + ")";
    }

    public int getDimension() {
        return hostCount;
    }

    /**
     * @return boolean result of test.
     */
    public double getStatisticValue(int dim) {
        return getBoolean(dim) ? 1.0 : 0.0;
    }

    /**  /**
     * @return true if the population tree is compatible with the species tree
     */
    public boolean getBoolean(int dim) {
        Set<Integer> incompatibleSet = new HashSet<Integer>();
        setupHosts();
        isCompatible(virusTree.getRoot(), incompatibleSet);

        return !incompatibleSet.contains(dim);
    }

    private int isCompatible(NodeRef node, Set<Integer> incompatibleSet) {

        double height = virusTree.getNodeHeight(node);
        int host;

        if (virusTree.isExternal(node)) {
            Taxon hostTaxon = (Taxon) virusTree.getTaxonAttribute(node.getNumber(), "host");
            if (transmissionHistoryModel != null) {
                host = transmissionHistoryModel.getHostIndex(hostTaxon);
            } else {
                host = hostTree.getTaxonIndex(hostTaxon);
            }

            if (host != -1 && height > transmissionTime[host]) {
                // This means that the sequence was sampled
                // before the host was infected so we should probably flag
                // this as an error before we get to this point...
                throw new RuntimeException("Sequence " + virusTree.getNodeTaxon(node) + ", was sampled ("+height+") before host, " + hostTaxon + ", was infected ("+transmissionTime[host]+")");
            }

        } else {

            // Tree should be bifurcating...
            int host1 = isCompatible(virusTree.getChild(node, 0), incompatibleSet);

            int host2 = isCompatible(virusTree.getChild(node, 1), incompatibleSet);

            if (host1 == host2) {

                host = host1;
                while (height > transmissionTime[host]) {
                    host = donorHost[host];
                }

            } else {
                while (height > transmissionTime[host1]) {
                    host1 = donorHost[host1];
                }

                while (height > transmissionTime[host2]) {
                    host2 = donorHost[host2];
                }

                if (host1 != host2) {
                    if (transmissionTime[host1] < transmissionTime[host2]) {
                        incompatibleSet.add(host1);
                        host = host2;
                    } else {
                        incompatibleSet.add(host2);
                        host = host1;
                    }
                } else {
                    host = host1;
                }
            }
        }

        return host;
    }

    // ****************************************************************
    // Private and protected stuff
    // ****************************************************************

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return TRANSMISSION_STATISTIC;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            String name = xo.getStringAttribute("name");
            Tree virusTree = (Tree) xo.getElementFirstChild("parasiteTree");

            if (xo.getChild(TransmissionHistoryModel.class) != null) {
                TransmissionHistoryModel history = (TransmissionHistoryModel) xo.getChild(TransmissionHistoryModel.class);
                return new TransmissionStatistic(name, history, virusTree);
            } else {
                Tree hostTree = (Tree) xo.getElementFirstChild("hostTree");
                return new TransmissionStatistic(name, hostTree, virusTree);
            }

        }

        public String getParserDescription() {
            return "A statistic that returns true if the given parasite tree is compatible with the host tree.";
        }

        public Class getReturnType() {
            return Statistic.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                new StringAttributeRule("name", "A name for this statistic for the purpose of logging"),
                new XORRule(
                        new ElementRule("hostTree",
                                new XMLSyntaxRule[]{new ElementRule(Tree.class)}),
                        new ElementRule(TransmissionHistoryModel.class,
                                "This describes the transmission history of the patients.")
                ),
                new ElementRule("parasiteTree",
                        new XMLSyntaxRule[]{new ElementRule(Tree.class)})
        };
    };

    /**
     * The host tree.
     */
    private Tree hostTree = null;

    private TransmissionHistoryModel transmissionHistoryModel = null;

    /**
     * The viruses tree.
     */
    private Tree virusTree = null;

    /**
     * The number of hosts.
     */
    private int hostCount;

    /**
     * The donor host for each recipient host (-1 for initial host).
     */
    private int[] donorHost;

    /**
     * The time of transmission into this host (POSITIVE_INFINITY for initial host).
     */
    private double[] transmissionTime;
}