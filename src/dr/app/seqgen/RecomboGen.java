/*
 * RecomboGen.java
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

package dr.app.seqgen;

import dr.evolution.alignment.Alignment;
import dr.evolution.alignment.SimpleAlignment;
import dr.evolution.coalescent.*;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.sequence.Sequence;
import dr.evolution.tree.*;
import dr.evolution.util.*;
import dr.evomodel.branchratemodel.*;
import dr.oldevomodel.sitemodel.GammaSiteModel;
import dr.oldevomodel.sitemodel.SiteModel;
import dr.oldevomodel.substmodel.FrequencyModel;
import dr.oldevomodel.substmodel.HKY;
import dr.inference.model.Parameter;
import dr.math.MathUtils;
import dr.xml.*;

import java.util.*;

/**
 * @author Andrew Rambaut
 * @version $Id$
 */
public class RecomboGen {

    public RecomboGen(final double recombinationRate, final int length,
                      final double ancestralPopulationSize,
                      final SiteModel siteModel,
                      final BranchRateModel branchRateModel,
                      final Taxa taxa) {
        this.recombinationRate = recombinationRate;
        this.length = length;
        this.ancestralPopulationSize = ancestralPopulationSize;
        this.siteModel = siteModel;
        this.branchRateModel = branchRateModel;
        this.taxa = taxa;
        calculateHeights();
    }

//    public RecomboGen(final double recombinationRate, final int length,
//                      final double ancestralPopulationSize,
//                      final SiteModel siteModel,
//                      final BranchRateModel branchRateModel,
//                      final String[] taxonNames, final double[] samplingTimes) {
//        this.recombinationRate = recombinationRate;
//        this.length = length;
//        this.ancestralPopulationSize = ancestralPopulationSize;
//        this.siteModel = siteModel;
//        this.branchRateModel = branchRateModel;
//        this.taxonNames = taxonNames;
//        this.samplingTimes = samplingTimes;
//    }

    private void calculateHeights() {

        dr.evolution.util.Date mostRecent = null;
        boolean usingDates = false;

        for (int i = 0; i < taxa.getTaxonCount(); i++) {
            if (TaxonList.Utils.hasAttribute(taxa, i, dr.evolution.util.Date.DATE)) {
                usingDates = true;
                dr.evolution.util.Date date = (dr.evolution.util.Date)taxa.getTaxonAttribute(i, dr.evolution.util.Date.DATE);
                if ((date != null) && (mostRecent == null || date.after(mostRecent))) {
                    mostRecent = date;
                }
            } else {
                // assume contemporaneous tips
                taxa.setTaxonAttribute(i, "height", 0.0);
            }
        }

        if (usingDates && mostRecent != null ) {
            TimeScale timeScale = new TimeScale(mostRecent.getUnits(), true, mostRecent.getAbsoluteTimeValue());

            for (int i =0; i < taxa.getTaxonCount(); i++) {
                dr.evolution.util.Date date = (dr.evolution.util.Date)taxa.getTaxonAttribute(i, dr.evolution.util.Date.DATE);

                if (date == null) {
                    throw new IllegalArgumentException("Taxon, " + taxa.getTaxonId(i) + ", is missing its date");
                }

                taxa.setTaxonAttribute(i, "height", timeScale.convertTime(date.getTimeValue(), date));
            }
        }
    }

    public Alignment generate() {
        List<Node> tips = new ArrayList<Node>();

        // add all the tips
        for (int i = 0; i < taxa.getTaxonCount(); i++) {
            tips.add(new Node(taxa.getTaxon(i)));
        }
        Collections.sort(tips);

        List<Node> unusedTips = new ArrayList<Node>(tips);

        double time = 0;
        List<Node> nodes = new ArrayList<Node>();

        // Add any tips with zero sampling time.
        List<Node> nodesToAdd = new ArrayList<Node>();
        for (Node tip : unusedTips) {
            if (tip.getTime() == 0.0) {
                nodesToAdd.add(tip);
                System.out.println("Time: " + time + " adding " + tip.getTaxon());
            }
        }
        nodes.addAll(nodesToAdd);
        unusedTips.removeAll(nodesToAdd);

        do {

            boolean tipsAdded;

            do {
                tipsAdded = false;

                final int lineageCount = nodes.size();

                // get time to next event...
                final double U = MathUtils.nextDouble(); // create unit uniform random variate
                final double interval = -Math.log(U)/ (lineageCount * recombinationRate);
                final double nextTime = time + interval;

                // Add any tips for which we have reached their sampling time.
                nodesToAdd.clear();
                for (Node tip : unusedTips) {
                    if (tip.getTime() <= nextTime) {
                        nodesToAdd.add(tip);
                        tipsAdded = true;

                        System.out.println("Time: " + tip.getTime() + " adding " + tip.getTaxon());

                        // reset the current time (the tips are sorted in time order
                        // so this will be the oldest added tip).
                        time = tip.getTime();
                    }
                }
                nodes.addAll(nodesToAdd);
                unusedTips.removeAll(nodesToAdd);

                if (!tipsAdded) {
                    time = nextTime;
                }

            } while(tipsAdded); // only continue when no tips are added

            int r = MathUtils.nextInt(nodes.size());
            Node node = nodes.get(r);

            // create two new parent nodes
            Node parent1 = new Node(node, time);
            Node parent2 = new Node(node, time);

            // select a break point in interval [1, length - 2] on
            // a zero-indexed line.
            int breakPoint = MathUtils.nextInt(length - 2) + 1;

            // setup child node with parents and break point
            node.setParent1(parent1);
            node.setParent2(parent2);
            node.setBreakPoint(breakPoint);

            System.out.println("Time: " + time + " recombining " + (node.getTaxon() != null ? node.getTaxon() : r) + " at breakpoint " + breakPoint);

            nodes.add(parent1);
            nodes.add(parent2);
            nodes.remove(node);

        } while (unusedTips.size() > 0);

        // Construct a taxon set for coalescent simulation of deep tree
        Taxa treeTaxa = new Taxa();
        int i = 0;

        Map<Node, Taxon> nodeMap = new HashMap<Node, Taxon>();

        for (Node node : nodes) {
            Taxon taxon = new Taxon("Taxon" + i);
            treeTaxa.addTaxon(taxon);
            nodeMap.put(node, taxon);
            i++;
        }

        CoalescentSimulator coalSim = new CoalescentSimulator();
        ConstantPopulation demo = new ConstantPopulation(Units.Type.YEARS);
        demo.setN0(ancestralPopulationSize);

        Tree tree = coalSim.simulateTree(treeTaxa, demo);

        System.out.println("Tree MRCA " + tree.getNodeHeight(tree.getRoot()) + time);

        SequenceSimulator seqSim = new SequenceSimulator(tree, siteModel, branchRateModel, length);

        Alignment ancestralAlignment = seqSim.simulate();

        SimpleAlignment alignment = new SimpleAlignment();

        // now construct the recombinant sequences from this alignment using the previously
        // generated recombinant history
        for (Node tip : tips) {
            String seqString = generateRecombinant(tip, nodeMap, ancestralAlignment);
            Sequence sequence = new Sequence(tip.getTaxon(), seqString);
//            System.out.println(">" + tip.getTaxon() + "\r" + sequence);

            alignment.addSequence(sequence);
        }

        return alignment;
    }

    private String generateRecombinant(final Node node, final Map<Node, Taxon> nodeMap, final Alignment alignment) {
        String seq = node.getSequence();

        if (seq == null) {
            // if the sequence hasn't already been cached, then construct it...
            if (node.getParent1() == null && node.getParent2() == null) {
                // no parents so must have a sequence
                Taxon taxon = nodeMap.get(node);
                int index = alignment.getTaxonIndex(taxon);
                Sequence sequence = alignment.getSequence(index);
                seq = sequence.getSequenceString();
            } else {
                String part1 = generateRecombinant(node.getParent1(), nodeMap, alignment);
                String part2 = generateRecombinant(node.getParent2(), nodeMap, alignment);

                int breakPoint = node.getBreakPoint();

                seq = part1.substring(0, breakPoint) + part2.substring(breakPoint);
            }

            // cache the sequence
            node.setSequence(seq);
        }

        return seq;
    }

    class Node implements Comparable<Node>{
        Node(final Taxon taxon) {
            this.parent1 = null;
            this.parent2 = null;
            this.child = null;
            this.breakPoint = -1;
            this.taxon = taxon;
            this.time = (Double)taxon.getAttribute("height");
        }

        Node(final Node child, final double time) {
            this.child = child;
            this.parent1 = null;
            this.parent2 = null;
            this.time = time;
            this.breakPoint = -1;
            this.taxon = null;
        }

        public Node getParent1() {
            return parent1;
        }

        public Node getParent2() {
            return parent2;
        }

        public void setParent1(final Node parent1) {
            this.parent1 = parent1;
        }

        public void setParent2(final Node parent2) {
            this.parent2 = parent2;
        }

        public Node getChild() {
            return child;
        }

        public double getTime() {
            return time;
        }

        public int getBreakPoint() {
            return breakPoint;
        }

        public void setBreakPoint(final int breakPoint) {
            this.breakPoint = breakPoint;
        }

        public Taxon getTaxon() {
            return taxon;
        }

        public String getSequence() {
            return sequence;
        }

        public void setSequence(final String sequence) {
            this.sequence = sequence;
        }

        public int compareTo(final Node node) {
            return Double.compare(time, node.time);
        }

        private Node parent1;
        private Node parent2;
        private final Node child;
        private final double time;
        private int breakPoint;
        private final Taxon taxon;
        private String sequence = null;
    }

    /* standard xml parser stuff follows */
    public static final String RECOMBINATION_SIMULATOR = "recombinationSimulator";
    public static final String SITE_MODEL = SiteModel.SITE_MODEL;
    public static final String TAXA = "taxa";
//    public static final String BRANCH_RATE_MODEL = "branchRateModel";
    public static final String RECOMBINATION_RATE = "recombinationRate";
    public static final String SUBSTITUTION_RATE = "substitutionRate";
    public static final String ANCESTRAL_POPULATION_SIZE = "ancestralPopulationSize";
    public static final String SEQUENCE_LENGTH = "sequenceLength";

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return RECOMBINATION_SIMULATOR;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            double recombinationRate = xo.getDoubleAttribute(RECOMBINATION_RATE);
            double ancestralPopulationSize = xo.getDoubleAttribute(ANCESTRAL_POPULATION_SIZE);
            int sequenceLength = xo.getIntegerAttribute(SEQUENCE_LENGTH);

            Taxa taxa = (Taxa) xo.getChild(Taxa.class);
            SiteModel siteModel = (SiteModel) xo.getChild(SiteModel.class);
            BranchRateModel rateModel = (BranchRateModel)xo.getChild(BranchRateModel.class);

            if (rateModel == null) {
                rateModel = new DefaultBranchRateModel();
            }

            RecomboGen s = new RecomboGen(
                    recombinationRate, sequenceLength,
                    ancestralPopulationSize, siteModel, rateModel, taxa);
            return s.generate();
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "A SequenceSimulator that generates random sequences for a given tree, siteratemodel and branch rate model";
        }

        public Class getReturnType() {
            return Alignment.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                new ElementRule(Taxa.class),
                new ElementRule(SiteModel.class),
                new ElementRule(BranchRateModel.class, true),
                AttributeRule.newDoubleRule(RECOMBINATION_RATE),
                AttributeRule.newDoubleRule(ANCESTRAL_POPULATION_SIZE),
                AttributeRule.newIntegerRule(SEQUENCE_LENGTH)
        };
    };

    public static void main(String[] argv) {
        // Simulate sequences on this tree to generate sequences at the top of the
        // recombination process.
        Parameter kappa = new Parameter.Default(1, 2);
        Parameter freqs = new Parameter.Default(new double[]{0.25, 0.25, 0.25, 0.25});
        FrequencyModel freqModel = new FrequencyModel(Nucleotides.INSTANCE, freqs);
        HKY hky = new HKY(kappa, freqModel);

        // create site model
        SiteModel siteModel = new GammaSiteModel(hky);
        // create branch rate model
        Parameter rate = new Parameter.Default(1, 1.0E-4);
        BranchRateModel branchRateModel = new StrictClockBranchRates(rate);

//        RecomboGen recomboGen = new RecomboGen(8.0E-2, 100, 1000, 1.0E-4,
//                new String[] { "taxon1_0", "taxon2_0", "taxon3_10", "taxon4_10", "taxon5_20", "taxon6_20"},
//                new double[] { 0, 0, 10, 10, 20, 20} );
//
//        recomboGen.generate();
    }

    private final double recombinationRate;
    private final int length;
    private final double ancestralPopulationSize;
    private final SiteModel siteModel;
    private final BranchRateModel branchRateModel;
    private final Taxa taxa;
//    private final String[] taxonNames;
//    private final double[] samplingTimes;
}