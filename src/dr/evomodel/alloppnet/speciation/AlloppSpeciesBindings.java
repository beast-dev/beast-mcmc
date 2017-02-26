/*
 * AlloppSpeciesBindings.java
 *
 * Copyright (c) 2002-2017 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodel.alloppnet.speciation;




import java.util.ArrayList;
import java.util.Arrays;
import java.util.Formatter;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Locale;
import java.util.Map;
import java.util.Set;
import java.util.Stack;

import jebl.util.FixedBitSet;

import dr.evolution.tree.NodeRef;
import dr.evolution.util.Taxon;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.alloppnet.parsers.AlloppSpeciesBindingsParser;
import dr.inference.loggers.LogColumn;
import dr.inference.loggers.Loggable;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.ModelListener;
import dr.inference.model.Variable;
import dr.inference.model.Variable.ChangeType;
import dr.math.MathUtils;
import dr.evomodel.alloppnet.util.AlloppMisc;

/**
 *
 * Knows how species (diploid and allopolyploid species) are made of individuals
 * and individuals are made of taxa (=sequences).
 *
 * @author Graham Jones
 *         Date: 30/04/2011
 */



/*
 *
 * AlloppSpeciesBindings knows how species are made of individuals
 * and how individuals are made of taxa (= diploid genomes within individuals).
 *
 * It also contains the list of gene trees - tree topologies and node
 * times, plus popfactors. Given a AlloppSpeciesNetworkModel
 * it can say if a gene tree is compatible.
 *
 * It is here that assignments of sequence copies within individuals
 * get permuted.
 *
 * The class GeneUnionTree, defined here, is used during calculations.
 *
 * ******************************
 *
 * geneTreeInfos. The array of gene trees. Each one has an array of
 * sequence assignments.
 *
 * apspecies. List of species containing Individuals
 *
 * indivs. A 'flattened' array of all Individuals from all species
 *
 * taxa. A 'flattened' array of all Taxons from all Individuals
 *
 * taxon2index, apspecies2index: maps to indices
 *
 * spsq used to map (species index, sequence index) to a single index.
 *
 * Eg of two individuals from a diploid "a" and two individuals from a tetraploid "b"
 * apspecies: a(01,02), b(03,04)
 * indivs:    01_a(0), 02_a(0), 03_b(0,1), 04_b(0,1)
 * taxa:      01_a0, 02_a0, 03_b0, 03_b1, 04_b0, 04_b1
 *
 */
public class AlloppSpeciesBindings extends AbstractModel implements Loggable {

    private final GeneTreeInfo[] geneTreeInfos;

    private final ApSpInfo[] apspecies;
    private final Taxon[] taxa;
    private final Map<Taxon, Integer> taxon2index = new HashMap<Taxon, Integer>();
    private final int spsq[][];
    private final int numberOfSpSeqs;
    private final double initialmingenenodeheight; // for starting network



    /* ************  subclasses ********************/

    /* class AlloppSpeciesBindings.Individual
     *
     * Individual is a list of Taxons. One for diploid individual, two  for
     * tetraploid, etc.
     */
    public static class Individual extends Taxon {
        final public String id; // individual ID, such as "02_Alpha" in AlloppSpeciesInfoParser XML example
        public final Taxon[] taxa;

        public Individual(String id, Taxon[] taxa) {
            super(id);
            this.id = id;
            this.taxa = taxa;
        }
    }


    private class SpeciesIndivPair {
        public int spIndex;
        public int ivIndex;

        public SpeciesIndivPair(int spIndex, int ivIndex) {
            this.spIndex = spIndex;
            this.ivIndex = ivIndex;
        }
    }


    /* class AlloppSpeciesBindings.ApSpInfo
     *
     * Information on one allopolyploid species
     *
     * name is species name, such as "Alpha" in AlloppSpeciesInfoParser XML example
     *
     * ploidylevel; 2 means diploid, 4 means allotetraploid, etc
     *
     * individuals[]
     */
    public static class ApSpInfo extends Taxon {

        final public String name;
        final public int ploidylevel; // 2 means diploid, 4 means allotetraploid, etc
        final Individual[] individuals;

        public ApSpInfo(String name, int ploidylevel, Individual[] individuals) {
            super(name);
            this.name = name;
            this.individuals = individuals;
            this.ploidylevel = ploidylevel;

            // check
            if (individuals != null) {
                int ntaxaperindiv = ploidylevel / 2;
                for (int ii = 0; ii < individuals.length; ++ii) {
                    assert(individuals[ii].taxa.length == ntaxaperindiv);
                    // may want to allow 3 as well as 2 for tetraploid for organelle DNA
                }
            }
        }


        public Taxon taxonFromIndSeq(int i, int sq) {
            return individuals[i].taxa[sq];
        }
    }


    /* class AlloppSpeciesBindings.GeneTreeInfo
     *
     * Adapted from SpeciesBindings
     *
     * tree is the gene tree
     *
     * seqassigns[]  is species index and sequence index. The latter
     * is what MCMC acts on when permute sequence assignments
     *
     * lineagesCount[] is count for each species, and applies to one or
     * more tips in MulLabTree. I assume this is the same for all genes
     * though the info is stored per gene tree
     *
     * popFactor (normal diploid=2, X chromosome=1.5, chloroplast=.5, etc)
     *
     * I don't like the CoalInfo method of determining compatibility.
     * Instead I copy the gene tree topology and heights into a new
     *  tree with unions in nodes. Maybe I'll find out why JH used CoalInfos...
     */
    private class GeneTreeInfo {
        private final TreeModel tree;
        private SequenceAssignment seqassigns[];
        private SequenceAssignment oldseqassigns[];
        private final int[] lineagesCount;
        private final double popFactor; // grjtodo-oneday will mul pops by this, eg for chloroplast data.


        /* class GeneTreeInfo.SequenceAssignments
         *
         * spIndex is an index for an allopolyploid species. For example, it identifies
         * a bit in a FixedBitSet (union) in a MulLabTree
         *
         * seqIndex identifies a sequence copy for this gene and for each individual.
         * seqIndex is 0 or 1 for tetraploids and it is these that get flipped to change
         * assignments of sequence copies to legs in AlloppSpeciesNetworkModel (or
         * equivalently to tips in a MulLabTree).
         *
         * 2011-06-23 spIndex is the same for all gene trees. Maybe
         * allow non-rectangular data later.
         */
        private class SequenceAssignment {
            public int spIndex;
            public int seqIndex;

            public SequenceAssignment(int spIndex, int seqIndex) {
                this.spIndex = spIndex;
                this.seqIndex = seqIndex;
            }

            public String toString() {
                String s = "" + seqIndex;
                return s;
            }
        }




        /* class GeneTreeInfo.GeneUnionNode
         *
         * Node for GeneTreeInfo.GeneUnionTree.
         */
        private class GeneUnionNode {
            private GeneUnionNode child[];
            private double height;
            private FixedBitSet union;
            private String name; // for debugging

            // Constructor makes a half-formed tip node. Tips need unions
            // and internal nodes need all fields filling in.
            public GeneUnionNode() {
                child = new GeneUnionNode[0];
                height = 0.0;
                union = new FixedBitSet(numberOfSpSeqs());
                name = "";
            }



            public String asText(int indentlen) {
                StringBuilder s = new StringBuilder();
                Formatter formatter = new Formatter(s, Locale.US);
                if (child.length == 0) {
                    formatter.format("%s ", name);
                } else {
                    formatter.format("%s ", "+");
                }
                while (s.length() < 30-indentlen) {
                    formatter.format("%s", " ");
                }
                formatter.format("%s ", AlloppMisc.nonnegIn8Chars(height));
                formatter.format("%20s ", AlloppMisc.FixedBitSetasText(union));

                return s.toString();
            }

        }

        /* class GeneTreeInfo.GeneUnionTree
         *
         * Copy of gene tree topology and heights with unions,
         * implemented as an array of GeneUnionNodes. This is
         * used during calculations fitsInNetwork(), treeLogLikelihood()
         * for one gene tree at a time, then discarded.
         */
        private class GeneUnionTree {
            private GeneUnionNode[] nodes;
            private int nextn;


            public GeneUnionTree() {
                nodes = new GeneUnionNode[tree.getNodeCount()];
                for (int i = 0; i < nodes.length; i++) {
                    nodes[i] = new GeneUnionNode();
                }
                genetree2geneuniontree(tree.getRoot());
            }


            public GeneUnionNode getRoot() {
                return nodes[nodes.length-1];
            }



            private boolean subtreeFitsInNetwork(GeneUnionNode node,
                                                 final AlloppSpeciesNetworkModel asnm) {
                for (int i = 0; i < node.child.length; i++) {
                    if (!subtreeFitsInNetwork(node.child[i], asnm)) {
                        return false;
                    }
                }
                return asnm.coalescenceIsCompatible(node.height, node.union);
            }


            private void subtreeRecordCoalescences(GeneUnionNode node,
                                                   final AlloppSpeciesNetworkModel asnm) {
                for (int i = 0; i < node.child.length; i++) {
                    subtreeRecordCoalescences(node.child[i], asnm);
                }
                if (node.child.length > 0) {
                    asnm.recordCoalescence(node.height, node.union);
                }
            }


            /*
            * Recursively copies the topology from subtree rooted at node into
            * GeneUnionTree implemented as array nodes[].
            * Fills in union fields.
            *
            */
            private void genetree2geneuniontree(NodeRef gnode) {
                if (tree.isExternal(gnode)) {
                    nodes[nextn].child = new GeneUnionNode[0];
                    int ti = taxon2index.get(tree.getNodeTaxon(gnode));
                    int spseq = spsq[seqassigns[ti].spIndex][seqassigns[ti].seqIndex];
                    nodes[nextn].union.set(spseq);
                    nodes[nextn].name = tree.getNodeTaxon(gnode).getId();
                } else {
                    genetree2geneuniontree(tree.getChild(gnode,0));
                    int c0 = nextn - 1;
                    genetree2geneuniontree(tree.getChild(gnode,1));
                    int c1 = nextn - 1;
                    nodes[nextn].child = new GeneUnionNode[2];
                    nodes[nextn].child[0] = nodes[c0];
                    nodes[nextn].child[1] = nodes[c1];
                    nodes[nextn].union.union(nodes[c0].union);
                    nodes[nextn].union.union(nodes[c1].union);
                }
                nodes[nextn].height = tree.getNodeHeight(gnode);
                nextn++;
            }





            public String asText() {
                String s = "";
                Stack<Integer> x = new Stack<Integer>();
                return subtreeAsText(getRoot(), s, x, 0, "");
            }


            private String subtreeAsText(GeneUnionNode node, String s, Stack<Integer> x, int depth, String b) {
                Integer[] y = x.toArray(new Integer[x.size()]);
                StringBuffer indent = new StringBuffer();
                for (int i = 0; i < depth; i++) {
                    indent.append("  ");
                }
                for (int i = 0; i < y.length; i++) {
                    indent.replace(2*y[i], 2*y[i]+1, "|");
                }
                if (b.length() > 0) {
                    indent.replace(indent.length()-b.length(), indent.length(), b);
                }
                s += indent;
                s += node.asText(indent.length());
                s += System.getProperty("line.separator");
                String subs = "";
                if (node.child.length > 0) {
                    x.push(depth);
                    subs += subtreeAsText(node.child[0], "", x, depth+1, "-");
                    x.pop();
                    subs += subtreeAsText(node.child[1], "", x, depth+1, "`-");
                }
                return s + subs;
            }

        } // end GeneTreeInfo.GeneUnionTree


        /*
         * GeneTreeInfo constructor
         *
         * JH's SpeciesBindings code has test
         *   if (tree.getTaxonIndex(t) >= 0) { add taxon to count }
         *
         * I am not clear about what happens if some gene trees don't have all taxa.
         * I insist all gene trees have all taxa, and use missing data
         * grjtodo-oneday make more efficient
         */
        GeneTreeInfo(TreeModel tree, double popFactor, boolean permuteSequenceAssignments) {
            this.tree = tree;
            this.popFactor = popFactor;
            seqassigns = new SequenceAssignment[taxa.length];
            oldseqassigns = new SequenceAssignment[taxa.length];

            // This uses taxa list for *all* gene trees, not this gene tree.
            for (int s = 0; s < apspecies.length; s++) {
                for (int i = 0; i < apspecies[s].individuals.length; i++) {
                    int nseqs = apspecies[s].individuals[i].taxa.length;
                    int asgns[] = new int [nseqs];
                    for (int x = 0; x < nseqs; x++) {
                        asgns[x] = x;
                    }
                    if (permuteSequenceAssignments) { MathUtils.permute(asgns); }
                    for (int x = 0; x < nseqs; x++) {
                        int t = taxon2index.get(apspecies[s].individuals[i].taxa[x]);
                        seqassigns[t] = new SequenceAssignment(s, asgns[x]);
                        oldseqassigns[t] = new SequenceAssignment(s, asgns[x]);
                    }
                }
            }

            lineagesCount = new int[apspecies.length];
            Arrays.fill(lineagesCount, 0);

            for (int nl = 0; nl < lineagesCount.length; ++nl) {
                for (Individual indiv : apspecies[nl].individuals) {
                    boolean got = false;
                    for (Taxon t : indiv.taxa) {
                        if (tree.getTaxonIndex(t) >= 0) {
                            got = true;
                        }
                    }
                    for (Taxon t : indiv.taxa) {
                        assert (tree.getTaxonIndex(t) >= 0) == got;
                    }
                    assert got;
                    if (got) {
                        ++lineagesCount[nl];
                    }
                }
            }
        }


        public String seqassignsAsText() {
            String s = "Sequence assignments" + System.getProperty("line.separator");
            for (int tx = 0; tx < seqassigns.length; tx++) {
                s += taxa[tx];
                s += ":";
                s += seqassigns[tx].seqIndex;
                if (tx+1 < seqassigns.length  &&  seqassigns[tx].spIndex != seqassigns[tx+1].spIndex) {
                    s += System.getProperty("line.separator");
                } else {
                    s += "  ";
                }
            }
            return s;
        }



        public String genetreeAsText() {
            GeneUnionTree gutree = new GeneUnionTree();
            return gutree.asText();
        }


        public boolean fitsInNetwork(final AlloppSpeciesNetworkModel asnm) {
            GeneUnionTree gutree = new GeneUnionTree();
            boolean fits = gutree.subtreeFitsInNetwork(gutree.getRoot(), asnm);
            if (AlloppSpeciesNetworkModel.DBUGTUNE) {
                if (!fits) {
                    System.err.println("INCOMPATIBLE");
                    System.err.println(seqassignsAsText());
                    System.err.println(gutree.asText());
                    System.err.println(asnm.mullabTreeAsText());
                }
            }
            return fits;
        }


        // returns log(P(g_i|S)) = probability that gene tree fits into species network
        public double treeLogLikelihood(final AlloppSpeciesNetworkModel asnm) {
            GeneUnionTree gutree = new GeneUnionTree();
            asnm.clearCoalescences();
            gutree.subtreeRecordCoalescences(gutree.getRoot(), asnm);
            asnm.sortCoalescences();
            asnm.recordLineageCounts();
            double llhood = asnm.geneTreeInNetworkLogLikelihood();
            if (AlloppSpeciesNetworkModel.DBUGTUNE) {
                System.err.println("COMPATIBLE: log-likelihood = " + llhood);
                System.err.println(seqassignsAsText());
                System.err.println(gutree.asText());
                System.err.println(asnm.mullabTreeAsText());
            }
            return llhood;
        }


        public void storeSequenceAssignments() {
            for (int i = 0; i < seqassigns.length; i++) {
                oldseqassigns[i].seqIndex = seqassigns[i].seqIndex;
            }
        }

        public void restoreSequenceAssignments() {
            for (int i = 0; i < seqassigns.length; i++) {
                seqassigns[i].seqIndex = oldseqassigns[i].seqIndex;
            }
        }




        public double spseqUpperBound(FixedBitSet spsq0, FixedBitSet spsq1) {
            GeneUnionTree gutree = new GeneUnionTree();
            return subtreeSpseqUpperBound(gutree.getRoot(), spsq0, spsq1, Double.MAX_VALUE);
        }



        public void permuteOneSpeciesOneIndiv() {
            int sp = MathUtils.nextInt(apspecies.length);
            int iv = MathUtils.nextInt(apspecies[sp].individuals.length);
            flipOneAssignment(sp, iv);
        }


        /* grjtodo-oneday
         * This is a bit odd. It collects individuals as (sp, iv) indices
         * that `belong' to a node in the sense that any taxon (sequence)
         * of an individual belongs to the clade of the node.
         * I've used a set but not made SpeciesIndivPair's comparable
         * so that if both sequences of an individual occurs in clade it appears
         * twice. Then flipOneAssignment() flips everything so that those
         * occurring twice get flipped twice and so not changed.
         *
         * Result is that individuals with one but not two sequences in
         * the clade of the node get flipped. Sometimes all individuals
         * are flipped, sometimes none, sometimes just one, the last is the
         * same as permuteOneSpeciesOneIndiv().
         *
         * 2011-07-29 it appears to work OK on minimal testing and I
         * don't have a good idea for a more rational or efficient version.
         *
         */
        public void permuteSetOfIndivs() {
            int num = tree.getInternalNodeCount();
            int i = MathUtils.nextInt(num);
            NodeRef node = tree.getInternalNode(i);
            Set<SpeciesIndivPair> spivs = new HashSet<SpeciesIndivPair>();
            collectIndivsOfNode(node, spivs);
            for (SpeciesIndivPair spiv : spivs) {
                flipOneAssignment(spiv.spIndex, spiv.ivIndex);
            }
        }




        public SequenceAssignment getSeqassigns(int tx) {
            return seqassigns[tx];
        }


        // called when a gene tree has changed, which affects likelihood.
        // 2011-08-12 I am not using dirty flags (yet). I return
        // false from getLikelihoodKnown() in AlloppMSCoalescent
        // and that seems to be sufficient.
        public void wasChanged() {
        }






        private void collectIndivsOfNode(NodeRef node, Set<SpeciesIndivPair> spivs) {
            if (tree.isExternal(node)) {
                SpeciesIndivPair x = apspeciesId2speciesindiv(tree.getNodeTaxon(node).getId());
                spivs.add(x);
            } else {
                collectIndivsOfNode(tree.getChild(node, 0), spivs);
                collectIndivsOfNode(tree.getChild(node, 1), spivs);
            }
        }





        // start at root of gutree and recurse.
        // A node which has one child which contains some of species spp0
        // and where the other contains some of species spp1, imposes a limit
        // on how early a speciation can occur.
        private double subtreeSpseqUpperBound(GeneUnionNode node,
                                              FixedBitSet spsq0, FixedBitSet spsq1, double bound) {
            if (node.child.length == 0) {
                return bound;
            }
            for (GeneUnionNode ch : node.child) {
                bound = Math.min(bound, subtreeSpseqUpperBound(ch, spsq0, spsq1, bound));
            }
            FixedBitSet genespp0 = node.child[0].union;
            int int00 = genespp0.intersectCardinality(spsq0);
            int int01 = genespp0.intersectCardinality(spsq1);
            FixedBitSet genespp1 = node.child[1].union;
            int int10 = genespp1.intersectCardinality(spsq0);
            int int11 = genespp1.intersectCardinality(spsq1);
            if ((int00 > 0 && int11 > 0)  ||  (int10 > 0 && int01 > 0)) {
                bound = Math.min(bound, node.height);
            }
            return bound;
        }



        private void flipOneAssignment(int sp, int iv) {
            // grjtodo-tetraonly
            int tx;
            if (apspecies[sp].individuals[iv].taxa.length == 2) {
                tx = taxon2index.get(apspecies[sp].individuals[iv].taxa[0]);
                seqassigns[tx].seqIndex = 1 - seqassigns[tx].seqIndex;
                tx = taxon2index.get(apspecies[sp].individuals[iv].taxa[1]);
                seqassigns[tx].seqIndex = 1 - seqassigns[tx].seqIndex;
            }
        }



        private void flipAssignmentsForSpecies(int sp) {
            for (int iv = 0; iv < apspecies[sp].individuals.length; iv++) {
                flipOneAssignment(sp, iv);
            }
        }



    }
    // end of GeneTreeInfo




    /* ******************** Constructor *****************************/

    /*
     * The standard constructor calls this with permuteSequenceAssignments==true.
     * permuteSequenceAssignments==false is for testing.
     */
    public AlloppSpeciesBindings(ApSpInfo[] apspecies, TreeModel[] geneTrees,
                                 double minheight, double[] popFactors, boolean permuteSequenceAssignments) {
        super(AlloppSpeciesBindingsParser.ALLOPPSPECIES);

        this.apspecies = apspecies;
        initialmingenenodeheight = minheight;
        // make the flattened arrays
        int n = 0;
        for (int s = 0; s < apspecies.length; s++) {
            n += apspecies[s].individuals.length;
        }
        Individual [] indivs = new Individual[n];
        n = 0;
        for (int s = 0; s < apspecies.length; s++) {
            for (int i = 0; i < apspecies[s].individuals.length; i++, n++) {
                indivs[n] =  apspecies[s].individuals[i];
            }
        }
        int t = 0;
        for (int i = 0; i < indivs.length; i++) {
            t += indivs[i].taxa.length;
        }
        taxa = new Taxon[t];
        t = 0;
        for (int i = 0; i < indivs.length; i++) {
            for (int j = 0; j < indivs[i].taxa.length; j++, t++) {
                taxa[t] =  indivs[i].taxa[j];
            }
        }
        // set up maps to indices
        for (int i = 0; i < taxa.length; i++) {
            taxon2index.put(taxa[i], i);
        }
        spsq = new int[apspecies.length][];
        int spsqindex = 0;
        for (int sp = 0; sp < apspecies.length; sp++) {
            spsq[sp] = new int[apspecies[sp].ploidylevel/2];
            for (int seq = 0; seq < spsq[sp].length; seq++, spsqindex++) {
                spsq[sp][seq] = spsqindex;
            }
        }
        numberOfSpSeqs = spsqindex;

        geneTreeInfos = new GeneTreeInfo[geneTrees.length];
        for (int i = 0; i < geneTrees.length; i++) {
            geneTreeInfos[i] = new GeneTreeInfo(geneTrees[i], popFactors[i], permuteSequenceAssignments);
        }

        for (GeneTreeInfo gti : geneTreeInfos) {
            NodeRef[] nodes = gti.tree.getNodes();
            for (NodeRef node : nodes) {
                if (!gti.tree.isExternal(node)) {
                    double height = gti.tree.getNodeHeight(node);
                    gti.tree.setNodeHeight(node, minheight+height);
                }
            }
        }
    }


    /*
     * The normal constructor
     */
    public AlloppSpeciesBindings(ApSpInfo[] apspecies, TreeModel[] geneTrees,
                                 double minheight, double[] popFactors) {
        this(apspecies, geneTrees,  minheight, popFactors, true);
    }

    /*
     * Minimal constructor for testing conversions network -> multree, diphist
     * and for testing likelihood in MUL tree.
     */
    public AlloppSpeciesBindings(ApSpInfo[] apspecies) {
        this(apspecies, new TreeModel[0], 0.0, new double[0]);
    }


    public double initialMinGeneNodeHeight() {
        return initialmingenenodeheight;
    }



    public FixedBitSet speciesseqEmptyUnion() {
        FixedBitSet union = new FixedBitSet(numberOfSpSeqs());
        return union;
    }


    // Taxons vs species.
    // Taxons may have a final "0", "1",... to distinguish sequences, while
    // species do not. AlloppLeggedTree uses a SimpleTree, which only has
    // Taxons, so same thing there. Multree needs distinguishable Taxons
    // so has suffices.
    public FixedBitSet taxonseqToTipUnion(Taxon tx, int seq) {
        FixedBitSet union = speciesseqEmptyUnion();
        int sp = apspeciesId2index(tx.getId());
        int spseq = spandseq2spseqindex(sp, seq);
        union.set(spseq);
        return union;
    }


    public FixedBitSet spsqunion2spunion(FixedBitSet spsqunion) {
        FixedBitSet spunion = new FixedBitSet(apspecies.length);
        for (int sp = 0; sp < apspecies.length; sp++) {
            boolean got = false;
            for (int seq = 0; seq < spsq[sp].length; seq++) {
                if (spsqunion.contains(spsq[sp][seq])) {
                    got = true;
                }
            }
            if (got) {
                spunion.set(sp);
            }
        }
        return spunion;
    }


    public int numberOfGeneTrees() {
        return geneTreeInfos.length;
    }


    public double maxGeneTreeHeight() {
        if (geneTreeInfos.length == 0) {
            return 999;   // for test code only
        }
        double maxheight = 0.0;
        for (GeneTreeInfo gti : geneTreeInfos) {
            double height = gti.tree.getNodeHeight(gti.tree.getRoot());
            if (height > maxheight) {
                maxheight = height;
            }
        }
        return maxheight;
    }




    public boolean geneTreeFitsInNetwork(int i, final AlloppSpeciesNetworkModel asnm) {
        return geneTreeInfos[i].fitsInNetwork(asnm);
    }

    public double geneTreeLogLikelihood(int i, final AlloppSpeciesNetworkModel asnm) {
        return geneTreeInfos[i].treeLogLikelihood(asnm);
    }


    public int numberOfSpecies() {
        return apspecies.length;
    }



    public String apspeciesName(int i) {
        return apspecies[i].name;
    }



    public Taxon[] SpeciesWithinPloidyLevel(int pl) {
        ArrayList<Taxon> names = new ArrayList<Taxon>();
        for (int i = 0; i < apspecies.length; i++) {
            if (apspecies[i].ploidylevel == pl) {
                names.add(new Taxon(apspecies[i].name));
            }
        }
        Taxon[] spp = new Taxon[names.size()];
        names.toArray(spp);
        return spp;
    }


    public int spandseq2spseqindex(int sp, int seq) {
        return spsq[sp][seq];
    }



    public int spseqindex2sp(int spsqindex) {
        return spseqindex2spandseq(spsqindex)[0];
    }

    public int spseqindex2seq(int spsqindex) {
        return spseqindex2spandseq(spsqindex)[1];
    }



    public int apspeciesId2index(String apspId) {
        int index = -1;
        for (int i = 0; i < apspecies.length; i++) {
            if (apspecies[i].name.compareTo(apspId) == 0) {
                assert index == -1;
                index = i;
            }
        }
        if (index == -1) {
            System.out.println("BUG in apspeciesId2index");
        }
        assert index != -1;
        return index;
    }


    public SpeciesIndivPair apspeciesId2speciesindiv(String apspId) {
        int sp = -1;
        int iv = -1;
        for (int s = 0; s < apspecies.length; s++) {
            for (int i = 0; i < apspecies[s].individuals.length; i++) {
                for (int t = 0; t < apspecies[s].individuals[i].taxa.length; t++) {
                    Taxon taxon = apspecies[s].individuals[i].taxa[t];
                    if (taxon.getId().compareTo(apspId) == 0) {
                        sp = s;
                        iv = i;
                    }
                }
            }
        }
        assert sp != -1;
        SpeciesIndivPair x = new SpeciesIndivPair(sp, iv);

        return x;
    }


    public int numberOfSpSeqs() {
        return numberOfSpSeqs;
    }


    int nLineages(int speciesIndex) {
        int n = geneTreeInfos[0].lineagesCount[speciesIndex];
        for (GeneTreeInfo gti : geneTreeInfos) {
            assert gti.lineagesCount[speciesIndex] == n;
        }
        return n;
    }



    public double spseqUpperBound(FixedBitSet left, FixedBitSet right) {
        double bound = Double.MAX_VALUE;
        for (GeneTreeInfo gti : geneTreeInfos) {
            bound = Math.min(bound, gti.spseqUpperBound(left, right));
        }
        return bound;
    }



    public void permuteOneSpeciesOneIndivForOneGene() {
        int i = MathUtils.nextInt(geneTreeInfos.length);
        geneTreeInfos[i].permuteOneSpeciesOneIndiv();
    }


    public void permuteSetOfIndivsForOneGene() {
        int i = MathUtils.nextInt(geneTreeInfos.length);
        geneTreeInfos[i].permuteSetOfIndivs();
    }


    public void flipAssignmentsForAllGenesOneSpecies(int sp) {
        for (GeneTreeInfo gti : geneTreeInfos) {
            gti.flipAssignmentsForSpecies(sp);
        }

    }



    public String seqassignsAsText(int g) {
        return geneTreeInfos[g].seqassignsAsText();
    }


    public String genetreeAsText(int g) {
        String s = "Gene tree " + g + "                     height             union" + System.getProperty("line.separator");
        s += geneTreeInfos[g].genetreeAsText();
        return s;
    }


    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
        for (GeneTreeInfo g : geneTreeInfos) {
            if (g.tree == model) {
                g.wasChanged();
                break;
            }
        }
        fireModelChanged(object, index);
        if (AlloppSpeciesNetworkModel.DBUGTUNE)
            System.err.println("AlloppSpeciesBindings.handleModelChangedEvent() " + model.getId());
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index,
                                              ChangeType type) {
        assert false; // copies SpeciesBindings; not understood
        if (AlloppSpeciesNetworkModel.DBUGTUNE)
            System.err.println("AlloppSpeciesBindings.handleVariableChangedEvent() " + variable.getId());
    }

    @Override
    protected void storeState() {
        for (GeneTreeInfo gti : geneTreeInfos) {
            gti.storeSequenceAssignments();
        }
        if (AlloppSpeciesNetworkModel.DBUGTUNE)
            System.err.println("AlloppSpeciesBindings.storeState()");
    }

    @Override
    protected void restoreState() {
        for (GeneTreeInfo gti : geneTreeInfos) {
            gti.restoreSequenceAssignments();
            if (AlloppSpeciesNetworkModel.DBUGTUNE)
                System.err.println("AlloppSpeciesBindings.restoreState()");
        }

    }

    @Override
    protected void acceptState() {
    }


    public void addModelListeners(ModelListener listener) {
        for (GeneTreeInfo gti : geneTreeInfos) {
            gti.tree.addModelListener(listener);
        }
        addModelListener(listener); // for sequence assignments
    }



    public LogColumn[] getColumns() {
        int ncols = geneTreeInfos.length * taxa.length;
        LogColumn[] columns = new LogColumn[ncols];
        for (int g = 0, i = 0; g < geneTreeInfos.length; g++) {
            for (int tx = 0; tx < taxa.length; tx++, i++) {
                GeneTreeInfo.SequenceAssignment sqa = geneTreeInfos[g].getSeqassigns(tx);
                String header = "Gene" + g + "taxon" + tx;
                columns[i] = new LogColumn.Default(header, sqa);
            }

        }

        return columns;
    }



    private int[] spseqindex2spandseq(int spsqindex) {
        int indexp = -1;
        int indexq = -1;
        for (int p = 0; p < spsq.length; p++) {
            for (int q = 0; q < spsq[p].length; q++) {
                if (spsq[p][q] == spsqindex) {
                    assert indexp == -1;
                    assert indexq == -1;
                    indexp = p;
                    indexq = q;
                }
            }
        }
        assert indexp != -1;
        assert indexq != -1;
        int[] pq = new int[2];
        pq[0] = indexp;
        pq[1] = indexq;
        return pq;
    }

}



