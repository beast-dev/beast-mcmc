package dr.evomodel.speciation;

import dr.evolution.tree.NodeRef;
import dr.evolution.util.Taxon;
import dr.evomodel.tree.TreeModel;
import dr.evomodelxml.speciation.PopsIOSpeciesBindingsParser;
import dr.inference.model.*;
import dr.util.AlloppMisc;
import jebl.util.FixedBitSet;

import java.util.*;

/**
 * User: Graham  Jones
 * Date: 10/05/12
 */



/*
 *
 * PopsIOSpeciesBindings knows how species are made of individuals (Taxons).
 * For a given gene an individual translates inso a sequence.
 *
 * It also contains the list of gene trees - tree topologies and node
 * times, plus popfactors. Given a PopsIOSpeciesTreeModel
 * it can say if a gene tree is compatible.
 *
 * The class GeneUnionTree, defined here, is used during calculations.
 */
public class PopsIOSpeciesBindings extends AbstractModel   {
    private final GeneTreeInfo[] geneTreeInfos;
    private final SpInfo[] spInfos;
    private final Taxon [] taxa;
    private final Map<Taxon, Integer> taxon2index = new HashMap<Taxon, Integer>();
    private final double initialmingenenodeheight; // for starting network



    /*
     * A SpInfo is a list of taxa = sequences
     */
    public static class SpInfo {
        final public String name;
        private final Taxon[] taxa;

        public SpInfo(String name, Taxon[] taxa) {
            this.name = name;
            this.taxa = taxa;
        }
    }


    public class GeneTreeInfo {
        public final TreeModel tree;
        private final int[] lineagesCount;
        private final double popFactor;

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
                union = new FixedBitSet(getSpecies().length);
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
                while (s.length() < 20-indentlen) {
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
                                                 final PopsIOSpeciesTreeModel piostm) {
                for (int i = 0; i < node.child.length; i++) {
                    if (!subtreeFitsInNetwork(node.child[i], piostm)) {
                        return false;
                    }
                }
                return piostm.coalescenceIsCompatible(node.height, node.union);
            }


            private void subtreeRecordCoalescences(GeneUnionNode node,
                                                   final PopsIOSpeciesTreeModel piostm) {
                for (int i = 0; i < node.child.length; i++) {
                    subtreeRecordCoalescences(node.child[i], piostm);
                }
                if (node.child.length > 0) {
                    piostm.recordCoalescence(node.height, node.union);
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
                    nodes[nextn].union.set(ti);
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
       * For now, 2012-05-10 I insist all gene trees have all taxa.
       */
        GeneTreeInfo(TreeModel tree, double popFactor) {
            this.tree = tree;
            this.popFactor = popFactor;

            lineagesCount = new int[spInfos.length];
            Arrays.fill(lineagesCount, 0);

            for (int nl = 0; nl < lineagesCount.length; ++nl) {
                for (Taxon tx : spInfos[nl].taxa) {
                        ++lineagesCount[nl];
                    }
                }
            }



        public String genetreeAsText() {
            GeneUnionTree gutree = new GeneUnionTree();
            return gutree.asText();
        }


        public boolean fitsInNetwork(final PopsIOSpeciesTreeModel piostm) {
            GeneUnionTree gutree = new GeneUnionTree();
            boolean fits = gutree.subtreeFitsInNetwork(gutree.getRoot(), piostm);
            return fits;
        }


        // returns log(P(g_i|S)) = probability that gene tree fits into species network
        public double treeLogLikelihood(final PopsIOSpeciesTreeModel piostm) {
            GeneUnionTree gutree = new GeneUnionTree();
            piostm.clearCoalescences();
            gutree.subtreeRecordCoalescences(gutree.getRoot(), piostm);
            piostm.sortCoalescences();
            piostm.recordLineageCounts();
            double llhood = piostm.geneTreeInSpeciesTreeLogLikelihood();
            return llhood;
        }

        public double coalescenceUpperBoundBetween(FixedBitSet spp0, FixedBitSet spp1) {
            GeneUnionTree gutree = new GeneUnionTree();
            return subtreeUpperBoundBetween(gutree.getRoot(), spp0, spp1, Double.MAX_VALUE);
        }


        // start at root of gutree and recurse.
        // A node which has one child which contains some of species spp0
        // and where the other contains some of species spp1, imposes a limit
        // on how early a speciation can occur.
        private double subtreeUpperBoundBetween(GeneUnionNode node,
                                              FixedBitSet spp0, FixedBitSet spp1, double bound) {
            if (node.child.length == 0) {
                return bound;
            }
            for (GeneUnionNode ch : node.child) {
                bound = Math.min(bound, subtreeUpperBoundBetween(ch, spp0, spp1, bound));
            }
            FixedBitSet genespp0 = node.child[0].union;
            int int00 = genespp0.intersectCardinality(spp0);
            int int01 = genespp0.intersectCardinality(spp1);
            FixedBitSet genespp1 = node.child[1].union;
            int int10 = genespp1.intersectCardinality(spp0);
            int int11 = genespp1.intersectCardinality(spp1);
            if ((int00 > 0 && int11 > 0)  ||  (int10 > 0 && int01 > 0)) {
                bound = Math.min(bound, node.height);
            }
            return bound;
        }


    }

    public PopsIOSpeciesBindings(SpInfo[] spInfos, TreeModel[] geneTrees,
                                 double minheight, double[] popFactors) {
        super(PopsIOSpeciesBindingsParser.PIO_SPECIES_BINDINGS);
        this.spInfos = spInfos;
        initialmingenenodeheight = minheight;

        int t = 0;
        for (SpInfo spi : spInfos) {
            t += spi.taxa.length;
        }
        taxa = new Taxon[t];
        t = 0;
        for (SpInfo spi : spInfos) {
            for (int j = 0; j < spi.taxa.length; j++, t++) {
                taxa[t] = spi.taxa[j];
            }
        }
        // set up maps to indices
        for (int i = 0; i < taxa.length; i++) {
            taxon2index.put(taxa[i], i);
        }

        geneTreeInfos = new GeneTreeInfo[geneTrees.length];
        for (int i = 0; i < geneTrees.length; i++) {
            geneTreeInfos[i] = new GeneTreeInfo(geneTrees[i], popFactors[i]);
        }

        for (GeneTreeInfo gti : geneTreeInfos) {
            NodeRef[] nodes = gti.tree.getNodes();
            for (NodeRef node : nodes) {
                if (!gti.tree.isExternal(node)) {
                    double height = gti.tree.getNodeHeight(node);
                    gti.tree.setNodeHeight(node, minheight + height);
                }
            }
        }

    }


    public int numberOfGeneTrees() {
        return geneTreeInfos.length;
    }


    public double maxGeneTreeHeight() {
        double maxheight = 0.0;
        for (GeneTreeInfo gti : geneTreeInfos) {
            double height = gti.tree.getNodeHeight(gti.tree.getRoot());
            if (height > maxheight) {
                maxheight = height;
            }
        }
        return maxheight;
    }



    public SpInfo [] getSpecies() {
        return spInfos;
    }


    public int speciesId2index(String spId) {
        int index = -1;
        for (int i = 0; i < spInfos.length; i++) {
            if (spInfos[i].name.compareTo(spId) == 0) {
                assert index == -1;
                index = i;
            }
        }
        if (index == -1) {
            System.out.println("BUG in speciesId2index");
        }
        assert index != -1;
        return index;
    }

    public double initialMinGeneNodeHeight() {
       return initialmingenenodeheight;
    }



    public boolean geneTreeFitsInNetwork(int i, PopsIOSpeciesTreeModel piostm) {
        return geneTreeInfos[i].fitsInNetwork(piostm);
    }

    public double geneTreeLogLikelihood(int i, PopsIOSpeciesTreeModel piostm) {
        return geneTreeInfos[i].treeLogLikelihood(piostm);
    }

    public FixedBitSet emptyUnion() {
        return new FixedBitSet(taxa.length);
    }


    public FixedBitSet tipUnionFromTaxon(Taxon tx) {
        int spi = speciesId2index(tx.getId());
        FixedBitSet x = new FixedBitSet(taxa.length);
        x.set(spi);
        return x;
    }


    public double coalescenceUpperBoundBetween(FixedBitSet left, FixedBitSet right) {
        double bound = Double.MAX_VALUE;
        for (GeneTreeInfo g : geneTreeInfos) {
            bound = Math.min(bound, g.coalescenceUpperBoundBetween(left, right));
        }
        return bound;
    }

    int nLineages(int speciesIndex) {
        int n = geneTreeInfos[0].lineagesCount[speciesIndex];
        for (GeneTreeInfo gti : geneTreeInfos) {
            assert gti.lineagesCount[speciesIndex] == n;
        }
        return n;
    }

    public GeneTreeInfo[] getGeneTrees() {
        return geneTreeInfos;
    }

    public String genetreeAsText(int g) {
        return geneTreeInfos[g].genetreeAsText();
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
        fireModelChanged(object, index);
        // grjtodo-oneday copied from elsewhere; not understood.
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Variable.ChangeType type) {
        assert false;
    }

    @Override
    protected void storeState() {
     }

    @Override
    protected void restoreState() {
    }

    @Override
    protected void acceptState() {
   }
}
