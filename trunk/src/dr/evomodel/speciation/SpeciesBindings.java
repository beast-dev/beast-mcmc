package dr.evomodel.speciation;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxon;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.ParameterChangeType;
import dr.util.HeapSort;
import dr.xml.*;
import jebl.util.FixedBitSet;

import java.util.*;

/**
 * Binds taxa in gene trees with species information.
 *
 * @author joseph
 *         Date: 25/05/2008
 */
public class SpeciesBindings extends AbstractModel {
    public static final String SPECIES = "species";
    public static final String SP = "sp";
    public static final String GENE_TREES = "geneTrees";

    // all gene trees
    GeneTreeInfo[] geneTrees;

    // convenience
    Map<Taxon, Integer> taxon2Species = new HashMap<Taxon, Integer>();

    // Species definition
    SPinfo[] species;

    private double[][] popTimesPair;
    private boolean dirty_pp;

    private double[][] popTimesSingle;
    private boolean dirty_sg;
    private boolean verbose = false;

    SpeciesBindings(SPinfo[] species, TreeModel[] geneTrees) {
        super(null);

        this.species = species;

        final int nsp = species.length;

        for (int ns = 0; ns < nsp; ++ns) {
            for (Taxon t : species[ns].taxa) {
                if (taxon2Species.containsKey(t)) {
                    throw new Error("Multiple assignments for taxon" + t);
                }
                taxon2Species.put(t, ns);
            }
        }

        this.geneTrees = new GeneTreeInfo[geneTrees.length];

        for (int i = 0; i < geneTrees.length; i++) {
            final TreeModel t = geneTrees[i];
            addModel(t);
            this.geneTrees[i] = new GeneTreeInfo(t);
        }

        for (GeneTreeInfo gt : this.geneTrees) {
            for (int ns = 0; ns < nsp; ++ns) {
                if (gt.nLineages(ns) == 0) {
                    throw new Error("Every gene tree must contain at least one tip from each species");
                }
            }
        }

        popTimesSingle = new double[nsp][];
        for (int ns = 0; ns < popTimesSingle.length; ++ns) {
            popTimesSingle[ns] = new double[allCoalPointsCount(ns)];
        }
        dirty_sg = true;

        popTimesPair = new double[(nsp * (nsp - 1)) / 2][];
        {
            final int nps = allPairCoalPointsCount();
            for (int ns = 0; ns < popTimesPair.length; ++ns) {
                popTimesPair[ns] = new double[nps];
            }
        }

        dirty_pp = true;
    }

    public int nSpecies() {
        return species.length;
    }

    /**
     * Per species coalecent times.
     * <p/>
     * Indexed by sp index, a list of coalescent times of taxa of this sp from all gene trees.
     *
     * @return Per species coalecent times
     */
    public double[][] getPopTimesSingle() {
        if (dirty_sg) {
            for (int ns = 0; ns < popTimesSingle.length; ++ns) {
                getAllCoalPoints(ns, popTimesSingle[ns]);
            }
            dirty_sg = false;
        }
        return popTimesSingle;
    }

    public double[][] getPopTimesPair() {
        if (dirty_pp) {
            final int nsp = nSpecies();
            for (int ns1 = 0; ns1 < nsp - 1; ++ns1) {
                final int z = (ns1 * (2 * nsp - ns1 - 3)) / 2 - 1;

                for (int ns2 = ns1 + 1; ns2 < nsp; ++ns2) {
                    getAllPairCoalPoints(ns1, ns2, popTimesPair[z + ns2]);
                }
            }
        }
        return popTimesPair;
    }

    private void getAllPairCoalPoints(int ns1, int ns2, double[] popTimes) {

        for (int i = 0; i < geneTrees.length; i++) {
            for (CoalInfo ci : geneTrees[i].getCoalInfo()) {
                if ((ci.sinfo[0].contains(ns1) && ci.sinfo[1].contains(ns2)) ||
                        (ci.sinfo[1].contains(ns1) && ci.sinfo[0].contains(ns2))) {
                    popTimes[i] = ci.ctime;
                    break;
                }
            }
        }
        HeapSort.sort(popTimes);
    }

    private int allCoalPointsCount(int spIndex) {
        int tot = 0;
        for (GeneTreeInfo t : geneTrees) {
            if (t.nLineages(spIndex) > 0) {
                tot += t.nLineages(spIndex) - 1;
            }
        }
        return tot;
    }

    void getAllCoalPoints(int spIndex, double[] points) {

        int k = 0;
        for (GeneTreeInfo t : geneTrees) {
            int k1 = t.nLineages(spIndex) - 1;
            int savek = k;
            for (CoalInfo ci : t.getCoalInfo()) {
//               if( ci == null ) {
//                assert ci != null;
//            }
                if (ci.allHas(spIndex)) {
                    points[k] = ci.ctime;
                    ++k;
                }
            }
            if (!(k1 >= 0 && savek + k1 == k) || (k1 < 0 && savek == k)) {
                System.err.println(k1);
            }
            assert (k1 >= 0 && savek + k1 == k) || (k1 < 0 && savek == k);
        }
        assert k == points.length;
        HeapSort.sort(points);
    }

    private int allPairCoalPointsCount() {
        return geneTrees.length;
    }

    public double speciationUpperBound(FixedBitSet sub1, FixedBitSet sub2) {
        //Determined by the last time any pair of sp's in sub1 x sub2 have been seen
        // together in any of the gene trees."""

        double bound = Double.MAX_VALUE;
        for (GeneTreeInfo g : getGeneTrees()) {
            for (CoalInfo ci : g.getCoalInfo()) {
                // if past time of current bound, can't change it anymore
                if (ci.ctime >= bound) {
                    break;
                }
                if ((ci.sinfo[0].intersectCardinality(sub1) > 0 && ci.sinfo[1].intersectCardinality(sub2) > 0)
                        ||
                        (ci.sinfo[0].intersectCardinality(sub2) > 0 && ci.sinfo[1].intersectCardinality(sub1) > 0)) {
                    bound = ci.ctime;
                    break;
                }
            }
        }
        return bound;
    }

    /**
     * Information on one species (sp)
     */
    static class SPinfo extends Taxon {
        // sp name
        final public String name;

        // all taxa belonging to sp
        private Taxon[] taxa;

        SPinfo(String name, Taxon[] taxa) {
            super(name);

            this.name = name;
            this.taxa = taxa;
        }
    }

    class CoalInfo implements Comparable<CoalInfo> {
        // zero based, 0 is taxa time, i.e. in tree branch units
        final double ctime;
        // sp info for each subtree
        FixedBitSet[] sinfo;

        CoalInfo(double t, int nc) {
            ctime = t;
            sinfo = new FixedBitSet[nc];
        }

        public int compareTo(CoalInfo o) {
            return o.ctime < ctime ? +1 : (o.ctime > ctime ? -1 : 0);
        }

        /**
         * @param s
         * @return true if all children have at least one taxa from sp 's'
         */
        public boolean allHas(int s) {
            for (FixedBitSet b : sinfo) {
                if (!b.contains(s)) {
                    return false;
                }
            }
            return true;
        }
    }

    /**
     * Collect coalescence information for sub-tree rooted at 'node'.
     *
     * @param tree
     * @param node
     * @param loc  Place node data in loc, sub-tree info before that.
     * @param info array to fill
     * @return location of next available location
     */
    private int collectCoalInfo(Tree tree, NodeRef node, int loc, CoalInfo[] info) {

        info[loc] = new CoalInfo(tree.getNodeHeight(node), tree.getChildCount(node));

        int newLoc = loc - 1;
        for (int i = 0; i < 2; i++) {
            NodeRef child = tree.getChild(node, i);
            info[loc].sinfo[i] = new FixedBitSet(nSpecies());

            if (tree.isExternal(child)) {
                info[loc].sinfo[i].set(taxon2Species.get(tree.getNodeTaxon(child)));
                assert tree.getNodeHeight(child) == 0;
            } else {
                final int used = collectCoalInfo(tree, child, newLoc, info);
                for (int j = 0; j < info[newLoc].sinfo.length; ++j) {
                    info[loc].sinfo[i].union(info[newLoc].sinfo[j]);
                }
                newLoc = used;
            }
        }
        return newLoc;
    }

    public class GeneTreeInfo {
        public final TreeModel tree;
        private int[] lineagesCount;
        private CoalInfo[] cList;
        private CoalInfo[] savedcList;
        private boolean dirty;
        private boolean wasBacked;

        GeneTreeInfo(TreeModel tree) {
            this.tree = tree;

            lineagesCount = new int[species.length];
            Arrays.fill(lineagesCount, 0);

            for (int nl = 0; nl < lineagesCount.length; ++nl) {
                for (Taxon t : species[nl].taxa) {
                    if (tree.getTaxonIndex(t) >= 0) {
                        ++lineagesCount[nl];
                    }
                }
            }

            cList = new CoalInfo[tree.getExternalNodeCount() - 1];
            savedcList = new CoalInfo[cList.length];
            wasChanged();
            getCoalInfo();
            wasBacked = false;
        }

        int nLineages(int speciesIndex) {
            return lineagesCount[speciesIndex];
        }

        public CoalInfo[] getCoalInfo() {
            if (dirty) {
                swap();

                collectCoalInfo(tree, tree.getRoot(), cList.length - 1, cList);
                HeapSort.sort(cList);
                dirty = false;
                wasBacked = true;
            }
            return cList;
        }

        private void swap() {
            CoalInfo[] tmp = cList;
            cList = savedcList;
            savedcList = tmp;
        }

        void wasChanged() {
            dirty = true;
            wasBacked = false;
        }

        boolean restore() {
            if (verbose) System.out.println(" SP binding: restore " + tree.getId() + " (" + wasBacked + ")");
            if (wasBacked) {
//                if( false ) {
//                    swap();
//                    dirty = true;
//                    getCoalInfo();
//                    for(int k = 0; k < cList.length; ++k) {
//                        assert cList[k].ctime == savedcList[k].ctime &&
//                                cList[k].sinfo[0].equals(savedcList[k].sinfo[0]) &&
//                                cList[k].sinfo[1].equals(savedcList[k].sinfo[1]);
//                    }
//                }
                swap();
                wasBacked = false;
                dirty = false;
                return true;
            }
            return false;
        }

        void accept() {
            if (verbose) System.out.println(" SP binding: accept " + tree.getId());

            wasBacked = false;
        }
    }

    public GeneTreeInfo[] getGeneTrees() {
        return geneTrees;
    }

    protected void handleModelChangedEvent(Model model, Object object, int index) {
        if (verbose) System.out.println(" SP binding: model changed " + model.getId());

        dirty_sg = true;
        dirty_pp = true;

        for (GeneTreeInfo g : geneTrees) {
            if (g.tree == model) {
                g.wasChanged();
                break;
            }
        }
        fireModelChanged(object, index);
    }

    protected final void handleParameterChangedEvent(Parameter parameter, int index, ParameterChangeType type) {
        assert false;
    }

    protected void storeState() {
        // do on a per need basis
    }

    protected void restoreState() {
        for (GeneTreeInfo g : geneTrees) {
            if (g.restore()) {
                dirty_sg = true;
                dirty_pp = true;
            }
        }
    }

    protected void acceptState() {
        for (GeneTreeInfo g : geneTrees) {
            g.accept();
        }
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return SPECIES;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            List<SPinfo> sp = new ArrayList<SPinfo>();
            for (int k = 0; k < xo.getChildCount(); ++k) {
                final Object child = xo.getChild(k);
                if (child instanceof SPinfo) {
                    sp.add((SPinfo) child);
                }
            }

            XMLObject xogt = (XMLObject) xo.getChild(GENE_TREES);
            TreeModel[] trees = new TreeModel[xogt.getChildCount()];
            for (int nt = 0; nt < trees.length; ++nt) {
                trees[nt] = (TreeModel) xogt.getChild(nt);
            }

            try {
                return new SpeciesBindings(sp.toArray(new SPinfo[sp.size()]), trees);
            } catch (Error e) {
                throw new XMLParseException(e.getMessage());
            }
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return new XMLSyntaxRule[]{
                    new ElementRule(SPinfo.class, 2, Integer.MAX_VALUE),
                    new ElementRule(GENE_TREES,
                            new XMLSyntaxRule[]{new ElementRule(TreeModel.class, 1, Integer.MAX_VALUE)}),
            };
        }

        public String getParserDescription() {
            return "Binds taxa in gene trees with species information.";
        }

        public Class getReturnType() {
            return SpeciesBindings.class;
        }
    };

    public static XMLObjectParser PPARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return SP;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            Taxon[] taxa = new Taxon[xo.getChildCount()];
            for (int nt = 0; nt < taxa.length; ++nt) {
                taxa[nt] = (Taxon) xo.getChild(nt);
            }
            return new SPinfo(xo.getId(), taxa);
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return new XMLSyntaxRule[]{
                    new ElementRule(Taxon.class, 1, Integer.MAX_VALUE)
            };
        }

        public String getParserDescription() {
            return "Taxon in a species tree";
        }

        public Class getReturnType() {
            return SPinfo.class;
        }
    };
}