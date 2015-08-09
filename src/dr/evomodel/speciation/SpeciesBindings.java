/*
 * SpeciesBindings.java
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

package dr.evomodel.speciation;

import dr.evolution.tree.MutableTree;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxon;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.*;
import dr.util.HeapSort;
import jebl.util.FixedBitSet;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

/**
 * Binds taxa in gene trees with species information.
 *
 * @author Joseph Heled
 *         Date: 25/05/2008
 */
public class SpeciesBindings extends AbstractModel {
    // all gene trees
    private final GeneTreeInfo[] geneTrees;

    // convenience
    private final Map<Taxon, Integer> taxon2Species = new HashMap<Taxon, Integer>();

    // Species definition
    final SPinfo[] species;

    private final double[][] popTimesPair;
    private boolean dirty_pp;

    private final double[][] popTimesSingle;
    private boolean dirty_sg;
    private final boolean verbose = false;

    public SpeciesBindings(SPinfo[] species, TreeModel[] geneTrees, double[] popFactors) {
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
            this.geneTrees[i] = new GeneTreeInfo(t, popFactors[i]);
        }

        for (GeneTreeInfo gt : this.geneTrees) {
            for (int ns = 0; ns < nsp; ++ns) {
                if (gt.nLineages(ns) == 0) {
                    throw new Error("Every gene tree must contain at least one tip from each species ("
                            + gt.tree.getId() + "," + species[ns].name + ")");
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

        addStatistic(new SpeciesLimits());
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

    // length of points must be right
    void getAllCoalPoints(int spIndex, double[] points) {

        int k = 0;
        for (GeneTreeInfo t : geneTrees) {
            final int totCoalEvents = t.nLineages(spIndex) - 1;
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
            if (!(totCoalEvents >= 0 && savek + totCoalEvents == k) || (totCoalEvents < 0 && savek == k)) {
                System.err.println(totCoalEvents);
            }
            assert (totCoalEvents >= 0 && savek + totCoalEvents == k) || (totCoalEvents < 0 && savek == k);
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

    public void makeCompatible(double rootHeight) {
        for( GeneTreeInfo t : getGeneTrees() ) {

            MutableTree tree = t.tree;

            for (int i = 0; i < tree.getExternalNodeCount(); i++) {
                final NodeRef node = tree.getExternalNode(i);
                final NodeRef p = tree.getParent(node);
                tree.setNodeHeight(p, rootHeight + tree.getNodeHeight(p));
            }
            MutableTree.Utils.correctHeightsForTips(tree);
             // (todo) ugly re-init - can I do something better?
            t.wasChanged();
            t.getCoalInfo();
            t.wasBacked = false;
            //t.wasChanged();
       }
    }

    /**
     * Information on one species (sp)
     */
    public static class SPinfo extends Taxon {
        // sp name
        final public String name;

        // all taxa belonging to sp
        private final Taxon[] taxa;

        public SPinfo(String name, Taxon[] taxa) {
            super(name);

            this.name = name;
            this.taxa = taxa;
        }
    }

    class CoalInfo implements Comparable<CoalInfo> {
        // zero based, 0 is taxa time, i.e. in tree branch units
        final double ctime;
        // sp info for each subtree
        final FixedBitSet[] sinfo;

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
        private final int[] lineagesCount;
        private CoalInfo[] cList;
        private CoalInfo[] savedcList;
        private boolean dirty;
        private boolean wasBacked;
        private final double popFactor;

        GeneTreeInfo(TreeModel tree, double popFactor) {
            this.tree = tree;
            this.popFactor = popFactor;

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

        public double popFactor() {
            return popFactor;
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

    protected final void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
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

    public class SpeciesLimits extends Statistic.Abstract {
        int nDim;
        int c[][];

        SpeciesLimits() {
            super("SpeciationBounds");

            nDim = 0;

            final int nsp = species.length;

            c = new int[nsp + 1][nsp + 1];
            for(int k = 0; k < nsp + 1; ++k) {
                c[k][0] = 1;
                c[k][k] = 1;
            }
            for(int k = 0; k < nsp + 1; ++k) {
                for(int j = 1; j < k; ++j) {
                    c[k][j] = c[k - 1][j - 1] + c[k - 1][j];
                }
            }

            for(int k = 0; k <= (int) (nsp / 2); ++k) {
                nDim += c[nsp][k];
            }

        }

        public int getDimension() {
            return nDim;
        }

        private double boundOnRoot() {
            double bound = Double.MAX_VALUE;
            final int nsp = species.length;
            for(GeneTreeInfo g : getGeneTrees()) {
                for(CoalInfo ci : g.getCoalInfo()) {
                    if( ci.sinfo[0].cardinality() == nsp || ci.sinfo[1].cardinality() == nsp ) {
                        bound = Math.min(bound, ci.ctime);
                        break;
                    }
                }
            }
            return bound;
        }

        public double getStatisticValue(int dim) {
            if( dim == 0 ) {
                return boundOnRoot();
            }

            final int nsp = species.length;
            int r = 0;
            int k;
            for(k = 0; k <= (int) (nsp / 2); ++k) {
                final int i = c[nsp][k];
                if( dim < r + i ) {
                    break;
                }
                r += i;
            }

            // Classic index -> select k of nsp subset

            // number of species in set is k
            int n = dim - r;
            FixedBitSet in = new FixedBitSet(nsp),
                    out = new FixedBitSet(nsp);
            int fr = nsp;
            for(int i = 0; i < nsp; ++i) {
                if( k == 0 ) {
                    out.set(i);
                } else {
                    if( n < c[fr - 1][k - 1] ) {
                        in.set(i);
                        k -= 1;
                    } else {
                        out.set(i);
                        n -= c[fr - 1][k];
                    }
                    fr -= 1;
                }
            }
            return speciationUpperBound(in, out);
        }
    }
}