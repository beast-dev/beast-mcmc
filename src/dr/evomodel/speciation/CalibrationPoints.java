/*
 * CalibrationPoints.java
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

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeUtils;
import dr.evolution.util.Taxa;
import dr.inference.model.Statistic;
import dr.math.distributions.Distribution;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * @author Joseph Heled
 *         Date: 8/06/2011
 */
public class CalibrationPoints {

    public static enum CorrectionType {
        EXACT("exact"),
        APPROXIMATED("approximated"),
        PEXACT("pexact"),
        NONE("none");

        CorrectionType(String name) {
            this.name = name;
        }

        public String toString() {
            return name;
        }

        private final String name;
    }

    public CalibrationPoints(Tree tree, boolean isYule, List<Distribution> dists, List<Taxa> clades,
                             List<Boolean> forParent, Statistic userPDF, CorrectionType correctionType) {

        this.densities = new Distribution[dists.size()];
        this.clades = new int[clades.size()][];
        this.forParent = new boolean[clades.size()];

        for (int k = 0; k < clades.size(); ++k) {
            final Taxa tk = clades.get(k);
            for (int i = k + 1; i < clades.size(); ++i) {
                final Taxa ti = clades.get(i);
                if (ti.containsAny(tk)) {
                    if (!(ti.containsAll(tk) || tk.containsAll(ti))) {
                        throw new IllegalArgumentException("Overlapping clades??");
                    }
                }
            }
        }

        Taxa[] taxaInOrder = new Taxa[clades.size()];

        {
            int loc = clades.size() - 1;
            while (loc >= 0) {
                //  place maximal clades at end one at a time
                int k = 0;
                for (/**/; k < clades.size(); ++k) {
                    if (isMaximal(clades, k)) {
                        break;
                    }
                }
                this.densities[loc] = dists.remove(k);
                this.forParent[loc] = forParent.remove(k);

                final Taxa tk = clades.get(k);
                final int tkcount = tk.getTaxonCount();
                this.clades[loc] = new int[tkcount];
                for (int nt = 0; nt < tkcount; ++nt) {
                    final int taxonIndex = tree.getTaxonIndex(tk.getTaxon(nt));
                    this.clades[loc][nt] = taxonIndex;
                    if (taxonIndex < 0) {
                        throw new IllegalArgumentException("Taxon not found in tree: " + tk.getTaxon(nt));
                    }
                }
                taxaInOrder[loc] = tk;
                clades.remove(k);

                --loc;
            }
        }

        List<Integer>[] tio = new List[taxaInOrder.length];
        for (int k = 0; k < taxaInOrder.length; ++k) {
            tio[k] = new ArrayList<Integer>();
        }

        for (int k = 0; k < taxaInOrder.length; ++k) {
            for (int i = k + 1; i < taxaInOrder.length; ++i) {
                if (taxaInOrder[i].containsAll(taxaInOrder[k])) {
                    tio[i].add(k);
                    break;
                }
            }
        }

        this.taxaPartialOrder = new int[taxaInOrder.length][];
        for (int k = 0; k < taxaInOrder.length; ++k) {
            List<Integer> tiok = tio[k];

            this.taxaPartialOrder[k] = new int[tiok.size()];
            for (int j = 0; j < tiok.size(); ++j) {
                this.taxaPartialOrder[k][j] = tiok.get(j);
            }
        }
        this.freeHeights = new int[this.clades.length];

        for (int k = 0; k < this.clades.length; ++k) {
            int taken = 0;
            for (int i : this.taxaPartialOrder[k]) {
                taken += this.clades[i].length - (this.forParent[i] ? 0 : 1);
            }
            this.freeHeights[k] = this.clades[k].length - (this.forParent[k] ? 1 : 2) - taken;
            assert this.freeHeights[k] >= 0;
        }

        // true if clade is not contained in any other clade
        boolean[] maximal = new boolean[this.clades.length];
        for (int k = 0; k < this.clades.length; ++k) {
            maximal[k] = true;
        }

        for (int k = 0; k < this.clades.length; ++k) {
            for (int i : this.taxaPartialOrder[k]) {
                maximal[i] = false;
            }
        }

        rootCorrection = this.clades[this.clades.length - 1].length < tree.getExternalNodeCount();

        this.calibrationLogPDF = userPDF;

        this.correctionType = correctionType;

        if (userPDF == null) {
            if (!isYule) {
                throw new IllegalArgumentException("Sorry, not implemented: conditional calibration prior for this non Yule models.");
            }

            if (correctionType == CorrectionType.EXACT) {
                if (densities.length == 1) {
                    // closed form formula
                } else {
                    boolean anyParent = false;
                    for (boolean in : this.forParent) {
                        if (in) {
                            anyParent = true;
                        }
                    }
                    if (anyParent) {
                        throw new IllegalArgumentException("Sorry, not implemented: calibration on parent for more than one clade.");
                    }
                    if (densities.length == 2 && taxaInOrder[1].containsAll(taxaInOrder[0])) {
                        // closed form formulas
                    } else {
                        setUpTables(tree);
                        linsIter = new CalibrationLineagesIterator(this.clades, this.taxaPartialOrder, maximal,
                                tree.getExternalNodeCount());
                        lastHeights = new double[this.clades.length];
                    }
                }
            } else if (correctionType == CorrectionType.PEXACT) {
                setUpTables(tree);
            }
        }
    }

    private void setUpTables(Tree tree) {
        final int MAX_N = tree.getExternalNodeCount() + 1;
        double[] lints = new double[MAX_N];
        lc2 = new double[MAX_N];
        lfactorials = new double[MAX_N];
        lNR = new double[MAX_N];

        lints[0] = Double.NEGATIVE_INFINITY; //-infinity, should never be used
        lints[1] = 0.0;
        for (int i = 2; i < MAX_N; ++i) {
            lints[i] = Math.log(i);
        }

        lc2[0] = lc2[1] = Double.NEGATIVE_INFINITY;
        for (int i = 2; i < MAX_N; ++i) {
            lc2[i] = lints[i] + lints[i - 1] - lg2;
        }

        lfactorials[0] = 0.0;
        for (int i = 1; i < MAX_N; ++i) {
            lfactorials[i] = lfactorials[i - 1] + lints[i];
        }

        lNR[0] = Double.NEGATIVE_INFINITY; //-infinity, should never be used
        lNR[1] = 0.0;

        for (int i = 2; i < MAX_N; ++i) {
            lNR[i] = lNR[i - 1] + lc2[i];
        }
    }

    private boolean isMaximal(List<Taxa> taxa, int k) {
        final Taxa tk = taxa.get(k);
        for (int i = 0; i < taxa.size(); ++i) {
            if (i != k) {
                final Taxa ti = taxa.get(i);
                if (ti.containsAll(tk)) {
                    return false;
                }
            }
        }
        return true;
    }

    public double getCorrection(Tree tree, final double lam) {
        double logL = 0.0;
        final int nDists = densities.length;
        double hs[] = new double[nDists];

        for (int k = 0; k < nDists; ++k) {
            NodeRef c;
            final int[] taxk = clades[k];
            if (taxk.length > 1) {
                // check if monophyly and find node
                c = TreeUtils.getCommonAncestor(tree, taxk);

                if (TreeUtils.getLeafCount(tree, c) != taxk.length) {
                    return Double.NEGATIVE_INFINITY;
                }
            } else {
                c = tree.getNode(taxk[0]);
                assert forParent[k];
            }

            if (forParent[k]) {
                c = tree.getParent(c);
            }

            final double h = tree.getNodeHeight(c);
            logL += densities[k].logPdf(h);

            hs[k] = h;
        }

        if (Double.isInfinite(logL)) {
            return logL;
        }

        if (correctionType == CorrectionType.NONE) {
            return logL;
        }


        if (calibrationLogPDF == null) {
            switch (correctionType) {
                case EXACT: {
                    if (nDists == 1) {
                        logL -= logMarginalDensity(lam, tree.getExternalNodeCount(), hs[0], clades[0].length, forParent[0]);
                    } else if (nDists == 2 && taxaPartialOrder[1].length == 1) {
                        assert !forParent[0] && !forParent[1];
                        logL -= logMarginalDensity(lam, tree.getExternalNodeCount(), hs[0], clades[0].length,
                                hs[1], clades[1].length);
                    } else {

                        if (lastLam == lam) {
                            int k = 0;
                            for (; k < hs.length; ++k) {
                                if (hs[k] != lastHeights[k]) {
                                    break;
                                }
                            }
                            if (k == hs.length) {
                                return lastValue;
                            }
                        }

                        // the slow and painful way
                        double[] hss = new double[hs.length];
                        int[] ranks = new int[hs.length];
                        for (int k = 0; k < hs.length; ++k) {
                            int r = 0;
                            for (double h : hs) {
                                r += (h < hs[k]) ? 1 : 0;
                            }
                            ranks[k] = r + 1;
                            hss[r] = hs[k];
                        }
                        logL -= logMarginalDensity(lam, hss, ranks, linsIter);

                        lastLam = lam;
                        System.arraycopy(hs, 0, lastHeights, 0, lastHeights.length);
                        lastValue = logL;
                    }
                    break;
                }
                case APPROXIMATED: {

                    final double loglam = Math.log(lam);

                    int maxh = 0;
                    for (int k = 0; k < nDists; ++k) {
                        final double v = -lam * hs[k];

                        if (freeHeights[k] > 0) {
                            logL -= Math.log1p(-Math.exp(v)) * freeHeights[k];
                        }

                        logL -= v + loglam;

                        if (hs[k] > hs[maxh]) {
                            maxh = k;
                        }
                    }

                    if (rootCorrection || true) {
                        logL -= -(forParent[maxh] ? 0 : 1) * lam * hs[maxh];
                    }

                    if (Double.isNaN(logL)) {
                        logL = Double.NEGATIVE_INFINITY;
                    }
                    break;
                }
                case PEXACT: {
                    Arrays.sort(hs);
                    int cs[] = new int[nDists + 1];

                    final int internalNodeCount = tree.getInternalNodeCount();
                    for (int k = 0; k < internalNodeCount; ++k) {
                        final double nhk = tree.getNodeHeight(tree.getInternalNode(k));
                        int i = 0;
                        for (/**/; i < hs.length; ++i) {
                            if (hs[i] >= nhk) {
                                break;
                            }
                        }
                        if (i == hs.length) {
                            cs[i]++;
                        } else {
                            if (nhk < hs[i]) {
                                cs[i]++;
                            }
                        }
                    }

                    if (false) {
                        int l = nDists;
                        for (int i = 0; i < cs.length; ++i) {
                            l += cs[i];
                        }
                        assert l == internalNodeCount;
                    }
                    double ll = 0;

                    ll += cs[0] * Math.log1p(-Math.exp(-lam * hs[0])) - lam * hs[0] - lfactorials[cs[0]];
                    for (int i = 1; i < cs.length - 1; ++i) {
                        int c = cs[i];
                        ll += c * (Math.log1p(-Math.exp(-lam * (hs[i] - hs[i - 1]))) - lam * hs[i - 1]);
                        ll += -lam * hs[i] - lfactorials[c];
                    }
                    ll += -lam * (cs[nDists] + 1) * hs[nDists - 1] - lfactorials[cs[nDists] + 1];
                    ll += Math.log(lam) * nDists;

                    logL -= ll;
                    break;
                }
            }
        } else {
            final double value = calibrationLogPDF.getStatisticValue(0);
            if (Double.isNaN(value) || Double.isInfinite(value)) {
                logL = Double.NEGATIVE_INFINITY;
            } else {
                logL -= value;
            }
        }
        return logL;
    }

    private double logMarginalDensity(final double lam, int nTaxa, final double h, int nClade, boolean forParent) {
        double lgp;

        final double lh = lam * h;

        if (forParent) {
            // n(n+1) factor left out

            lgp = -2 * lh + Math.log(lam);
            if (nClade > 1) {
                lgp += (nClade - 1) * Math.log(1 - Math.exp(-lh));
            }
        } else {
            assert nClade > 1;

            lgp = -3 * lh + (nClade - 2) * Math.log(1 - Math.exp(-lh)) + Math.log(lam);

            // root is a special case
            if (nTaxa == nClade) {
                // n(n-1) factor left out
                lgp += lh;
            } else {
                // (n^3-n)/2 factor left out
            }
        }

        return lgp;
    }

    private double logMarginalDensity(final double lam, final int nTaxa, double h2, final int n, double h1, int nm) {

        assert h2 <= h1 && n < nm;

        final int m = nm - n;

        final double elh2 = Math.exp(-lam * h2);
        final double elh1 = Math.exp(-lam * h1);

        double lgl = 2 * Math.log(lam);

        lgl += (n - 2) * Math.log(1 - elh2);
        lgl += (m - 3) * Math.log(1 - elh1);

        lgl += Math.log(1 - 2 * m * elh1 + 2 * (m - 1) * elh2
                - m * (m - 1) * elh1 * elh2 + (m * (m + 1) / 2.) * elh1 * elh1
                + ((m - 1) * (m - 2) / 2.) * elh2 * elh2);

        if (nm < nTaxa) {
            /* lgl += Math.log(0.5*(n*(n*n-1))*(n+1+m)) */
            lgl -= lam * (h2 + 3 * h1);
        } else {
            /* lgl += Math.log(lam) /* + Math.log(n*(n*n-1)) */
            lgl -= lam * (h2 + 2 * h1);
        }

        return lgl;
    }

    private double logMarginalDensity(final double lam, double[] hs, int[] ranks, CalibrationLineagesIterator cli) {

        final int ni = cli.setup(ranks);

        final int nHeights = hs.length;

        double[] lehs = new double[nHeights + 1];
        lehs[0] = 0.0;
        for (int i = 1; i < lehs.length; ++i) {
            lehs[i] = -lam * hs[i - 1];
        }

        // assert maxRank == len(sit)
        boolean noRoot = ni == lehs.length;

        int nLevels = nHeights + (noRoot ? 1 : 0);

        double[] lebase = new double[nLevels];

        for (int i = 0; i < nHeights; ++i) {
            lebase[i] = lehs[i] + Math.log1p(-Math.exp(lehs[i + 1] - lehs[i]));
        }

        if (noRoot) {
            lebase[nHeights] = lehs[nHeights];
        }

        int[] linsAtLevel = new int[nLevels];

        int[][] joiners = cli.allJoiners();

        double val = 0;
        boolean first = true;

        int[][] linsInLevels;
        int ccc = 0;
        while ((linsInLevels = cli.next()) != null) {
            ccc++;
            double v = countRankedTrees(nLevels, linsInLevels, joiners, linsAtLevel);
            // 1 for root formula, 1 for kludge in iterator which sets root as 2 lineages
            if (noRoot) {
                final int ll = linsAtLevel[nLevels - 1] + 2;
                linsAtLevel[nLevels - 1] = ll;

                v -= lc2[ll] + lg2;
            }

            for (int i = 0; i < nLevels; ++i) {
                v += linsAtLevel[i] * lebase[i];
            }

            if (first) {
                val = v;
                first = false;
            } else {
                if (val > v) {
                    val += Math.log1p(Math.exp(v - val));
                } else {
                    val = v + Math.log1p(Math.exp(val - v));
                }
            }
        }

        double logc0 = 0.0;
        int totLin = 0;
        for (int i = 0; i < ni; ++i) {
            final int l = cli.nStart(i);
            if (l > 0) {
                logc0 += lNR[l];
                totLin += l;
            }
        }

        final double logc1 = lfactorials[totLin];

        double logc2 = nHeights * Math.log(lam);

        for (int i = 1; i < nHeights + 1; ++i) {
            logc2 += lehs[i];
        }

        if (!noRoot) {
            // we dont have an iterator for 0 free lineages
            logc2 += 1 * lehs[nHeights];
        }

        // Missing scale by total of all possible trees over all ranking orders.
        // Add it outside if needed for comparison.

        val += logc0 + logc1 + logc2;

        return val;
    }

    private double
    countRankedTrees(final int nLevels, final int[][] linsAtCrossings, final int[][] joiners, int[] linsAtLevel) {
        double logCount = 0;

        for (int i = 0; i < nLevels; ++i) {
            int sumLins = 0;
            for (int k = i; k < nLevels; ++k) {
                int[] lack = linsAtCrossings[k];
                int cki = lack[i];
                if (joiners[k][i] > 0) {
                    ++cki;
                    if (cki > 1) {
                        // can be 1 if iterator without lins - for joiners only - need to check this is correct
                        logCount += lc2[cki];
                    } //assert(cki >= 2);
                }
                final int l = cki - lack[i + 1];   //assert(l >= 0);
                logCount -= lfactorials[l];
                sumLins += l;
            }
            linsAtLevel[i] = sumLins;
        }

        return logCount;
    }

    // Flavour of marginal computation.
    private final CorrectionType correctionType;

    // Calibrated clades, each as a list of node ids.
    // Clades are partially ordered by inclusion - if X <= Y then X appears before Y.
    private final int[][] clades;

    // One calibration density associated with each clade
    private final Distribution[] densities;

    // true if density is for clade parent
    private final boolean[] forParent;

    // For each clade, lists the clades contained in it (using their index in clades)
    private final int[][] taxaPartialOrder;

    private final int[] freeHeights;

    private final boolean rootCorrection;

    // User provided function to calculate the marginal density
    private final Statistic calibrationLogPDF;

    // speedup constants
    private final double lg2 = Math.log(2.0);
    private double[] lc2;
    private double[] lNR;
    private double[] lfactorials;

    private CalibrationLineagesIterator linsIter = null;

    // simple cache of last result can go a long way in a big tree with a few calibration nodes, for non-global tree operators which do
    // not change the calibration nodes heights.

    double lastLam = Double.NEGATIVE_INFINITY;
    double[] lastHeights;
    double lastValue = Double.NEGATIVE_INFINITY;
}
