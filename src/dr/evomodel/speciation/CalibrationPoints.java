package dr.evomodel.speciation;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxa;
import dr.inference.model.Statistic;
import dr.math.distributions.Distribution;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Joseph Heled
 *         Date: 8/06/2011
 */
public class CalibrationPoints {

    public CalibrationPoints(Tree tree, boolean isYule, List<Distribution> dists, List<Taxa> taxa,
                           List<Boolean> forParent, Statistic userPDF, boolean approxOK) {

        this.distributions = new Distribution[dists.size()];
        this.taxa = new int[taxa.size()][];
        this.forParent = new boolean[taxa.size()];

        for(int k = 0; k < taxa.size(); ++k) {
            final Taxa tk = taxa.get(k);
            for(int i = k+1; i < taxa.size(); ++i) {
                final Taxa ti = taxa.get(i);
                if( ti.containsAny(tk) ) {
                    if( ! (ti.containsAll(tk) || tk.containsAll(ti) ) ) {
                        throw new IllegalArgumentException("Overlapping clades??");
                    }
                }
            }
        }

        Taxa[] taxaInOrder = new Taxa[taxa.size()];

        {
            int loc = taxa.size() - 1;
            while( loc >= 0 ) {
                //  place maximal clades at end one at a time
                int k = 0;
                for(/**/; k < taxa.size(); ++k) {
                    if( isMaximal(taxa, k) ) {
                        break;
                    }
                }
                this.distributions[loc] = dists.remove(k);
                this.forParent[loc] = forParent.remove(k);


                final Taxa tk = taxa.get(k);
                final int tkcount = tk.getTaxonCount();
                this.taxa[loc] = new int[tkcount];
                for(int nt = 0; nt < tkcount; ++nt) {
                    this.taxa[loc][nt] = tree.getTaxonIndex(tk.getTaxon(nt));
                }
                taxaInOrder[loc] = tk;
                taxa.remove(k);

                --loc;
            }
        }

        List<Integer>[] tio = new List[taxaInOrder.length];
        for(int k = 0; k < taxaInOrder.length; ++k) {
            tio[k] = new ArrayList<Integer>();
        }

        for(int k = 0; k < taxaInOrder.length; ++k) {
            for(int i = k+1; i < taxaInOrder.length; ++i) {
                if( taxaInOrder[i].containsAll(taxaInOrder[k]) ) {
                    tio[i].add(k);
                    break;
                }
            }
        }
        this.taxaPartialOrder = new int[taxaInOrder.length][];
        for(int k = 0; k < taxaInOrder.length; ++k) {
            List<Integer> tiok = tio[k];

            this.taxaPartialOrder[k] = new int[tiok.size()];
            for(int j = 0; j < tiok.size(); ++j) {
                this.taxaPartialOrder[k][j] = tiok.get(j);
            }
        }
        this.freeHeights = new int[this.taxa.length];

        for(int k = 0; k < this.taxa.length; ++k) {
            int taken = 0;
            for( int i : this.taxaPartialOrder[k] ) {
                taken += this.taxaPartialOrder[i].length - (this.forParent[i] ? 0 : 1);
            }
            this.freeHeights[k] = this.taxa[k].length - (this.forParent[k] ? 1 : 2) - taken;
            assert this.freeHeights[k] >= 0;
        }

        this.maximal = new boolean[this.taxa.length];
        for(int k = 0; k < this.taxa.length; ++k) {
            maximal[k] = true;
        }
        for(int k = 0; k < this.taxa.length; ++k) {
           for( int i : this.taxaPartialOrder[k] ) {
               maximal[i] = false;
           }
        }
        
        this.calibrationLogPDF = userPDF;

        if( userPDF == null ) {
            if( ! isYule ) {
                throw new IllegalArgumentException("Sorry, not implemented: conditional calibration prior for this non Yule models.");
            }

            if( ! approxOK ) {
                if( distributions.length > 2 ) {
                    throw new IllegalArgumentException("Sorry, not implemented: multiple internal calibrations - please provide the " +
                            "log marginal explicitly.");
                }

                if( distributions.length == 2 ) {

                    if( ! taxaInOrder[1].containsAll(taxaInOrder[0]) ) {
                        throw new IllegalArgumentException( "Sorry, not implemented: two non-nested clades.");
                    }

                    if( this.forParent[0] || this.forParent[1] ) {
                        throw new IllegalArgumentException("Sorry, not implemented: calibration on parent for more than one clade.");
                    }
                }
            }
        }
    }

    private boolean isMaximal(List<Taxa> taxa, int k) {
        final Taxa tk = taxa.get(k);
        for(int i = 0; i < taxa.size(); ++i) {
            if( i != k ) {
                final Taxa ti = taxa.get(i);
                if( ti.containsAll(tk) ) {
                    return false;
                }
            }
        }
        return true;
    }

    public double getCorrection(Tree tree, final double lam) {
        double logL = 0.0;
        final int nDists = distributions.length;
        double hs[] = new double[nDists];

        for(int k = 0; k < nDists; ++k) {
            NodeRef c;
            final int[] taxk = taxa[k];
            if( taxa.length > 1 ) {
                // check if monophyly and find node
                c = Tree.Utils.getCommonAncestor(tree, taxk);

                if( Tree.Utils.getLeafCount(tree, c) != taxk.length ) {
                    return Double.NEGATIVE_INFINITY;
                }
            } else {
                c = tree.getNode(taxk[0]);
                assert forParent[k];
            }

            if( forParent[k] ) {
                c = tree.getParent(c);
            }

            final double h = tree.getNodeHeight(c);
            logL += distributions[k].logPdf(h);

            hs[k] = h;
        }

        if( calibrationLogPDF == null ) {
            if( nDists == 1 ) {
                logL -= logMarginalDensity(lam, tree.getExternalNodeCount(), hs[0], taxa[0].length, forParent[0]);
            } else if( nDists == 2 && taxaPartialOrder[1].length == 1
                    && !forParent[0] && !forParent[1] ) {
                logL -= logMarginalDensity(lam, tree.getExternalNodeCount(), hs[0], taxa[0].length,
                                           hs[1], taxa[1].length);
            } else {
                final double loglam = Math.log(lam);
                if( false ) {
                    for(int k = 0; k < nDists; ++k) {
                        final double v = -lam * hs[k];
                        logL -= v  + loglam;
                    }
                } else {
                    int maxh = 0;
                    for(int k = 0; k < nDists; ++k) {
                        final double v = -lam * hs[k];

                        logL -= Math.log1p(-Math.exp(v)) * freeHeights[k];

                        logL -= v * ((forParent[k] ? 1 : 2) - (maximal[k] ? 0 : 1)) + loglam;
                        if( hs[k] > hs[maxh] ) {
                            maxh = k;
                        }
                    }
                    logL -= -(forParent[maxh] ? 0 : 1) * lam * hs[maxh];
                }
            }
        } else {
            final double value = calibrationLogPDF.getStatisticValue(0);
            if( Double.isNaN(value) || Double.isInfinite(value) )  {
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

        if( forParent ) {
            // n(n+1) factor left out

            lgp = -2 * lh + Math.log(lam);
            if( nClade > 1 ) {
                lgp += (nClade-1) * Math.log(1 - Math.exp(-lh));
            }
        } else {
            assert nClade > 1;

            lgp = -3 * lh + (nClade-2) * Math.log(1 - Math.exp(-lh)) + Math.log(lam);

            // root is a special case
            if( nTaxa == nClade ) {
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

        final double elh2 = Math.exp(-lam*h2);
        final double elh1 = Math.exp(-lam*h1);

        double lgl= 2 * Math.log(lam);

        lgl += (n-2) * Math.log(1-elh2);
        lgl += (m-3) * Math.log(1-elh1);

        lgl += Math.log(1 - 2*m*elh1 + 2*(m-1)*elh2
                - m*(m-1)*elh1*elh2 + (m*(m+1)/2.)*elh1*elh1
                + ((m-1)*(m-2)/2.)*elh2*elh2);

        if( nm < nTaxa ) {
            /* lgl += Math.log(0.5*(n*(n*n-1))*(n+1+m)) */
            lgl -= lam*(h2+3*h1);
        } else {
            /* lgl += Math.log(lam) /* + Math.log(n*(n*n-1)) */
            lgl -= lam*(h2+2*h1);
        }

        return lgl;
    }

    private  final Distribution[] distributions;
    private  final int[][] taxa;
    private  final boolean[] forParent;
    private  final int[][] taxaPartialOrder;
    private  final int[] freeHeights;

    private  final boolean[] maximal;

    private  final Statistic calibrationLogPDF;
}
