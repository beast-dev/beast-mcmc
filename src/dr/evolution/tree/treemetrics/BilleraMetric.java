package dr.evolution.tree.treemetrics;


import dr.evolution.tree.Clade;
import dr.evolution.tree.CladeSet;
import dr.evolution.tree.Tree;
import jebl.evolution.trees.TreeBiPartitionInfo;

import java.util.ArrayList;
import java.util.List;
import java.util.Set;

import static dr.evolution.tree.treemetrics.TreeMetric.Utils.checkTreeTaxa;

/**
 * Billera tree distance - sum of change in branch lengths required to transform one tree to the second
 *
 * Note that this interface is not optimal for a large set where all pairs are required.
 * Creating TreeBiPartitionInfo's as a pre step is better unless memory is an issue.
 * 
 * @author Joseph Heled
 * @version $Id$
 */
public class BilleraMetric implements TreeMetric {

    @Override
    public double getMetric(Tree tree1, Tree tree2) {
        checkTreeTaxa(tree1, tree2);

        Set<Clade> clades1 = Clade.getCladeSet(tree1);
        Set<Clade> clades2 = Clade.getCladeSet(tree2);
        return distance(clades1, clades2, true);
    }

    public double distance(Set<Clade> clades1, Set<Clade> clades2, boolean norm) {

//        for( BiPartiotionInfo k : t2.all.values() ) {
//            k.has = false;
//        }
//        double din = 0;
//        double dout = 0;
//
//        for( Map.Entry<FixedBitSet, BiPartiotionInfo> k : t1.all.entrySet() ) {
//            final BiPartiotionInfo info = t2.all.get(k.getKey());
//            final double b1 = t1.t.getLength(k.getValue().n);
//            double dif;
//            if( info != null ) {
//
//                final double b2 = t2.t.getLength(info.n);
//                info.has = true;
//
//                dif = Math.abs(b1 - b2);
//            } else {
//                dif = b1;
//            }
//            if( norm == DistanceNorm.NORM1 ) {
//                //d += dif;
//                din += dif;
//            } else {
//                //d += dif * dif;
//                din += dif * dif;
//            }
//        }
//
//        for( BiPartiotionInfo info : t2.all.values() ) {
//            if( !info.has ) {
//                final double dif = t2.t.getLength(info.n);
//                if( norm == DistanceNorm.NORM1 ) {
//                    //d += dif;
//                    dout += dif;
//                } else {
//                    //d += dif * dif;
//                    dout += dif * dif;
//                }
//            }
//        }
//        double d = din + dout;
//        return ( norm == DistanceNorm.NORM1 ) ? d : Math.sqrt(d);
        return 0;
    }
}
