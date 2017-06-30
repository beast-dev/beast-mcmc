package dr.evolution.tree.treemetrics;

import dr.evolution.tree.Clade;
import dr.evolution.tree.CladeSet;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;

import java.util.*;

import static dr.evolution.tree.treemetrics.TreeMetric.Utils.checkTreeTaxa;

/**
 * For each clade_j in treeB, find the MRCA_j of the taxa in clade_j in treeA.
 *
 * obviously, if clade_i exists in treeB then clade_i == MRCA_i.
 *
 * Then find the sum of squares of the differences in height:
 *
 * d = sqrt( sum across i[ (height(clade_i) - height(MRCA_i))^2 ] +
 * 		sum across j[ (height(clade_j) - height(MRCA_j))^2 ] )
 *
 * The rationale is that if a clade moves then this includes the size of movement
 * across the MRCA node and down again. I.e., a clade that moves from one side
 * of the tree to the other scores the difference between the height of that node
 * and the root and back down to the height of the node in the other tree.
 *
 *   +---------A               +---------A
 * +-+                       +-+
 * | +---------B             | |+--------B
 * +                 ====>   + ++
 * | +---------C             |  +--------C
 * +-+                       |
 *   +---------D             +-----------D
 *
 * 	height(tree1)	tmrca(tree2)	diff
 * AB 	10		10		0
 * CD	10		12		2
 * ABCD	12		12		0
 *
 * 	height(tree1)	tmrca(tree2)	diff
 * BC 	9		12		3
 * ABC	10		10		0
 * ABCD	12		12		0
 *
 * So the score is sqrt(2^2 + 3^2) = sqrt(13)
 *
 * Scores much less than this:
 *
 *   +----------A               +----------A
 * + +                        +-+
 * | +----------B             | |        +-B
 * +                 ====>    + +--------+
 * |        +---C             |          +-C
 * +--------+                 |
 *          +---D             +------------D
 *
 * 	height(tree1)	tmrca(tree2)	diff
 * AB 	10		10		0
 * CD	4		12		8
 * ABCD	12		12		0
 *
 * 	height(tree1)	tmrca(tree2)	diff
 * BC 	2		12		10
 * ABC	10		10		0
 * ABCD	12		12		0
 *
 * So the score is sqrt(8^2 + 10^2) = sqrt(164)
 *
 * @author Andrew Rambaut
 * @version $Id$
 */
public class CladeHeightMetric implements TreeMetric {
	public static Type TYPE = Type.CLADE_HEIGHT;

	public CladeHeightMetric() {

	}

	@Override
	public double getMetric(Tree tree1, Tree tree2) {

		checkTreeTaxa(tree1, tree2);

		Set<Clade> clades1 = Clade.getCladeSet(tree1);
		Set<Clade> clades2 = Clade.getCladeSet(tree2);

		return getDistance(clades1, clades2);
	}

	private double getDistance(Set<Clade> clades1, Set<Clade> clades2) {

	    double distance = 0.0;

	    for (Clade clade1 : clades1) {
	        double height1 = clade1.getHeight();

	        Clade clade2 = findMRCA(clade1, clades2);
	        double height2 = clade2.getHeight();

	        distance += (height1 - height2) * (height1 - height2);
	    }

	    for (Clade clade2 : clades2) {
	        double height2 = clade2.getHeight();

	        Clade clade1 = findMRCA(clade2, clades1);
	        double height1 = clade1.getHeight();

	        distance += (height1 - height2) * (height1 - height2);
	    }

	    return Math.sqrt(distance);
	}

	private Clade findMRCA(Clade clade1, Set<Clade> clades) {

	    for (Clade clade2 : clades) {
	        if (isMRCA(clade1, clade2)) {
	            return clade2;
	        }
	    }

	    return null;
	}

	private boolean isMRCA(Clade clade1, Clade clade2) {
	    if (clade1.getSize() > clade2.getSize()) {
	        return false;
	    }

	    tmpBits.clear();
	    tmpBits.or(clade1.getBits());
	    tmpBits.and(clade2.getBits());

	    return tmpBits.cardinality() == clade1.getSize();
	}

	@Override
	public Type getType() {
		return TYPE;
	}

	@Override
	public String toString() {
		return getType().getShortName();
	}

	BitSet tmpBits = new BitSet();

}
