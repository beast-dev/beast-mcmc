package dr.evolution.tree.treemetrics;

import dr.evolution.tree.Clade;
import dr.evolution.tree.Tree;

import java.util.*;

import static dr.evolution.tree.treemetrics.TreeMetric.Utils.checkTreeTaxa;

/**
 * @author Andrew Rambaut
 * @version $Id$
 */
public class RobinsonFouldsMetric implements TreeMetric {
	public static Type TYPE = Type.ROBINSON_FOULDS;

	public RobinsonFouldsMetric() {

	}

	@Override
	public double getMetric(Tree tree1, Tree tree2) {

		checkTreeTaxa(tree1, tree2);

		Set<Clade> clades1 = Clade.getCladeSet(tree1);
		Set<Clade> clades2 = Clade.getCladeSet(tree2);

		clades1.removeAll(clades2);

		// Technically RF would be twice this because it doesn't assume
		// the same set of tips in both trees (so may have a different
		// number of clades missing from each).
		return clades1.size();
	}

	@Override
	public Type getType() {
		return TYPE;
	}

	@Override
	public String toString() {
		return getType().getShortName();
	}

	// todo - add in Citable:
	// Robinson, D. R.; Foulds, L. R. (1981). "Comparison of phylogenetic trees". Mathematical Biosciences. 53: 131-147. doi:10.1016/0025-5564(81)90043-2.
}
