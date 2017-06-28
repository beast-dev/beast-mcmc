package dr.evolution.tree.treemetrics;

import dr.evolution.tree.Clade;
import dr.evolution.tree.CladeSet;
import dr.evolution.tree.Tree;

import java.util.*;

import static dr.evolution.tree.treemetrics.TreeMetric.Utils.checkTreeTaxa;

/**
 * @author Andrew Rambaut
 * @version $Id$
 */
public class RobinsonFouldsMetric implements TreeMetric {

	public RobinsonFouldsMetric() {

	}

	@Override
	public double getMetric(Tree tree1, Tree tree2) {

		checkTreeTaxa(tree1, tree2);

		Set<Clade> clades1 = Clade.getCladeSet(tree1);
		Set<Clade> clades2 = Clade.getCladeSet(tree2);

		clades1.removeAll(clades2);

		return clades1.size();
	}

}
