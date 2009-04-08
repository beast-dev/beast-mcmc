/**
 * 
 */
package dr.evolution.tree;

/**
 * @author Sebastian Hoehna
 *
 */
public interface ImportanceDistribution {
	
	/**
	 * 
	 * @param tree
	 *            - the tree to be added
	 */
	public void addTree(Tree tree);	
	
	/**
	 * 
	 * Splits a clade into two sub-clades according to the importance distribution
	 * 
	 * @param parent - the clade which is split
	 * @param children - a call by reference parameter which is an empty, two element array of clades at time of call and contains the to sub clades afterwards
	 * @return the chance for this split
	 */
	public double splitClade(Clade parent, Clade[] children);
	
	/**
	 * 
	 * Calculates the probability of a given tree.
	 * 
	 * @param tree
	 *            - the tree to be analyzed
	 * @return estimated posterior probability in log
	 */
	public double getTreeProbability(Tree tree);

}
