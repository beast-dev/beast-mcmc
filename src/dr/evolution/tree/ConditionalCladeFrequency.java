/**
 * 
 */
package dr.evolution.tree;

import java.io.IOException;
import java.io.Reader;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.HashMap;
import java.util.List;

import dr.evolution.io.Importer;
import dr.evomodel.tree.TreeTrace;

/**
 * @author shhn001
 * 
 * This class calculates the conditional clade probabilities for a set of trees.
 * It can be used to estimate the marginal posterior for a given tree. 
 *
 */
/**
 * @author Sebastian Hoehna
 *
 */
public class ConditionalCladeFrequency {
	
	private double EPSILON;
	
	private long samples = 0;
	
	private HashMap<BitSet, Double> cladeProbabilities;
	
	private HashMap<BitSet, HashMap<BitSet, Double>> cladeCoProbabilities;

	private HashMap<String, Integer> taxonMap;
	
	private TreeTrace[] traces;
	
	private int burnin;
	
	/**
	 * @param epsilon - the default number of occurences for each clade which wasn't observed to guarantee non-zero probabilities
	 */
	public ConditionalCladeFrequency(Tree tree, double epsilon) {

		// initializing global variables
		cladeProbabilities = new HashMap<BitSet, Double>();
		cladeCoProbabilities = new HashMap<BitSet, HashMap<BitSet,Double>>();
		
		// setting global variables
		EPSILON = epsilon;
		
		// extract the taxon
		taxonMap = getTaxonMap(tree);
		
	}
	
	/**
	 * @param traces - samples of trees in a tree traces array.
	 * @param epsilon - the default number of occurences for each clade which wasn't observed to guarantee non-zero probabilities
	 * @param burnIn - number of trees discarded from the trace
	 * @param verbose - hide the runtime status and outputs 
	 */
	public ConditionalCladeFrequency(TreeTrace[] traces, double epsilon, int burnIn, boolean verbose) {

		// initializing global variables
		cladeProbabilities = new HashMap<BitSet, Double>();
		cladeCoProbabilities = new HashMap<BitSet, HashMap<BitSet,Double>>();
		
		// setting global variables
		EPSILON = epsilon;
		this.traces = traces;
		
		// calculates the burn-in to 10% if it was set out of the boundaries
		int minMaxState = Integer.MAX_VALUE;
		for (TreeTrace trace : traces) {
			if (trace.getMaximumState() < minMaxState) {
				minMaxState = trace.getMaximumState();
			}
		}

		if (burnIn < 0 || burnIn >= minMaxState) {
			this.burnin = minMaxState / (10 * traces[0].getStepSize());
			if (verbose)
				System.out.println("WARNING: Burn-in larger than total number of states - using 10% of smallest trace");
		} else {
			this.burnin = burnIn;
		}

		// analyzing the whole trace -> reading the trees
		analyzeTrace(verbose);
	}
	
	
	/**
	 * Actually analyzes the trace given the burn-in.
	 * Each tree from the trace is read and the conditional clade frequencies incremented.
	 *
	 * @param verbose if true then progress is logged to stdout
	 */
	public void analyzeTrace(boolean verbose) {

		if (verbose) {
			if (traces.length > 1) System.out.println("Combining " + traces.length + " traces.");
		}

		// get first tree to extract the taxon
		Tree tree = getTree(0);
		taxonMap = getTaxonMap(tree);

		// read every tree from the trace
		for (TreeTrace trace : traces) {
			// do some output stuff
			int treeCount = trace.getTreeCount(burnin * trace.getStepSize());
			double stepSize = treeCount / 60.0;
			int counter = 1;

			if (verbose) {
				System.out.println("Analyzing " + treeCount + " trees...");
				System.out.println("0              25             50             75            100");
				System.out.println("|--------------|--------------|--------------|--------------|");
				System.out.print("*");
			}
			for (int i = 1; i < treeCount; i++) {
				// get the next tree
				tree = trace.getTree(i, burnin * trace.getStepSize());

				// add the tree and its clades to the frequencies
				addTree(tree);
				
				// some more output stuff
				if (i >= (int) Math.round(counter * stepSize) && counter <= 60) {
					if (verbose) {
						System.out.print("*");
						System.out.flush();
					}
					counter += 1;
				}
			}
			if (verbose) {
				System.out.println("*");
			}
		}
	}
	
	
	/**
	 * 
	 * Creates the report. The estimated posterior of the given tree is printed.
	 * 
	 * @param minCladeProbability clades with at least this posterior probability will be included in report.
	 * @throws IOException if general I/O error occurs
	 */
	public void report(Tree tree) throws IOException {

		System.err.println("making report");

		SimpleTree sTree = new SimpleTree(tree);
		System.out.println("Estimated marginal posterior by condiational clade frequencies:");
		System.out.println(getTreeProbability(sTree));
		
		System.out.flush();


	}
	
	/**
	 * 
	 * Calculates the probability of a given tree.
	 * 
	 * @param tree - the tree to be analyzed
	 * @return estimated posterior probability in log
	 */
	public double getTreeProbability(SimpleTree tree) {
		return calculateTreeProbabilityLog(tree);
	}

	/**
	 * 
	 * Calculates the probability of a given tree.
	 * 
	 * @param tree - the tree to be analyzed
	 * @return estimated posterior probability in log
	 */
	private double calculateTreeProbabilityLog(SimpleTree tree) {
		double prob = 0.0;

		List<Clade> clades = new ArrayList<Clade>();
		// get clades contained in the tree
		getNonComplementaryClades(taxonMap, tree, (SimpleNode) tree.getRoot(), clades, null, true);
		List<Clade> parent_clades = new ArrayList<Clade>();
		// get clades contained in the tree
		getClades(taxonMap, tree, (SimpleNode) tree.getRoot(), parent_clades, null);

		// for every clade multiply its conditional clade probability to the tree probability
		for (Clade c : clades) {
			// get the bits of the clade
			BitSet bs = c.getBits();
			BitSet parent = getParentClade(parent_clades, c);
			
			// set the occurrences to epsilon
			double tmp = EPSILON;
			double parentOccurrences = EPSILON;
			if (cladeProbabilities.containsKey(parent)) {
				// if we observed this clade in the trace, add the occurrences to epsilon
				parentOccurrences += cladeProbabilities.get(parent);
			}
			
			if (cladeCoProbabilities.containsKey(parent)) {
				// if we observed the parent clade
				HashMap<BitSet, Double> conditionalProbs = cladeCoProbabilities.get(parent);
				if (conditionalProbs.containsKey(bs)){
					// if we observed this conditional clade in the trace, add the occurrences to epsilon
					tmp += conditionalProbs.get(bs);					
				}
			}
			// multiply the conditional clade probability to the tree probability
			prob += Math.log(tmp / parentOccurrences);
		}
		
		return prob;
	}
	
	/**
	 * 
	 * get the i'th tree of the trace
	 * 
	 * @param index
	 * @return the ith tree of the trace
	 */
	public final Tree getTree(int index) {

		int oldTreeCount = 0;
		int newTreeCount = 0;
		for (TreeTrace trace : traces) {
			newTreeCount += trace.getTreeCount(burnin * trace.getStepSize());

			if (index < newTreeCount) {
				return trace.getTree(index - oldTreeCount, burnin * trace.getStepSize());
			}
			oldTreeCount = newTreeCount;
		}
		throw new RuntimeException("Couldn't find tree " + index);
	}

	
	/**
	 * 
	 * increments the number of occurrences for all conditional clades
	 * 
	 * @param tree - the tree to be added
	 */
	public void addTree(Tree tree) {
		
		samples++;

		List<Clade> clades = new ArrayList<Clade>();
		SimpleTree sTree = new SimpleTree(tree);
		// get all clades contained in the tree
		getClades(taxonMap, sTree, (SimpleNode) sTree.getRoot(), clades, null);

		// increment the occurrences of the clade and the conditional clade
		for (Clade c : clades) {
			// get the bits representing this clade
			BitSet bs = c.getBits();
			// find the parent clade
			BitSet parentClade = getParentClade(clades, c);
			
			Double frequency = 1.0;
			Double coFrequency = 1.0;
			HashMap<BitSet, Double> coFreqs;
			// increment the clade occurrences
			if (cladeProbabilities.containsKey(bs)) {
				frequency += cladeProbabilities.get(bs);
			}
			cladeProbabilities.put(bs, frequency);

			// increment the conditional clade occurrences
			if (cladeCoProbabilities.containsKey(parentClade)) {
				coFreqs = cladeCoProbabilities.get(parentClade);
			} else {
				// if it's the first time we observe the parent then we need a new list for its conditional clades
				coFreqs = new HashMap<BitSet, Double>();
				cladeCoProbabilities.put(parentClade, coFreqs);
			}

			// add the previous observed occurrences for this conditional clade
			if (coFreqs.containsKey(bs)) {
				coFrequency += coFreqs.get(bs);
			}
			coFreqs.put(bs, coFrequency);
		}
	}
	
	/**
	 * 
	 * Finds the parent of a given clade in a list of clades.
	 * The parent is the direct parent and not the grandparent or so.
	 * 
	 * @param clades - list of clades in which we are searching the parent
	 * @param child - the child of whom we are searching the parent
	 * @return the parent clade if found, otherwise itself
	 */
	private BitSet getParentClade(List<Clade> clades, Clade child) {
		BitSet parent = null;
		BitSet itself = child.getBits();

		// look in all clades of the list which contains the child and has the minimum cardinality (least taxa) -> that's the parent :-)
		for (int i = 0; i < clades.size(); i++) {
			BitSet tmp = clades.get(i).getBits();
			if (!itself.equals(tmp)
					&& containsClade(tmp, itself)) {
				if (parent == null || parent.cardinality() > tmp.cardinality()) {
					parent = tmp;
				}
			}
		}
		// if there isn't a parent, then you probably asked for the whole tree
		if (parent == null){
			parent = itself;
		}

		return parent;
	}
	
	/**
	 * 
	 * Checks if clade i contains clade j.
	 * 
	 * @param i - the parent clade
	 * @param j -  the child clade
	 * @return true, if i contains j
	 */
	private boolean containsClade(BitSet i, BitSet j) {
		BitSet tmpI = (BitSet) i.clone();
		BitSet tmpJ = (BitSet) j.clone();

		// just set the bits which are either in j but not in i or in i but not in j
		tmpI.xor(tmpJ);
		int numberOfBitsInEither = tmpI.cardinality();
		// which bits are just in i
		tmpI.and(i);
		int numberIfBitJustInContaining = tmpI.cardinality();

		// if the number of bits just in i is equal to the number of bits just in one of i or j
		// then i contains j
		return numberOfBitsInEither == numberIfBitJustInContaining;
	}

	/**
	 * 
	 * Creates a list with all clades of the tree
	 * 
	 * @param taxonMap - the lookup map for the taxon representing an index
	 * @param tree - the tree from which the clades are extracted
	 * @param node - the starting node. All clades below starting at this branch are added
	 * @param clades - the list in which the clades are stored
	 * @param bits - a bit set to which the current bits of the clades are added
	 */
	private void getClades(HashMap<String, Integer> taxonMap, SimpleTree tree,
			SimpleNode node, List<Clade> clades, BitSet bits) {

		// create a new bit set for this clade
		BitSet bits2 = new BitSet();

		// check if the node is external
		if (tree.isExternal(node)) {
			
			// if so, the only taxon in the clade is I
			int index = taxonMap.get(node.getTaxon().getId());
			bits2.set(index);

		} else {

			// otherwise, call all children and add its taxon together to one clade
			for (int i = 0; i < tree.getChildCount(node); i++) {
				SimpleNode child = (SimpleNode) tree.getChild(node, i);
				getClades(taxonMap, tree, child, clades, bits2);
			}
			// add my bit set to the list
			clades.add(new Clade(bits2, tree.getNodeHeight(node)));
		}

		// add my bit set to the bit set I was given
		// this is needed for adding all children clades together
		if (bits != null) {
			bits.or(bits2);
		}
	}
	
	/**
	 * 
	 * creates a list with all clades but just the non-complementary ones.
	 * ((A,B),(C,D)) just {A,B,C,D} and {A,B} are inserted.
	 * {A,B} is complementary to {C,D}
	 * 
	 * @param taxonMap - the lookup map for the taxon representing an index
	 * @param tree - the tree from which the clades are extracted
	 * @param node - the starting node. All clades below starting at this branch are added
	 * @param clades - the list in which the clades are stored
	 * @param bits - a bit set to which the current bits of the clades are added
	 * @param add - if this clade starting at the given node should be added
	 */
	private void getNonComplementaryClades(HashMap<String, Integer> taxonMap, SimpleTree tree,
			SimpleNode node, List<Clade> clades, BitSet bits, boolean add) {

		// create a new bit set for this clade
		BitSet bits2 = new BitSet();

		// check if the node is external
		if (tree.isExternal(node)) {
			
			// if so, the only taxon in the clade is I
			int index = taxonMap.get(node.getTaxon().getId());
			bits2.set(index);

		} else {

			// otherwise, call all children and add its taxon together to one clade
			for (int i = 0; i < tree.getChildCount(node); i++) {
				SimpleNode child = (SimpleNode) tree.getChild(node, i);
				// add just my first child to the list
				// the second child is complementary to the first
				getNonComplementaryClades(taxonMap, tree, child, clades, bits2, i==0);
			}
			if (add){
				clades.add(new Clade(bits2, tree.getNodeHeight(node)));
			}
		}

		// add my bit set to the bit set I was given
		// this is needed for adding all children clades together
		if (bits != null) {
			bits.or(bits2);
		}
	}
	
	/**
	 * 
	 * @return a mapping between the taxon and indices
	 */
	private HashMap<String, Integer> getTaxonMap(Tree tree) {
		HashMap<String, Integer> tm = new HashMap<String, Integer>();

		for (int i = 0; i < tree.getTaxonCount(); i++) {
			tm.put(tree.getTaxon(i).getId(), new Integer(i));
		}

		return tm;
	}
	
	/**
	 * @param reader  the readers to be analyzed
	 * @param burnin  the burnin in states
	 * @param verbose true if progress should be logged to stdout
	 * @return an analyses of the trees in a log file.
	 * @throws java.io.IOException if general I/O error occurs
	 */
	public static ConditionalCladeFrequency analyzeLogFile(Reader[] reader, double e, int burnin, boolean verbose) throws IOException {

		TreeTrace[] trace = new TreeTrace[reader.length];
		for (int i = 0; i < reader.length; i++) {
			try {
				trace[i] = TreeTrace.loadTreeTrace(reader[i]);
			} catch (Importer.ImportException ie) {
				throw new RuntimeException(ie.toString());
			}
			reader[i].close();

		}

		return new ConditionalCladeFrequency(trace, e, burnin, verbose);
	}
}
