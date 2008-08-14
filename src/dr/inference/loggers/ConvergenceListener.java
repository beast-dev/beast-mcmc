/**
 * 
 */
package dr.inference.loggers;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxon;
import dr.inference.markovchain.MarkovChain;
import dr.inference.markovchain.MarkovChainListener;
import dr.inference.model.Model;

/**
 * @author shhn001
 *
 */
public class ConvergenceListener implements MarkovChainListener {
	
	private MarkovChain chain;
	
	private final int LOG_EVERY = 100;
	
	private final int CHECK_EVERY = 10000;
	
	private final double THRESHOLD = 0.001;

	private Tree tree = null;

	private BufferedWriter bw = null;

	private long trees = 0;

	private HashMap<BitSet, Double> referenceCladeFrequencies = null;

	private HashMap<BitSet, Double> cladeOccurences = null;

	private String outputFilename;

	private HashMap<String, Integer> taxonMap;
	
	private String [] taxonArray;



	/**
	 * 
	 */
	public ConvergenceListener(MarkovChain chain, String outputFilename, Tree tree) {
		this.chain = chain;		
		this.outputFilename = outputFilename;
		this.tree = tree;
		
		taxonMap = new HashMap<String, Integer>();
		cladeOccurences = new HashMap<BitSet, Double>();
		referenceCladeFrequencies = new HashMap<BitSet, Double>();
		
		startListening();
	}
	
	private void startListening() {

		File f = new File(outputFilename);
		if (f.exists()) {
			f.delete();
		}
		
		try {
			bw = new BufferedWriter(new FileWriter(f));
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	private void fillTaxonMap(Tree tree){
		addTaxon(tree, tree.getRoot());
	}
	
	private void addTaxon(Tree tree, NodeRef node){
		if (tree.isExternal(node)){
			Taxon taxon = tree.getNodeTaxon(node);
			String name = taxon.getId();
			taxonMap.put(name, taxonMap.size());
		}
		else {
			int children = tree.getChildCount(node);
			for (int i=0; i<children; i++){
				addTaxon(tree, tree.getChild(node, i));
			}
		}
	}

	private void addTree(Tree tree) {
		if (trees == 0){
			fillTaxonMap(tree);
			taxonArray = new String[taxonMap.size()];
			Set<String> taxon = taxonMap.keySet();
			for (String s : taxon){
				taxonArray[taxonMap.get(s)] = s;
			}
		}
		
		trees++;
		List<BitSet> clades = getClades(tree);

		for (BitSet c : clades) {
			if (!cladeOccurences.containsKey(c)) {
				cladeOccurences.put(c, 1.0);
			} else {
				cladeOccurences.put(c, cladeOccurences.get(c) + 1);
			}
		}
	}

	private List<BitSet> getClades(Tree tree) {
		List<BitSet> clades = new ArrayList<BitSet>();

		NodeRef root = tree.getRoot();
		fillClades(clades, root, tree);

		return clades;
	}

	private BitSet fillClades(List<BitSet> clades, NodeRef root, Tree tree) {
		BitSet clade = new BitSet();

		if (!tree.isExternal(root)) {
			clade.or(fillClades(clades, tree.getChild(root, 0), tree));
			clade.or(fillClades(clades, tree.getChild(root, 1), tree));
			clades.add(clade);
		} else {
			Taxon taxon = tree.getNodeTaxon(root);
			String name = taxon.getId();
			clade.set(taxonMap.get(name));
		}

		return clade;

	}

	private double getMaxCladeDistance() {

		return getMaxDeviation();
	}

	private double getMaxDeviation() {
		double[] deviation = getDeviations();

		double max = 0;
		for (double m : deviation) {
			if (m > max) {
				max = m;
			}
		}
		return max;
	}
	
	private void setReferenceClades(){
		referenceCladeFrequencies.clear();
		Set<BitSet> keys = null;
		keys = cladeOccurences.keySet();
		
		for (BitSet key : keys) {
			double value = (cladeOccurences.get(key)* 100.0) / (double) trees;
			referenceCladeFrequencies.put(key, value);
		}
	}

	private double[] getDeviations() {
		double[] deviations = new double[cladeOccurences.size()];

		Set<BitSet> keys = null;
		keys = cladeOccurences.keySet();

		Iterator<BitSet> it = keys.iterator();
		for (int i = 0; i < deviations.length; i++) {
			BitSet k = it.next();
			double r = 0.0;
			if (referenceCladeFrequencies.containsKey(k)) {
				r = referenceCladeFrequencies.get(k);
			}
			double o = 0.0;
			if (cladeOccurences.containsKey(k)) {
				o = cladeOccurences.get(k);
			}
			deviations[i] = Math.abs(r - ((o * 100.0) / (double) trees));
		}

		return deviations;
	}

	/* (non-Javadoc)
	 * @see dr.inference.markovchain.MarkovChainListener#bestState(int, dr.inference.model.Model)
	 */
	public void bestState(int state, Model bestModel) {
		// do nothing
	}

	/* (non-Javadoc)
	 * @see dr.inference.markovchain.MarkovChainListener#currentState(int, dr.inference.model.Model)
	 */
	public void currentState(int state, Model currentModel) {
		if (state % LOG_EVERY == 0) {
			addTree(tree);
			if (state % (LOG_EVERY * CHECK_EVERY) == 0) {
				double distance = getMaxCladeDistance();
				System.out.println("Distance to the last sample: " + distance);
				if (distance < THRESHOLD){
					chain.pleaseStop();
				}
				setReferenceClades();
			}
		}
	}

	/* (non-Javadoc)
	 * @see dr.inference.markovchain.MarkovChainListener#finished(int)
	 */
	public void finished(int chainLength) {
		printCladeFrequencies(chainLength);
		try {
			bw.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	private void printCladeFrequencies(int chainLength) {

		try {
			bw.write("Finished after " + chainLength + " iterations.");
			bw.newLine();
			bw.newLine();
			
			Set<BitSet> keys = cladeOccurences.keySet();
			
			for (BitSet key : keys) {
				double value = cladeOccurences.get(key);
				double percent = ((int)((value * 10000.0) / trees))/100.0;
				String summary = value + "\t" + percent + "\t" + "(";
				int oldIndex = 0;
				for (int i=0; i<key.cardinality(); i++){
					int index = key.nextSetBit(oldIndex);
					oldIndex = index+1;
					summary += taxonArray[index];
					
					if (i+1 < key.cardinality()){
						summary += ",";
					}
				}
				summary += ")";
				bw.write(summary);
				bw.newLine();
				bw.flush();
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

}
