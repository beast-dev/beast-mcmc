package dr.app.bss;

import java.io.Serializable;

import dr.evolution.tree.Tree;
import dr.evolution.util.Taxa;

@SuppressWarnings("serial")
public class TreesTableRecord implements Serializable {

	private String name = "";
	private Taxa taxa = null;
	private Tree tree = null;

	private boolean treeSet = false;
	private boolean taxaSet = false;

	public TreesTableRecord() {
	}// END: Constructor

	public TreesTableRecord(String name, Tree tree) {
		this.name = name;
		this.tree = tree;

		treeSet = true;
		applyTreeName();
	}// END: Constructor

	public TreesTableRecord(String name, Taxa taxa) {
		this.name = name;
		this.taxa = taxa;

		taxaSet = true;
		applyTaxaName();
	}// END: Constructor

	private void applyTaxaName() {

		if (taxaSet && name != null) {

			for (int i = 0; i < taxa.getTaxonCount(); i++) {
				taxa.setTaxonAttribute(i, Utils.TREE_FILENAME, name);
			}
		}
	}// END: applyTaxaName

	private void applyTreeName() {

		if (treeSet && name != null) {

			for (int i = 0; i < tree.getTaxonCount(); i++) {
				tree.getTaxon(i).setAttribute(Utils.TREE_FILENAME, name);
			}
		}
	}// END: applyTreeName

	public int getTaxaCount() {

		int taxaCount = 0;

		if (taxaSet) {

			taxaCount = taxa.getTaxonCount();

		} else if (treeSet) {

			taxaCount = tree.getTaxonCount();

		} else {
			// do nothing
		}

		return taxaCount;
	}// END: getTaxaCount

	public String getName() {
		return name;
	}

	public Tree getTree() {
		return tree;
	}

	public Taxa getTaxa() {
		return taxa;
	}

	public String toString() {
		return getName();
	}

	public boolean isTreeSet() {
		return this.treeSet;
	}
	
	public boolean isTaxaSet() {
		return this.taxaSet;
	}
	
}// END: class
