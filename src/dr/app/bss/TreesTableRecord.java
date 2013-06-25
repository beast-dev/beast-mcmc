package dr.app.bss;

import dr.evolution.tree.Tree;
import dr.evolution.util.Taxa;

public class TreesTableRecord {

	private String name = "";
	private Taxa taxa = null;// new Taxa();
	private Tree tree = null;

	public boolean treeSet = false;
	public boolean taxaSet = false;

	public TreesTableRecord() {
	}// END: Constructor

	public TreesTableRecord(String name, Tree tree
	// , Taxa taxa
	) {
		this.name = name;
		this.tree = tree;
		// this.taxa = null;

		treeSet = true;
		applyTreeName();
	}// END: Constructor

	public TreesTableRecord(String name, Taxa taxa) {
		this.name = name;
		this.taxa = taxa;
		// this.tree = null;

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

	// public void setName(String name) {
	// this.name = name;
	// applyTaxaName();
	// }

	public Tree getTree() {
		return tree;
	}

	// public void setTree(Tree tree) {
	// this.tree = tree;
	// }

	public Taxa getTaxa() {
		return taxa;
	}

	// public void setTaxa(Taxa taxa) {
	// this.taxa = taxa;
	// }

	public String toString() {
		return getName();
	}

}// END: class
