package dr.app.bss;

import dr.evolution.tree.Tree;
import dr.evolution.util.Taxa;

public class TreesTableRecord {

	String name = "";
	Taxa taxa = new Taxa();
	Tree tree = null;

	public TreesTableRecord() {
	}

	public TreesTableRecord(String name, Tree tree, Taxa taxa) {
		this.name = name;
		this.tree = tree;
		this.taxa = taxa;
	}

	public TreesTableRecord(String name, Taxa taxa) {
		this.name = name;
		this.taxa = taxa;
		this.tree = null;
	}

	public int getTaxaCount() {
		return taxa.getTaxonCount();
	}

	public String getName() {
		return name;
	}

	public void setName(String name) {
		this.name = name;
	}

	public Tree getTree() {
		return tree;
	}

	public void setTree(Tree tree) {
		this.tree = tree;
	}

	public Taxa getTaxa() {
		return taxa;
	}

	public void setTaxa(Taxa taxa) {
		this.taxa = taxa;
	}

	public String toString() {
		return getName();
	}

}
