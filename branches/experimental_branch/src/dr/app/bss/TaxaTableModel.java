package dr.app.bss;

import javax.swing.table.AbstractTableModel;

/**
 * @author Filip Bielejec
 * @version $Id$
 */
@SuppressWarnings("serial")
public class TaxaTableModel extends AbstractTableModel {

	private PartitionDataList dataList;

	private String[] COLUMN_NAMES = { "Name", "Height", "Data-set" };
	private double[] heights = null;
	private String[] trees = null;

	public final static int NAME_INDEX = 0;
	public final static int HEIGHT_INDEX = 1;
	public final static int TAXA_SET_INDEX = 2;

	public TaxaTableModel(PartitionDataList dataList) {
		this.dataList = dataList;
	}// END: Constructor

	public int getColumnCount() {
		return COLUMN_NAMES.length;
	}

	public int getRowCount() {
		return this.dataList.allTaxa.getTaxonCount();
	}

	public Class<? extends Object> getColumnClass(int c) {
		return getValueAt(0, c).getClass();
	}

	public boolean isCellEditable(int row, int col) {
		return false;
	}

	public String getColumnName(int column) {
		return COLUMN_NAMES[column];
	}

	public Object getValueAt(int row, int col) {
		switch (col) {

		case NAME_INDEX:
			return this.dataList.allTaxa.getTaxonId(row);

		case HEIGHT_INDEX:

			if (heights != null) {
				return heights[row];
			} else {
				return 0.0;
			}

		case TAXA_SET_INDEX:

			if (trees != null) {
				return trees[row];
			} else {
				return "";
			}

		default:
			return null;

		}// END: switch
	}// END: getValueAt

	private void getHeights() {

		heights = new double[dataList.allTaxa.getTaxonCount()];
		for (int i = 0; i < dataList.allTaxa.getTaxonCount(); i++) {

			heights[i] = (Double) dataList.allTaxa.getTaxon(i).getAttribute(
					Utils.ABSOLUTE_HEIGHT);

		}// END: taxon loop

	}// END: getHeights

	private void getTrees() {

		trees = new String[dataList.allTaxa.getTaxonCount()];
		for (int i = 0; i < dataList.allTaxa.getTaxonCount(); i++) {

			trees[i] = (String) dataList.allTaxa.getTaxon(i).getAttribute(
					Utils.TREE_FILENAME);

		}// END: taxon loop

	}// END: getHeights

	public String toString() {
		StringBuffer buffer = new StringBuffer();

		buffer.append(getColumnName(0));
		for (int j = 1; j < getColumnCount(); j++) {
			buffer.append("\t");
			buffer.append(getColumnName(j));
		}
		buffer.append("\n");

		for (int i = 0; i < getRowCount(); i++) {
			buffer.append(getValueAt(i, 0));
			for (int j = 1; j < getColumnCount(); j++) {
				buffer.append("\t");
				buffer.append(getValueAt(i, j));
			}
			buffer.append("\n");
		}

		return buffer.toString();
	}// END: toString

	public void fireTaxaChanged() {
		getHeights();
		getTrees();
		fireTableDataChanged();
	}// END: fireTaxaChanged

	public void setDataList(PartitionDataList dataList) {
		this.dataList = dataList;
	}// END: setDataList

}// END: TaxaTableModel class
