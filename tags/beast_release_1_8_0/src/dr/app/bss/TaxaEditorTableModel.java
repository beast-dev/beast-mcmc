package dr.app.bss;

import javax.swing.table.AbstractTableModel;

import dr.evolution.util.Taxa;
import dr.evolution.util.Taxon;

@SuppressWarnings("serial")
public class TaxaEditorTableModel extends AbstractTableModel {

	private String[] COLUMN_NAMES = { "Name", "Height" };
	public final static int NAME_INDEX = 0;
	public final static int HEIGHT_INDEX = 1;

	private Taxa taxaSet;

	public TaxaEditorTableModel() {

		taxaSet = new Taxa();

	}// END: Constructor

	public void addEmptyRow() {
		String name = "";
		Taxon taxon = new Taxon(name);
		taxon.setAttribute(Utils.ABSOLUTE_HEIGHT, (Double) 0.0);
		taxaSet.addTaxon(taxon);
		this.fireTableDataChanged();
	}// END: addEmptyRow

	public void removeLastRow() {
		int lastIndex = taxaSet.getTaxonCount() - 1;
		Taxon taxon = taxaSet.getTaxon(lastIndex);
		taxaSet.removeTaxon(taxon);
		this.fireTableDataChanged();
	}// END: removeLastRow

	public String getColumnName(int column) {
		return COLUMN_NAMES[column];
	}

	@Override
	public int getColumnCount() {
		return COLUMN_NAMES.length;
	}// END: getColumnCount

	@Override
	public int getRowCount() {
		return taxaSet.getTaxonCount();
	}// END: getRowCount

	public boolean isCellEditable(int row, int col) {
		return true;
	}

	public Class<? extends Object> getColumnClass(int c) {
		return getValueAt(0, c).getClass();
	}

	@Override
	public Object getValueAt(int row, int col) {

		switch (col) {

		case NAME_INDEX:

			return taxaSet.getTaxonId(row);

		case HEIGHT_INDEX:

			return taxaSet.getTaxon(row).getAttribute(Utils.ABSOLUTE_HEIGHT);

		default:
			return null;

		}// END: switch

	}// END: getValueAt

	public void setValueAt(Object value, int row, int col) {

		switch (col) {

		case NAME_INDEX:
			taxaSet.setTaxonId(row, (String) value);
			break;

		case HEIGHT_INDEX:
			taxaSet.setTaxonAttribute(row, Utils.ABSOLUTE_HEIGHT,
					(Double) value);
			break;

		default:
			break;

		}// END: switch
	}// END: setValueAt

	public String toString(boolean header) {
		StringBuffer buffer = new StringBuffer();

		if (header) {
			buffer.append(getColumnName(0));
			for (int j = 1; j < getColumnCount(); j++) {
				buffer.append("\t");
				buffer.append(getColumnName(j));
			}
			buffer.append("\n");
		}

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

	public void fireTaxonListChanged() {
		this.fireTableDataChanged();
	}// END: fireTaxaChanged

	public void setTaxaSet(Taxa taxaSet) {
		this.taxaSet = taxaSet;
	}// END: setTaxonList

	public Taxa getTaxaSet() {
		return taxaSet;
	}// END: getTaxonList

}// END: TaxaEditorTableModel class