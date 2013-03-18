package dr.app.bss;

import jam.framework.Exportable;
import jam.table.HeaderRenderer;
import jam.table.TableRenderer;

import java.awt.BorderLayout;
import java.awt.Insets;

import javax.swing.JComponent;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTable;
import javax.swing.SwingConstants;
import javax.swing.plaf.BorderUIResource;

import dr.app.gui.table.TableEditorStopper;

@SuppressWarnings("serial")
public class TaxaPanel extends JPanel implements Exportable {

	private PartitionDataList dataList = null;

	private JScrollPane scrollPane = new JScrollPane();
	private JTable taxaTable = null;
	private TaxaTableModel taxaTableModel = null;

//	private double[] heights = null;

	public TaxaPanel(PartitionDataList dataList) {

		this.dataList = dataList;

		taxaTable = new JTable();

		taxaTableModel = new TaxaTableModel(this.dataList);
		taxaTable.setModel(taxaTableModel);

		taxaTable.getTableHeader().setReorderingAllowed(false);

		taxaTable.getTableHeader()
				.setDefaultRenderer(
						new HeaderRenderer(SwingConstants.LEFT, new Insets(0,
								2, 0, 2)));

		taxaTable
				.getColumnModel()
				.getColumn(0)
				.setCellRenderer(
						new TableRenderer(SwingConstants.LEFT, new Insets(0, 2,
								0, 2)));

		taxaTable.getColumnModel().getColumn(0).setPreferredWidth(80);

		taxaTable
				.getColumnModel()
				.getColumn(1)
				.setCellRenderer(
						new TableRenderer(SwingConstants.LEFT, new Insets(0, 2,
								0, 2)));

		taxaTable.getColumnModel().getColumn(1).setPreferredWidth(80);

		TableEditorStopper.ensureEditingStopWhenTableLosesFocus(taxaTable);

		scrollPane = new JScrollPane(taxaTable,
				JScrollPane.VERTICAL_SCROLLBAR_ALWAYS,
				JScrollPane.HORIZONTAL_SCROLLBAR_ALWAYS);

		RowNumberTable rowNumberTable = new RowNumberTable(taxaTable);
		scrollPane.setRowHeaderView(rowNumberTable);
		scrollPane.setCorner(JScrollPane.UPPER_LEFT_CORNER,
				rowNumberTable.getTableHeader());

		scrollPane.getViewport().setOpaque(false);

		setOpaque(false);
		setBorder(new BorderUIResource.EmptyBorderUIResource(
				new java.awt.Insets(12, 12, 12, 12)));
		setLayout(new BorderLayout(0, 0));

		add(scrollPane, "Center");

	}// END: Constructor

	public JComponent getExportableComponent() {
		return taxaTable;
	}

//	private void getHeights() {
//
//		heights = new double[dataList.taxonList.getTaxonCount()];
//		for (int i = 0; i < dataList.taxonList.getTaxonCount(); i++) {
//
//			heights[i] = (Double) dataList.taxonList.getTaxon(i).getAttribute(
//					Utils.ABSOLUTE_HEIGHT);
//
//		}// END: taxon loop
//
//	}// END: getHeights

//	private class TaxaTableModel extends AbstractTableModel {
//
//		private PartitionDataList dataList;
//		private String[] columnNames = { "Name", "Height" };
//
//		public TaxaTableModel(PartitionDataList dataList) {
//			this.dataList = dataList;
//		}
//
//		public void setDataList(PartitionDataList dataList) {
//			this.dataList = dataList;
//		}
//
//		public int getColumnCount() {
//			return columnNames.length;
//		}
//
//		public int getRowCount() {
//			return this.dataList.taxonList.getTaxonCount();
//		}
//
//		public Object getValueAt(int row, int col) {
//			switch (col) {
//
//			case 0:
//				return this.dataList.taxonList.getTaxonId(row);
//
//			case 1:
//
//				if (heights != null) {
//					return heights[row];
//				} else {
//					return 0.0;
//				}
//
//			default:
//				return null;
//
//			}// END: switch
//		}// END: getValueAt
//
//		public void setValueAt(Object value, int row, int col) {
//
//			switch (col) {
//
//			case 0:
//				this.dataList.taxonList.getTaxon(row).setId(value.toString());
//				break;
//
//			case 1:
//				// dataList.get(0).taxonList.getTaxon(row).getHeight();
//				break;
//
//			default:
//				break;
//
//			}// END: switch
//		}// END: setValueAt
//
//		public boolean isCellEditable(int row, int col) {
//			return false;
//		}
//
//		public String getColumnName(int column) {
//			return columnNames[column];
//		}
//
//		public Class<? extends Object> getColumnClass(int c) {
//			return getValueAt(0, c).getClass();
//		}
//
//		public String toString() {
//			StringBuffer buffer = new StringBuffer();
//
//			buffer.append(getColumnName(0));
//			for (int j = 1; j < getColumnCount(); j++) {
//				buffer.append("\t");
//				buffer.append(getColumnName(j));
//			}
//			buffer.append("\n");
//
//			for (int i = 0; i < getRowCount(); i++) {
//				buffer.append(getValueAt(i, 0));
//				for (int j = 1; j < getColumnCount(); j++) {
//					buffer.append("\t");
//					buffer.append(getValueAt(i, j));
//				}
//				buffer.append("\n");
//			}
//
//			return buffer.toString();
//		}
//
//	}// END: TaxaTableModel class

	public void updateTaxaTable(PartitionDataList dataList) {
		taxaTableModel.setDataList(dataList);
		setDataList(dataList);
//		taxaTableModel.
		fireTableDataChanged();
	}// END: updateTaxaTable

	public void fireTableDataChanged() {
//		getHeights();
		taxaTableModel.fireTaxaChanged();
	}// END: fireTableDataChanged

	public void setDataList(PartitionDataList dataList) {
		this.dataList = dataList;
	}

}// END: class
