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
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;
import javax.swing.plaf.BorderUIResource;
import javax.swing.table.AbstractTableModel;

import dr.app.gui.table.DateCellEditor;
import dr.app.gui.table.TableEditorStopper;
import dr.app.gui.table.TableSorter;

@SuppressWarnings("serial")
public class TaxaPanel extends JPanel implements Exportable {

	private BeagleSequenceSimulatorFrame frame = null;
	private PartitionDataList dataList = null;

	private JScrollPane scrollPane = new JScrollPane();
	private JTable dataTable = null;
	private TaxaTableModel taxaTableModel = null;

	private double[] heights = null;

	public TaxaPanel(BeagleSequenceSimulatorFrame frame,
			PartitionDataList dataList) {

		this.frame = frame;
		this.dataList = dataList;

		taxaTableModel = new TaxaTableModel();
		TableSorter sorter = new TableSorter(taxaTableModel);
		dataTable = new JTable(sorter);

		sorter.setTableHeader(dataTable.getTableHeader());

		dataTable.getTableHeader().setReorderingAllowed(false);
		dataTable.getTableHeader()
				.setDefaultRenderer(
						new HeaderRenderer(SwingConstants.LEFT, new Insets(0,
								4, 0, 4)));

		dataTable
				.getColumnModel()
				.getColumn(0)
				.setCellRenderer(
						new TableRenderer(SwingConstants.LEFT, new Insets(0, 4,
								0, 4)));
		dataTable.getColumnModel().getColumn(0).setPreferredWidth(80);

		dataTable
				.getColumnModel()
				.getColumn(1)
				.setCellRenderer(
						new TableRenderer(SwingConstants.LEFT, new Insets(0, 4,
								0, 4)));
		dataTable.getColumnModel().getColumn(1).setPreferredWidth(80);
		dataTable.getColumnModel().getColumn(1)
				.setCellEditor(new DateCellEditor());

		TableEditorStopper.ensureEditingStopWhenTableLosesFocus(dataTable);

		dataTable.getSelectionModel().addListSelectionListener(
				new ListSelectionListener() {
					public void valueChanged(ListSelectionEvent evt) {
						selectionChanged();
					}
				});

		scrollPane = new JScrollPane(dataTable,
				JScrollPane.VERTICAL_SCROLLBAR_ALWAYS,
				JScrollPane.HORIZONTAL_SCROLLBAR_ALWAYS);
		scrollPane.setOpaque(false);

		setOpaque(false);
		setBorder(new BorderUIResource.EmptyBorderUIResource(
				new java.awt.Insets(12, 12, 12, 12)));
		setLayout(new BorderLayout(0, 0));

		add(scrollPane, "Center");

	}

	public JComponent getExportableComponent() {
		return dataTable;
	}

	public void selectionChanged() {

		int[] selRows = dataTable.getSelectedRows();
		if (selRows == null || selRows.length == 0) {
			frame.dataSelectionChanged(false);
		} else {
			frame.dataSelectionChanged(true);
		}
	}

	// TODO: multiple partitions
	private void getHeights() {

		System.out.println("Parsing heights");

		heights = new double[dataList.taxonList.getTaxonCount()];;
		for (int i = 0; i < dataList.taxonList.getTaxonCount(); i++) {

			heights[i] = dataList.taxonList.getTaxon(i).getHeight();

		}// END: taxon loop
		
		taxaTableModel.fireTableDataChanged();
	}//END: getHeights

	private class TaxaTableModel extends AbstractTableModel {

		String[] columnNames = { "Name", "Height" };

		public TaxaTableModel() {
		}

		public int getColumnCount() {
			return columnNames.length;
		}

		public int getRowCount() {
			return dataList.taxonList.getTaxonCount();
		}

		public Object getValueAt(int row, int col) {
			switch (col) {
			case 0:
				return dataList.taxonList.getTaxonId(row);
			case 1:
				
//				getHeights();
				
				if (heights != null) {
					return heights[row];
				} else {
					
					return "0.0";
					
				}
				default:
					return null;
			}
		}//END: getValueAt

		public void setValueAt(Object value, int row, int col) {
			
			getHeights();
			
			if (col == 0) {
				dataList.taxonList.getTaxon(row).setId(value.toString());
			}
			
//			if (col == 1) {
//				
//				dataList.get(0).taxonList.getTaxon(row).getHeight();
//				
//			}
			
		}

		public boolean isCellEditable(int row, int col) {
			return false;
		}

		public String getColumnName(int column) {
			return columnNames[column];
		}

		public Class<? extends Object> getColumnClass(int c) {
			return getValueAt(0, c).getClass();
		}

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
		}

	}// END: TaxaTableModel class

}// END: class
