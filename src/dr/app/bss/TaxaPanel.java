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

	public void updateTaxaTable(PartitionDataList dataList) {
		taxaTableModel.setDataList(dataList);
		setDataList(dataList);
		fireTaxaChanged();
	}// END: updateTaxaTable

	public void fireTaxaChanged() {
		taxaTableModel.fireTaxaChanged();
	}// END: fireTableDataChanged

	public void setDataList(PartitionDataList dataList) {
		this.dataList = dataList;
	}

}// END: class
