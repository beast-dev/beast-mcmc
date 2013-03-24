package dr.app.bss;

import jam.framework.Exportable;
import jam.panels.ActionPanel;
import jam.table.TableRenderer;

import java.awt.BorderLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;

import javax.swing.AbstractAction;
import javax.swing.Action;
import javax.swing.JComponent;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTable;
import javax.swing.SwingConstants;
import javax.swing.table.TableColumn;

import dr.evolution.util.Taxon;

@SuppressWarnings("serial")
public class TreesPanel extends JPanel implements Exportable {

	private PartitionDataList dataList;
	private MainFrame frame;

	private JTable treesTable = null;
	private TreesTableModel treesTableModel = null;

	private JScrollPane scrollPane;
	private TableColumn column;
	private int treesCount;

	private Action addTreeAction = new AbstractAction("+") {
		public void actionPerformed(ActionEvent ae) {

			treesTableModel.addDefaultRow();
			setTrees();

		}// END: actionPerformed
	};

	private Action removeTreeAction = new AbstractAction("-") {
		public void actionPerformed(ActionEvent ae) {
			if (treesCount > 1) {

				String value = dataList.treeFileList.get(treesCount - 1).getName();
				removeTaxonWithAttributeValue(dataList, Utils.TREE_FILENAME, value);
				
				treesTableModel.deleteRow(treesCount - 1);
				
				frame.fireTaxaChanged();
				setTrees();
				
			}
		}// END: actionPerformed
	};

	public TreesPanel(MainFrame frame, PartitionDataList dataList) {

		this.frame = frame;
		this.dataList = dataList;

		treesTable = new JTable();
		treesTable.getTableHeader().setReorderingAllowed(false);
		treesTable.addMouseListener(new JTableButtonMouseListener(treesTable));

		treesTableModel = new TreesTableModel(this.dataList, this.frame);
		treesTable.setModel(treesTableModel);

		setLayout(new BorderLayout());

		scrollPane = new JScrollPane(treesTable,
				JScrollPane.VERTICAL_SCROLLBAR_ALWAYS,
				JScrollPane.HORIZONTAL_SCROLLBAR_ALWAYS);
		RowNumberTable rowNumberTable = new RowNumberTable(treesTable);
		scrollPane.setRowHeaderView(rowNumberTable);
		scrollPane.setCorner(JScrollPane.UPPER_LEFT_CORNER,
				rowNumberTable.getTableHeader());

		scrollPane.getViewport().setOpaque(false);

		add(scrollPane, BorderLayout.CENTER);

		column = treesTable.getColumnModel().getColumn(
				TreesTableModel.TREE_FILE_INDEX);
		// pass dataList for labels on buttons
		column.setCellRenderer(new JTableButtonCellRenderer(this.dataList));
		column.setCellEditor(new JTableButtonCellEditor());

		column = treesTable.getColumnModel().getColumn(
				TreesTableModel.TAXA_INDEX);
		column.setCellRenderer(new TableRenderer(SwingConstants.LEFT,
				new Insets(0, 2, 0, 2)));

		column = treesTable.getColumnModel().getColumn(
				TreesTableModel.TREES_INDEX);
		column.setCellRenderer(new TableRenderer(SwingConstants.LEFT,
				new Insets(0, 2, 0, 2)));

		ActionPanel actionPanel = new ActionPanel(false);
		actionPanel.setAddAction(addTreeAction);
		actionPanel.setRemoveAction(removeTreeAction);
		add(actionPanel, BorderLayout.SOUTH);

		setTrees();

	}// END: Constructor

	private void setTrees() {

		treesCount = dataList.treeFileList.size();

		addTreeAction.setEnabled(true);
		if (treesCount == 1) {
			removeTreeAction.setEnabled(false);
		} else {
			removeTreeAction.setEnabled(true);
		}

		ColumnResizer.adjustColumnPreferredWidths(treesTable);
	}// END: setPartitions

	private void removeTaxonWithAttributeValue(PartitionDataList dataList,
			String attribute, String value) {

		synchronized (dataList.taxonList) {
			for (int i = 0; i < dataList.taxonList.getTaxonCount(); i++) {

				Taxon taxon = dataList.taxonList.getTaxon(i);
				if (taxon.getAttribute(attribute).toString()
						.equalsIgnoreCase(value)) {
					dataList.taxonList.removeTaxon(taxon);
					i--;
				}
			}
		}

	}// END: removeTaxaWithAttributeValue
	
	@Override
	public JComponent getExportableComponent() {
		return this;
	}// END: getExportableComponent

}// END: class
