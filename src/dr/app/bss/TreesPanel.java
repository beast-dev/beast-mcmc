package dr.app.bss;

import java.awt.BorderLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.io.File;

import jam.framework.Exportable;
import jam.panels.ActionPanel;
import jam.table.TableRenderer;

import javax.swing.AbstractAction;
import javax.swing.Action;
import javax.swing.JComponent;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTable;
import javax.swing.SwingConstants;
import javax.swing.event.TableModelEvent;
import javax.swing.event.TableModelListener;
import javax.swing.table.TableColumn;




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

				treesTableModel.deleteRow(treesCount - 1);
				setTrees();

			}
		}// END: actionPerformed
	};
	
	public TreesPanel(MainFrame frame,
			PartitionDataList dataList) {
		
		this.frame = frame;
		this.dataList = dataList;
		
		treesTable = new JTable();
		treesTable.getTableHeader().setReorderingAllowed(false);
		treesTable.addMouseListener(new JTableButtonMouseListener(
				treesTable));
		
		treesTableModel = new TreesTableModel(this.dataList);
		treesTableModel
				.addTableModelListener(new TreesTableModelListener());
		treesTable.setModel(treesTableModel);
		
		column = treesTable.getColumnModel().getColumn(
				TreesTableModel.TREE_FILE_INDEX);
		column.setCellRenderer(new JTableButtonCellRenderer());
		column.setCellEditor(new JTableButtonCellEditor());

		column = treesTable.getColumnModel().getColumn(
				TreesTableModel.TAXA_INDEX);
		column.setCellRenderer(
				new TableRenderer(SwingConstants.LEFT, new Insets(0, 2,
						0, 2)));
		
		column = treesTable.getColumnModel().getColumn(
				TreesTableModel.TREES_INDEX);
		column.setCellRenderer(
				new TableRenderer(SwingConstants.LEFT, new Insets(0, 2,
						0, 2)));
		
		
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
		
		ActionPanel actionPanel = new ActionPanel(false);
		actionPanel.setAddAction(addTreeAction);
		actionPanel.setRemoveAction(removeTreeAction);
		add(actionPanel, BorderLayout.SOUTH);
		
		setTrees();
		
	}//END: Constructor
	
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
	
	// Listen to tree choices, set tree model in partition data
	private class TreesTableModelListener implements TableModelListener {

		public void tableChanged(TableModelEvent ev) {

			if (ev.getType() == TableModelEvent.UPDATE) {
				int row = ev.getFirstRow();
				int column = ev.getColumn();

				if (column == TreesTableModel.TREE_FILE_INDEX) {

					File value = (File) treesTableModel.getValueAt(row,
							column);
					
					//TODO
//					dataList.get(row).treeFile = value;
//					dataList.treeFileList = dataList.treeFileList.add(value);
					
				}

			}// END: event check

		}// END: tableChanged
	}// END: InteractiveTableModelListener
	
	@Override
	public JComponent getExportableComponent() {
		return this;
	}// END: getExportableComponent

}// END: class
