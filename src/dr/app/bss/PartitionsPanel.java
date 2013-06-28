package dr.app.bss;

import jam.framework.Exportable;
import jam.panels.ActionPanel;
import jam.table.TableRenderer;

import java.awt.BorderLayout;
import java.awt.Component;
import java.awt.Insets;
import java.awt.event.ActionEvent;

import javax.swing.AbstractAction;
import javax.swing.Action;
import javax.swing.DefaultCellEditor;
import javax.swing.DefaultListCellRenderer;
import javax.swing.JComboBox;
import javax.swing.JComponent;
import javax.swing.JList;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTable;
import javax.swing.SwingConstants;
import javax.swing.event.TableModelEvent;
import javax.swing.event.TableModelListener;
import javax.swing.table.TableCellRenderer;
import javax.swing.table.TableColumn;

/**
 * @author Filip Bielejec
 * @version $Id$
 */
@SuppressWarnings("serial")
public class PartitionsPanel extends JPanel implements Exportable {

	private PartitionDataList dataList = null;

	private JTable partitionTable = null;
	private PartitionTableModel partitionTableModel = null;
	
	private TableColumnHider hider;
	private JScrollPane scrollPane;
	private TableColumn column;
	private int partitionsCount;

	private Action addPartitionAction = new AbstractAction("+") {
		public void actionPerformed(ActionEvent ae) {

			// partitionTableModel.copyPreviousRow();
			partitionTableModel.addDefaultRow();
			setPartitions();
			
		}// END: actionPerformed
	};

	private Action removePartitionAction = new AbstractAction("-") {
		public void actionPerformed(ActionEvent ae) {
			if (partitionsCount > 1) {
				
				partitionTableModel.deleteRow(partitionsCount - 1);
				setPartitions();
			
			}
		}// END: actionPerformed
	};

	public PartitionsPanel(PartitionDataList dataList) {

		this.dataList = dataList;

		partitionTable = new JTable();
		partitionTable.getTableHeader().setReorderingAllowed(false);
		partitionTable.addMouseListener(new JTableButtonMouseListener(
				partitionTable));

		partitionTableModel = new PartitionTableModel(this.dataList);
		partitionTableModel
				.addTableModelListener(new PartitionTableModelListener());
		partitionTable.setModel(partitionTableModel);

		hider = new TableColumnHider(partitionTable);

		setLayout(new BorderLayout());

		scrollPane = new JScrollPane(partitionTable,
				JScrollPane.VERTICAL_SCROLLBAR_ALWAYS,
				JScrollPane.HORIZONTAL_SCROLLBAR_ALWAYS);
		RowNumberTable rowNumberTable = new RowNumberTable(partitionTable);
		scrollPane.setRowHeaderView(rowNumberTable);
		scrollPane.setCorner(JScrollPane.UPPER_LEFT_CORNER,
				rowNumberTable.getTableHeader());

		scrollPane.getViewport().setOpaque(false);

		add(scrollPane, BorderLayout.CENTER);

		column = partitionTable.getColumnModel().getColumn(
				PartitionTableModel.DATA_INDEX);
		column.setCellRenderer(new JTableComboBoxCellRenderer());
		column.setCellEditor(new JTableComboBoxCellEditor());
		column.setMinWidth(100);
		
//		column = partitionTable.getColumnModel().getColumn(
//				PartitionTableModel.DATA_TYPE_INDEX);
//		column.setCellEditor(new JTableComboBoxCellEditor());
//		column.setCellRenderer(new JTableComboBoxCellRenderer());

		column = partitionTable.getColumnModel().getColumn(
				PartitionTableModel.FROM_INDEX);
		column.setCellRenderer(
				new TableRenderer(SwingConstants.LEFT, new Insets(0, 2,
						0, 2)));
		column.setPreferredWidth(80);
		
		column = partitionTable.getColumnModel().getColumn(
				PartitionTableModel.TO_INDEX);
		column.setCellRenderer(
				new TableRenderer(SwingConstants.LEFT, new Insets(0, 2,
						0, 2)));
		column.setPreferredWidth(80);
		
		column = partitionTable.getColumnModel().getColumn(
				PartitionTableModel.EVERY_INDEX);
		column.setCellRenderer(
				new TableRenderer(SwingConstants.LEFT, new Insets(0, 2,
						0, 2)));
		column.setPreferredWidth(80);
		
		column = partitionTable.getColumnModel().getColumn(
				PartitionTableModel.DEMOGRAPHIC_MODEL_INDEX);
		column.setCellRenderer(new JTableButtonCellRenderer());
		column.setCellEditor(new JTableButtonCellEditor());
		
		column = partitionTable.getColumnModel().getColumn(
				PartitionTableModel.BRANCH_SUBSTITUTION_MODEL_INDEX);
		column.setCellRenderer(new JTableButtonCellRenderer());
		column.setCellEditor(new JTableButtonCellEditor());

		column = partitionTable.getColumnModel().getColumn(
				PartitionTableModel.SITE_RATE_MODEL_INDEX);
		column.setCellRenderer(new JTableButtonCellRenderer());
		column.setCellEditor(new JTableButtonCellEditor());

		column = partitionTable.getColumnModel().getColumn(
				PartitionTableModel.CLOCK_RATE_MODEL_INDEX);
		column.setCellRenderer(new JTableButtonCellRenderer());
		column.setCellEditor(new JTableButtonCellEditor());

		column = partitionTable.getColumnModel().getColumn(
				PartitionTableModel.FREQUENCY_MODEL_INDEX);
		column.setCellRenderer(new JTableButtonCellRenderer());
		column.setCellEditor(new JTableButtonCellEditor());

		column = partitionTable.getColumnModel().getColumn(
				PartitionTableModel.ANCESTRAL_SEQUENCE_INDEX);
		column.setCellRenderer(new JTableButtonCellRenderer());
		column.setCellEditor(new JTableButtonCellEditor());

		ActionPanel actionPanel = new ActionPanel(false);
		actionPanel.setAddAction(addPartitionAction);
		actionPanel.setRemoveAction(removePartitionAction);
		add(actionPanel, BorderLayout.SOUTH);

		setPartitions();

	}// END: Constructor

	private void setPartitions() {

		partitionsCount = dataList.size();

		addPartitionAction.setEnabled(true);
		if (partitionsCount == 1) {
			removePartitionAction.setEnabled(false);
		} else {
			removePartitionAction.setEnabled(true);
		}

		ColumnResizer.adjustColumnPreferredWidths(partitionTable);
	}// END: setPartitions

	// Listen to tree choices, set tree model in partition data
	private class PartitionTableModelListener implements TableModelListener {

		public void tableChanged(TableModelEvent ev) {

			if (ev.getType() == TableModelEvent.UPDATE) {
				int row = ev.getFirstRow();
				int column = ev.getColumn();

				if (column == PartitionTableModel.DATA_INDEX) {

					TreesTableRecord value = (TreesTableRecord) partitionTableModel.getValueAt(row,
							column);
					
					dataList.get(row).record = (TreesTableRecord) value;

				}// END: column check

			}// END: event check

			// frame.collectAllSettings();

		}// END: tableChanged
	}// END: InteractiveTableModelListener

	public class JTableComboBoxCellRenderer extends JComboBox implements
			TableCellRenderer {

		private DefaultListCellRenderer comboBoxRenderer = new DefaultListCellRenderer() {

			@Override
			public Component getListCellRendererComponent(JList list,
					Object value, int index, boolean isSelected,
					boolean cellHasFocus) {

				 super.getListCellRendererComponent(list, value, index, isSelected, cellHasFocus);
				
				if (value != null) {
					
					TreesTableRecord record = (TreesTableRecord) value;
					this.setText(record.getName());

				}

				return this;
			}

		};

		public JTableComboBoxCellRenderer() {

			super();
			setOpaque(true);
			this.setRenderer(comboBoxRenderer);

		}// END: Constructor

		public Component getTableCellRendererComponent(JTable table,
				Object value, boolean isSelected, boolean hasFocus, int row,
				int column) {

			if (isSelected) {

				this.setForeground(table.getSelectionForeground());
				this.setBackground(table.getSelectionBackground());

			} else {

				this.setForeground(table.getForeground());
				this.setBackground(table.getBackground());

			}

			// Select the current value
			setSelectedItem(value);

			if (value != null) {
				removeAllItems();
				addItem(value);
			}

			return this;
		}
	}// END: JTableComboBoxCellRenderer class
	
	private class JTableComboBoxCellEditor extends DefaultCellEditor {

		public JTableComboBoxCellEditor() {
			super(new JComboBox());
		}

		public Component getTableCellEditorComponent(JTable table,
				Object value, boolean isSelected, int row, int column) {

			((JComboBox) editorComponent).removeAllItems();

			if (column == PartitionTableModel.DATA_INDEX) {

				for (TreesTableRecord record : dataList.recordsList) {
					((JComboBox) editorComponent).addItem(record);
				}// END: fill loop
				
//			} else if (column == PartitionTableModel.DATA_TYPE_INDEX) {
//
//				for (String dataType : PartitionData.dataTypes) {
//					((JComboBox) editorComponent).addItem(dataType);
//				}// END: fill loop

			} else {

				// do nothing

			}// END: column check

			((JComboBox) editorComponent).setSelectedItem(value);
			delegate.setValue(value);

			return editorComponent;
		}// END: getTableCellEditorComponent

	}// END: JTableComboBoxCellEditor class

	public void hideTreeColumn() {
		hider.hide(PartitionTableModel.COLUMN_NAMES[PartitionTableModel.DATA_INDEX]);
	}

	public void showTreeColumn() {
		hider.show(PartitionTableModel.COLUMN_NAMES[PartitionTableModel.DATA_INDEX]);
	}

	public JComponent getExportableComponent() {
		return this;
	}// END: getExportableComponent

	public void updatePartitionTable(PartitionDataList dataList) {
		partitionTableModel.setDataList(dataList);
		setDataList(dataList);
		setPartitions();
		partitionTableModel.fireTableDataChanged();
	}// END: updatePartitionTable

	public void setDataList(PartitionDataList dataList) {
		this.dataList = dataList;
	}

}// END: class
