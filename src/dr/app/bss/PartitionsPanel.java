package dr.app.bss;

import jam.framework.Exportable;

import java.awt.BorderLayout;
import java.awt.Component;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.io.File;

import javax.swing.AbstractAction;
import javax.swing.Action;
import javax.swing.DefaultCellEditor;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JComponent;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTable;
import javax.swing.event.TableModelEvent;
import javax.swing.event.TableModelListener;
import javax.swing.table.TableCellRenderer;
import javax.swing.table.TableColumn;

import org.virion.jam.panels.ActionPanel;

@SuppressWarnings("serial")
public class PartitionsPanel extends JPanel implements Exportable {

	private BeagleSequenceSimulatorFrame frame = null;
	private PartitionDataList dataList = null;

	private JTable partitionTable = null;
	private PartitionTableModel partitionTableModel = null;
	private JScrollPane scrollPane;

	private TableColumn column;

	private int partitionsCount = 1;

	private Action addPartitionAction = new AbstractAction("+") {
		public void actionPerformed(ActionEvent ae) {

			partitionTableModel.addDefaultRow();

			partitionsCount++;
			setPartitions();
		}// END: actionPerformed
	};

	private Action removePartitionAction = new AbstractAction("-") {
		public void actionPerformed(ActionEvent ae) {
			if (partitionsCount > 1) {
				partitionTableModel.deleteRow(partitionsCount - 1);
				partitionsCount--;
				setPartitions();
			}
		}// END: actionPerformed
	};

	public PartitionsPanel(final BeagleSequenceSimulatorFrame frame,
			final PartitionDataList dataList) {

		super();

		this.frame = frame;
		this.dataList = dataList;

		partitionTableModel = new PartitionTableModel(this.dataList);
		partitionTableModel
				.addTableModelListener(new PartitionTableModelListener());

		partitionTable = new JTable(partitionTableModel);
		partitionTable.setFillsViewportHeight(true);
		partitionTable.getTableHeader().setReorderingAllowed(false);
		partitionTable.addMouseListener(new JTableButtonMouseListener(
				partitionTable));

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
				PartitionTableModel.PARTITION_TREE_INDEX);
		column.setCellEditor(new JTableComboBoxCellEditor());
		column.setCellRenderer(new JTableComboBoxCellRenderer());

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

		ActionPanel actionPanel = new ActionPanel(false);
		actionPanel.setAddAction(addPartitionAction);
		actionPanel.setRemoveAction(removePartitionAction);
		add(actionPanel, BorderLayout.SOUTH);

		setPartitions();

	}// END: Constructor

	public JComponent getExportableComponent() {
		return this;
	}// END: getExportableComponent

	private void setPartitions() {
		addPartitionAction.setEnabled(true);
		if (partitionsCount == 1) {
			removePartitionAction.setEnabled(false);
		} else {
			removePartitionAction.setEnabled(true);
		}

		ColumnResizer.adjustColumnPreferredWidths(partitionTable);
		partitionTable.revalidate();
		this.updateUI();
	}// END: setPartitions

	// Listen to tree choices, set tree model in partition data
	private class PartitionTableModelListener implements TableModelListener {
		public void tableChanged(TableModelEvent ev) {

			if (ev.getType() == TableModelEvent.UPDATE) {
				int row = ev.getFirstRow();
				int column = ev.getColumn();

				if (column == PartitionTableModel.PARTITION_TREE_INDEX) {

					File value = (File) partitionTableModel.getValueAt(row, column);
					// System.out.println(value);
					dataList.get(row).treeModel = dataList.forestMap.get(value);

				}// END: column check
			}// END: event check

			frame.fireModelChanged();

		}// END: tableChanged
	}// END: InteractiveTableModelListener

	public class JTableComboBoxCellRenderer extends JComboBox implements
			TableCellRenderer {

		public JTableComboBoxCellRenderer() {
			super();
			setOpaque(true);
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

			if (column == PartitionTableModel.PARTITION_TREE_INDEX) {
				for (File file : dataList.forestMap.keySet()) {
					((JComboBox) editorComponent).addItem(file);
				}
			}

			((JComboBox) editorComponent).setSelectedItem(value);
			delegate.setValue(value);

			return editorComponent;
		}

	}// END: JTableComboBoxCellEditor class

	private class JTableButtonCellRenderer extends JButton implements
			TableCellRenderer {

		public JTableButtonCellRenderer() {
			super();
			setOpaque(true);
		}// END: Constructor

		@Override
		public Component getTableCellRendererComponent(JTable table,
				Object value, boolean isSelected, boolean hasFocus, int row,
				int column) {

			setEnabled(table.isEnabled());

			JButton button = (JButton) value;
			setBackground(isSelected ? table.getSelectionBackground() : table
					.getBackground());
			// button.setText((value == null) ? "" : value.toString());

			return button;
		}
	}// END: JTableButtonRenderer class

	private class JTableButtonCellEditor extends DefaultCellEditor {

		protected JButton button;
		private String label;
		private boolean isPushed;

		public JTableButtonCellEditor() {
			super(new JCheckBox());
			button = new JButton();
			button.setOpaque(true);
			button.addActionListener(new ActionListener() {
				public void actionPerformed(ActionEvent e) {
					fireEditingStopped();
				}
			});
		}

		public Component getTableCellEditorComponent(JTable table,
				Object value, boolean isSelected, int row, int column) {

			if (isSelected) {
				button.setForeground(table.getSelectionForeground());
				button.setBackground(table.getSelectionBackground());
			} else {
				button.setForeground(table.getForeground());
				button.setBackground(table.getBackground());
			}

			// TODO: set label from model name
			label = (value == null) ? "" : value.toString();
			button.setText(label);
			isPushed = true;

			return button;
		}

		public Object getCellEditorValue() {
			if (isPushed) {
				//
			}
			isPushed = false;
			return new String(label);
		}

		public boolean stopCellEditing() {
			isPushed = false;
			return super.stopCellEditing();
		}

		protected void fireEditingStopped() {
			super.fireEditingStopped();
		}

	}// END: JTableButtonCellEditor

	private class JTableButtonMouseListener extends MouseAdapter {

		private final JTable table;

		public JTableButtonMouseListener(JTable table) {
			this.table = table;
		}// END: Constructor

		public void mouseClicked(MouseEvent e) {

			int column = table.getColumnModel().getColumnIndexAtX(e.getX());
			int row = e.getY() / table.getRowHeight();

			if (row < table.getRowCount() && row >= 0 && column < table.getColumnCount() && column >= 0) {

				Object value = table.getValueAt(row, column);
				if (value instanceof JButton) {

					((JButton) value).doClick();

				}// END: JButton check

			}// END: placement check
		}// END: mouseClicked

	}// END: JTableButtonMouseListener class

}// END: class
