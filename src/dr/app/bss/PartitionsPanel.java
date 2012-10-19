package dr.app.bss;

import jam.framework.Exportable;

import java.awt.BorderLayout;
import java.awt.Component;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.util.ArrayList;

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
	private ArrayList<PartitionData> dataList = null;

	private JTable partitionTable = null;
	private PartitionTableModel partitionTableModel = null;
	private JScrollPane scrollPane;

	private TableColumn column;

	private int partitionsCount = 1;
	
	private Action addPartitionAction = new AbstractAction("+") {
		public void actionPerformed(ActionEvent ae) {
			partitionsCount++;
			setPartitions();
		}// END: actionPerformed
	};

	private Action removePartitionAction = new AbstractAction("-") {
		public void actionPerformed(ActionEvent ae) {
			if (partitionsCount > 1) {
				partitionsCount--;
				setPartitions();
			}
		}// END: actionPerformed
	};
	
	public PartitionsPanel(final BeagleSequenceSimulatorFrame frame,
			final ArrayList<PartitionData> dataList) {

		super();

		this.frame = frame;
		this.dataList = dataList;

		partitionTableModel = new PartitionTableModel(this.dataList);
		partitionTableModel
				.addTableModelListener(new PartitionTableModelListener());

		partitionTable = new JTable(partitionTableModel);

		partitionTable.setFillsViewportHeight(true);
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

		// set columns
		column = partitionTable.getColumnModel().getColumn(
				partitionTableModel.PARTITION_TREE_INDEX);
		column.setCellRenderer(new JTableComboBoxCellRenderer(
				new String[] { "" }));
		column.setCellEditor(new JTableComboBoxCellEditor(new String[] { "" }));

		column = partitionTable.getColumnModel().getColumn(
				partitionTableModel.BRANCH_SUBSTITUTION_MODEL_INDEX);
		column.setCellRenderer(new JTableButtonCellRenderer());
		column.setCellEditor(new JTableButtonCellEditor(new JCheckBox()));

		column = partitionTable.getColumnModel().getColumn(
				partitionTableModel.SITE_RATE_MODEL_INDEX);
		column.setCellRenderer(new JTableButtonCellRenderer());
		column.setCellEditor(new JTableButtonCellEditor(new JCheckBox()));

		column = partitionTable.getColumnModel().getColumn(
				partitionTableModel.CLOCK_RATE_MODEL_INDEX);
		column.setCellRenderer(new JTableButtonCellRenderer());
		column.setCellEditor(new JTableButtonCellEditor(new JCheckBox()));

		column = partitionTable.getColumnModel().getColumn(
				partitionTableModel.FREQUENCY_MODEL_INDEX);
		column.setCellRenderer(new JTableButtonCellRenderer());
		column.setCellEditor(new JTableButtonCellEditor(new JCheckBox()));

		ActionPanel actionPanel = new ActionPanel(false);
		actionPanel.setAddAction(addPartitionAction);
		actionPanel.setRemoveAction(removePartitionAction);
		add(actionPanel, BorderLayout.SOUTH);
		 
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
        
		this.updateUI();
	}//END: setPartitions
	
	private class PartitionTableModelListener implements TableModelListener {
		public void tableChanged(TableModelEvent ev) {

			System.out.println("TODO: PartitionTableModelListener");

		}
	}// END: InteractiveTableModelListener

	private class JTableComboBoxCellRenderer extends JComboBox implements
			TableCellRenderer {

		public JTableComboBoxCellRenderer(String[] items) {
			super(items);
		}

		public Component getTableCellRendererComponent(JTable table,
				Object value, boolean isSelected, boolean hasFocus, int row,
				int column) {

			if (isSelected) {

				setForeground(table.getSelectionForeground());
				super.setBackground(table.getSelectionBackground());

			} else {

				setForeground(table.getForeground());
				setBackground(table.getBackground());
			}

			// Select the current value
			setSelectedItem(value);

			return this;
		}
	}// END: JTableComboBoxCellRenderer class

	private class JTableComboBoxCellEditor extends DefaultCellEditor {
		public JTableComboBoxCellEditor(String[] items) {
			super(new JComboBox(items));
		}
	}// END: JTableComboBoxCellEditor class

	private class JTableButtonCellRenderer extends JButton implements
			TableCellRenderer {

		public JTableButtonCellRenderer() {
			super();
		}

		@Override
		public Component getTableCellRendererComponent(JTable table,
				Object value, boolean isSelected, boolean hasFocus, int row,
				int column) {

			JButton button = (JButton) value;
			setBackground(isSelected ? table.getSelectionBackground() : table
					.getBackground());
			// setText((value == null) ? "" : value.toString());
			return button;
		}
	}// END: JTableButtonRenderer

	private class JTableButtonCellEditor extends DefaultCellEditor {
		protected JButton button;

		private String label;

		private boolean isPushed;

		public JTableButtonCellEditor(JCheckBox checkBox) {
			super(checkBox);
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
		}

		public void mouseClicked(MouseEvent e) {
			int column = table.getColumnModel().getColumnIndexAtX(e.getX());
			int row = e.getY() / table.getRowHeight();

			if (row < table.getRowCount() && row >= 0
					&& column < table.getColumnCount() && column >= 0) {
				Object value = table.getValueAt(row, column);
				if (value instanceof JButton) {
					((JButton) value).doClick();
				}
			}
		}
	}// END: JTableButtonMouseListener

	//

}// END: class
