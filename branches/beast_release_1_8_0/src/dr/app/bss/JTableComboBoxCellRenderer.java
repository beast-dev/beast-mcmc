package dr.app.bss;

import java.awt.Component;

import javax.swing.DefaultListCellRenderer;
import javax.swing.JComboBox;
import javax.swing.JList;
import javax.swing.JTable;
import javax.swing.table.TableCellRenderer;


@SuppressWarnings("serial")
public class JTableComboBoxCellRenderer extends JComboBox implements
		TableCellRenderer {

	private int columnIndex;
	
	private DefaultListCellRenderer comboBoxRenderer = new DefaultListCellRenderer() {

		@Override
		public Component getListCellRendererComponent(JList list, Object value,
				int index, boolean isSelected, boolean cellHasFocus) {

			super.getListCellRendererComponent(list, value, index, isSelected,
					cellHasFocus);

			if (value != null) {

				switch (columnIndex) {

				case PartitionTableModel.DATA_INDEX:

					TreesTableRecord record = (TreesTableRecord) value;
					this.setText(record.getName());
					break;

				default:
					break;

				}// END: switch				
				
			}//END: null check

			return this;
		}

	};

	public JTableComboBoxCellRenderer(int columnIndex) {

		super();
		setOpaque(true);
		this.setRenderer(comboBoxRenderer);
		this.columnIndex = columnIndex;
		
	}// END: Constructor

	public Component getTableCellRendererComponent(JTable table, Object value,
			boolean isSelected, boolean hasFocus, int row, int column) {

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