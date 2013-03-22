package dr.app.bss;

import java.awt.Component;

import javax.swing.JButton;
import javax.swing.JTable;
import javax.swing.table.TableCellRenderer;

@SuppressWarnings("serial")
public class JTableButtonCellRenderer extends JButton implements
		TableCellRenderer {

	public JTableButtonCellRenderer() {
		super();
		setOpaque(true);
	}// END: Constructor

	@Override
	public Component getTableCellRendererComponent(JTable table, Object value,
			boolean isSelected, boolean hasFocus, int row, int column) {

		setEnabled(table.isEnabled());

		JButton button = (JButton) value;
		setBackground(isSelected ? table.getSelectionBackground() : table
				.getBackground());
		
		//TODO: label on button
		// button.setText((value == null) ? "" : value.toString());

		return button;
	}// END: getTableCellRendererComponent

}// END: JTableButtonRenderer class
