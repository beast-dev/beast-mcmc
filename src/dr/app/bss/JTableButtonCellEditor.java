package dr.app.bss;

import java.awt.Component;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.DefaultCellEditor;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JTable;

/**
 * @author Filip Bielejec
 * @version $Id$
 */
@SuppressWarnings("serial")
public class JTableButtonCellEditor extends DefaultCellEditor {

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

	}// END: Constructor

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

}// END: class
