package dr.app.bss;

import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;

import javax.swing.JButton;
import javax.swing.JTable;

/**
 * @author Filip Bielejec
 * @version $Id$
 */
public class JTableButtonMouseListener extends MouseAdapter {

	private final JTable table;

	public JTableButtonMouseListener(JTable table) {
		this.table = table;
	}// END: Constructor

	public void mouseClicked(MouseEvent e) {

		int column = table.getColumnModel().getColumnIndexAtX(e.getX());
		int row = e.getY() / table.getRowHeight();

		if (row < table.getRowCount() && row >= 0
				&& column < table.getColumnCount() && column >= 0) {

			Object value = table.getValueAt(row, column);
			if (value instanceof JButton) {

				((JButton) value).doClick();

			}// END: JButton check

		}// END: placement check
	}// END: mouseClicked

}// END: class