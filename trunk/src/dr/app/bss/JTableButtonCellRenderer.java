package dr.app.bss;

import java.awt.Color;
import java.awt.Component;

import javax.swing.JButton;
import javax.swing.JTable;
import javax.swing.UIManager;
import javax.swing.border.LineBorder;
import javax.swing.table.TableCellRenderer;

/**
 * @author Filip Bielejec
 * @version $Id$
 */
@SuppressWarnings("serial")
public class JTableButtonCellRenderer extends JButton implements
		TableCellRenderer {

	private PartitionDataList dataList = null;
	
	public JTableButtonCellRenderer() {
		super();
		setOpaque(true);
	}// END: Constructor

	public JTableButtonCellRenderer(PartitionDataList dataList) {
		super();
		setOpaque(true);
		this.dataList = dataList;
	}// END: Constructor
	
	@Override
	public Component getTableCellRendererComponent(JTable table, Object value,
			boolean isSelected, boolean hasFocus, int row, int column) {

		JButton button = (JButton) value;
		
		if (isSelected) {
			button.setForeground(table.getSelectionForeground());
			button.setBackground(table.getSelectionBackground());
		} else {
			button.setForeground(table.getForeground());
			button.setBackground(UIManager.getColor("Button.background"));
		}

		if (hasFocus) {
			button.setBorder( new LineBorder(Color.BLUE) );
		} else {
			button.setBorder( button.getBorder() );
		}
		
		// label on button
		if (dataList != null) {

			String label;
			if (dataList.recordsList.get(row).getName().equalsIgnoreCase("")) {
				
				label = Utils.CHOOSE_FILE;
				
			} else {
				
				label = dataList.recordsList.get(row).getName();
				
			}

			button.setText(label);
		}//END: null check

		return button;
	}// END: getTableCellRendererComponent

}// END: JTableButtonRenderer class
