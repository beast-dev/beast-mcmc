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
	private int columnIndex;
	private String toolTipText = null;
    
	public JTableButtonCellRenderer() {
		super();
		setOpaque(true);
	}// END: Constructor

	public JTableButtonCellRenderer(PartitionDataList dataList, int columnIndex) {
		super();

		setOpaque(true);
		this.dataList = dataList;
		this.columnIndex = columnIndex;

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
			button.setBorder(new LineBorder(Color.BLUE));
		} else {
			button.setBorder(button.getBorder());
		}

		// label on button
		if (dataList != null) {

			// switch between edit taxa and choose file labels
			String label = null;
			switch (columnIndex) {

			case TreesTableModel.TREE_FILE_INDEX:

				if (!dataList.recordsList.get(row).isTreeSet()) {
					label = Utils.CHOOSE_FILE;
				} else {
					label = dataList.recordsList.get(row).getName();
				}

				break;

			case TreesTableModel.TAXA_SET_INDEX:

				if (!dataList.recordsList.get(row).isTaxaSet()) {
					label = Utils.EDIT_TAXA_SET;
				} else {
					label = dataList.recordsList.get(row).getName();
				}

				break;

			default:

				break;
			}// END: switch

			button.setText(label);
		}// END: null check

		// tooltips
		if (toolTipText != null) {
			button.setToolTipText(toolTipText);
		}
		
		return button;
	}// END: getTableCellRendererComponent

	public void setColumnToolTipText(String toolTipText) {
		this.toolTipText = toolTipText;
	}// END: setColumnToolTipText
	
//	public void setDataList(PartitionDataList dataList) {
//		this.dataList = dataList;
//	}//END: setDataList
	
}// END: class
