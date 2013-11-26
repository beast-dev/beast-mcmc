package dr.app.bss;

import java.awt.Component;

import javax.swing.DefaultCellEditor;
import javax.swing.JComboBox;
import javax.swing.JTable;


@SuppressWarnings("serial")
public class JTableComboBoxCellEditor extends DefaultCellEditor {

	private PartitionDataList dataList = null;
	
	public JTableComboBoxCellEditor(PartitionDataList dataList) {
		super(new JComboBox());
		this.dataList = dataList;
	}// END: Constructor
	
	@SuppressWarnings("unchecked")
	public Component getTableCellEditorComponent(JTable table,
			Object value, boolean isSelected, int row, int column) {

		((JComboBox)editorComponent).removeAllItems();

		if (column == PartitionTableModel.DATA_INDEX) {

			for (TreesTableRecord record : dataList.recordsList) {
				((JComboBox)editorComponent).addItem(record);
			}// END: fill loop
			
		} else if (column == PartitionTableModel.DATA_TYPE_INDEX) {

			for (String dataType : PartitionData.dataTypes) {
				((JComboBox)editorComponent).addItem(dataType);
			}// END: fill loop

		} else {

			// do nothing

		}// END: column check

		((JComboBox)editorComponent).setSelectedItem(value);
		delegate.setValue(value);

		return editorComponent;
	}// END: getTableCellEditorComponent

	public void setDataList(PartitionDataList dataList) {
		this.dataList = dataList;
	}// END: setDataList
	
}// END: JTableComboBoxCellEditor class