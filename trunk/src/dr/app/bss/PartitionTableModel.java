package dr.app.bss;

import java.util.ArrayList;
import javax.swing.table.AbstractTableModel;
import dr.evomodel.tree.TreeModel;

@SuppressWarnings("serial")
public class PartitionTableModel extends AbstractTableModel {

	public static final int PARTITION_TREE_INDEX = 0;
	public static final int BRANCH_SUBSTITUTION_MODEL_INDEX = 1;
	public static final int SITE_RATE_MODEL_INDEX = 2;
	public static final int CLOCK_RATE_MODEL_INDEX = 3;
	public static final int FREQUENCY_MODEL_INDEX = 4;

	private String[] columnNames = { "Partition Tree",
			"Branch Substitution Model", "Site Rate Model", "Clock Rate Model",
			"Frequency Model" };

	private ArrayList<BeagleSequenceSimulatorData> dataList;

	public PartitionTableModel(ArrayList<BeagleSequenceSimulatorData> dataList) {
		this.dataList = dataList;
	}// END: constructor

	@Override
	public int getColumnCount() {
		return columnNames.length;
	}//END: getColumnCount

	@Override
	public int getRowCount() {
		return dataList.size();
	}//END: getRowCount

	@Override
	public Object getValueAt(int row, int column) {
		switch (column) {
		case PARTITION_TREE_INDEX:
			return dataList.get(row).treeModel;
		case BRANCH_SUBSTITUTION_MODEL_INDEX:
			return dataList.get(row).substitutionModel;
		case SITE_RATE_MODEL_INDEX:
			return dataList.get(row).siteModel;
		case CLOCK_RATE_MODEL_INDEX:
			return dataList.get(row).clockModel;
		case FREQUENCY_MODEL_INDEX:
			return dataList.get(row).frequencyModel;
		default:
			return new Object();
		}
	}// END: getValueAt
	
	public void setValueAt(Object value, int row, int column) {
		
		switch (column) {
		case PARTITION_TREE_INDEX:
			 dataList.get(row).treeModel = (TreeModel) value;
		case BRANCH_SUBSTITUTION_MODEL_INDEX:
			 dataList.get(row).substitutionModel = (Integer) value;
		case SITE_RATE_MODEL_INDEX:
			 dataList.get(row).siteModel= (Integer) value;
		case CLOCK_RATE_MODEL_INDEX:
			 dataList.get(row).clockModel= (Integer) value;
		case FREQUENCY_MODEL_INDEX:
			 dataList.get(row).frequencyModel= (Integer) value;
		default:
			System.out.println("invalid index");
		}

		fireTableCellUpdated(row, column);
	}

	public String getColumnName(int column) {
		return columnNames[column];
	}//END: getColumnName

	public String[] getColumn(int index) {

		String[] column = new String[dataList.size()];

		for (int i = 0; i < dataList.size(); i++) {
			column[i] = String.valueOf(getValueAt(i, index));
		}

		return column;
	}// END: getColumn
	
	public boolean isCellEditable(int row, int column) {
			return true;
	}

}// END: class
