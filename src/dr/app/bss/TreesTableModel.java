package dr.app.bss;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;

import javax.swing.JButton;
import javax.swing.table.AbstractTableModel;

@SuppressWarnings("serial")
public class TreesTableModel extends AbstractTableModel {

	private PartitionDataList dataList;

	public static String[] COLUMN_NAMES = { "Tree File", "Taxa", "Trees" };
	private static final Class<?>[] COLUMN_TYPES = new Class<?>[] {
			JButton.class, Integer.class, Integer.class };
	
	public final static int TREE_FILE_INDEX = 0;
	public final static int TAXA_INDEX = 1;
	public final static int TREES_INDEX = 2;
	
	
	
	public TreesTableModel(PartitionDataList dataList) {
		this.dataList = dataList;
	}// END: Constructor
	
	public void addRow(File treeFile) {
		dataList.treeFileList.add(treeFile);
		fireTableDataChanged();
	}
	
	public void deleteRow(int row) {
		dataList.treeFileList.remove(row);
		fireTableDataChanged();
	}
	
	
	@Override
	public Class<?> getColumnClass(int columnIndex) {
		return COLUMN_TYPES[columnIndex];
	}// END: getColumnClass
	
	public boolean isCellEditable(int row, int column) {
		switch (column) {
		case TREE_FILE_INDEX:
			return false;
		case TAXA_INDEX:
			return false;
		case TREES_INDEX:
			return false;
		default:
			return false;
		}
	}// END: isCellEditable
	
	public String getColumnName(int column) {
		return COLUMN_NAMES[column];
	}// END: getColumnName
	
	@Override
	public int getColumnCount() {
		return COLUMN_NAMES.length;
	}

	@Override
	public int getRowCount() {
		//TODO
		return 1;//dataList.treeFileList.size();
	}

	@Override
	public Object getValueAt(int row, int column) {
		switch (column) {
		case TREE_FILE_INDEX:
			
			 JButton treeFileButton = new JButton(
					COLUMN_NAMES[column]);
			 treeFileButton
					.addActionListener(new ListenLoadTreeFile(
							row));
			return treeFileButton;
			
		case TAXA_INDEX:
			return 0;
			
		case TREES_INDEX:
			return 0;
			
		default:
			return "Error";
		}
	}// END: getValueAt


	//TODO
	private class ListenLoadTreeFile implements ActionListener {

		private int row;

		public ListenLoadTreeFile(int row) {
			this.row = row;
		}// END: Constructor

		public void actionPerformed(ActionEvent ev) {

			System.out.println("TODO");
			
			// branchSubstitutionModelEditor = new
			// BranchSubstitutionModelEditor(
			// dataList, row);
			// branchSubstitutionModelEditor.launch();

		}// END: actionPerformed
	}// END: ListenOpenBranchSubstitutionModelEditor
	
	public void setDataList(PartitionDataList dataList) {
		this.dataList = dataList;
	}
	
}//END: class
