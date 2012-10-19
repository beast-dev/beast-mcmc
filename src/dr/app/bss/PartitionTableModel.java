package dr.app.bss;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.ArrayList;

import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JOptionPane;
import javax.swing.table.AbstractTableModel;

@SuppressWarnings("serial")
public class PartitionTableModel extends AbstractTableModel {

	public final int PARTITION_TREE_INDEX = 0;
	public final int FROM_INDEX = 1;
	public final int TO_INDEX = 2;
	public final int EVERY_INDEX = 3;
	public final int BRANCH_SUBSTITUTION_MODEL_INDEX = 4;
	public final int SITE_RATE_MODEL_INDEX = 5;
	public final int CLOCK_RATE_MODEL_INDEX = 6;
	public final int FREQUENCY_MODEL_INDEX = 7;

	private String[] COLUMN_NAMES = { "Partition Tree", "From", "To", "Every",
			"Branch Substitution Model", "Site Rate Model", "Clock Rate Model",
			"Frequency Model" };

	private static final Class<?>[] COLUMN_TYPES = new Class<?>[] {
			String.class, Integer.class, Integer.class, Integer.class,
			JComboBox.class, JButton.class, JButton.class, JButton.class };
	
	private ArrayList<PartitionData> dataList;

	public PartitionTableModel(ArrayList<PartitionData> dataList) {
		this.dataList = dataList;
	}// END: constructor

	@Override
	public int getColumnCount() {
		return COLUMN_NAMES.length;
	}//END: getColumnCount

	@Override
	public int getRowCount() {
		return dataList.size();
	}//END: getRowCount

	@Override
	public Class<?> getColumnClass(int columnIndex) {
		return COLUMN_TYPES[columnIndex];
	}//END: getColumnClass
	
	@Override
	public Object getValueAt(final int row, final int column) {
		switch (column) {
		case PARTITION_TREE_INDEX:
			return "TODO";
		case FROM_INDEX:
			return dataList.get(row).from;
		case TO_INDEX:
			return dataList.get(row).to;
		case EVERY_INDEX:
			return dataList.get(row).every;
		case BRANCH_SUBSTITUTION_MODEL_INDEX: // fall through
		case SITE_RATE_MODEL_INDEX:
		case CLOCK_RATE_MODEL_INDEX:
		case FREQUENCY_MODEL_INDEX:
			final JButton button = new JButton(COLUMN_NAMES[column]);
			button.addActionListener(new ActionListener() {
				public void actionPerformed(ActionEvent arg0) {
					JOptionPane.showMessageDialog(
							JOptionPane.getFrameForComponent(button),
							"Button clicked for row " + row + " column "
									+ column);
				}
			});
			return button;
		default:
			return "Error";
		}

	}// END: getValueAt
	
	public void setValueAt(Object value, int row, int column) {
		
		switch (column) {
		case PARTITION_TREE_INDEX:
//			 dataList.get(row).treeModel = (TreeModel) value;
			System.out.println("FUBAR1");
		case FROM_INDEX:
			System.out.println("FUBAR2");
		case TO_INDEX:
			System.out.println("FUBAR3");
		case EVERY_INDEX:
			System.out.println("FUBAR4");
		case BRANCH_SUBSTITUTION_MODEL_INDEX:
//			 dataList.get(row).substitutionModel = (Integer) value;
		case SITE_RATE_MODEL_INDEX:
//			 dataList.get(row).siteModel= (Integer) value;
		case CLOCK_RATE_MODEL_INDEX:
//			 dataList.get(row).clockModel= (Integer) value;
		case FREQUENCY_MODEL_INDEX:
//			 dataList.get(row).frequencyModel= (Integer) value;
		default:
			System.out.println("invalid index");
		}

		fireTableCellUpdated(row, column);
	}

	public String getColumnName(int column) {
		return COLUMN_NAMES[column];
	}//END: getColumnName

	public String[] getColumn(int index) {

		String[] column = new String[dataList.size()];

		for (int i = 0; i < dataList.size(); i++) {
			column[i] = String.valueOf(getValueAt(i, index));
		}

		return column;
	}// END: getColumn
	
	public boolean isCellEditable(int row, int column) {
		switch (column) {
		case PARTITION_TREE_INDEX:
			return false;
		case FROM_INDEX:
			return true;
		case TO_INDEX:
			return true;
		case EVERY_INDEX:
			return true;
		case BRANCH_SUBSTITUTION_MODEL_INDEX:
			return false;
		case SITE_RATE_MODEL_INDEX:
			return false;
		case CLOCK_RATE_MODEL_INDEX:
			return false;
		case FREQUENCY_MODEL_INDEX:
			return false;
		default:
			return false;
		}
	}// END: isCellEditable

}// END: class
