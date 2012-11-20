package dr.app.bss;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;

import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.table.AbstractTableModel;
import javax.swing.text.BadLocationException;

@SuppressWarnings("serial")
public class PartitionTableModel extends AbstractTableModel {

	private BranchSubstitutionModelEditor branchSubstitutionModelEditor;
	private SiteRateModelEditor siteRateModelEditor;
	private ClockRateModelEditor clockRateModelEditor;
	private FrequencyModelEditor frequencyModelEditor;

	public final static int PARTITION_TREE_INDEX = 0;
	public final static int FROM_INDEX = 1;
	public final static int TO_INDEX = 2;
	public final static int EVERY_INDEX = 3;
	public final static int BRANCH_SUBSTITUTION_MODEL_INDEX = 4;
	public final static int SITE_RATE_MODEL_INDEX = 5;
	public final static int CLOCK_RATE_MODEL_INDEX = 6;
	public final static int FREQUENCY_MODEL_INDEX = 7;

	private String[] COLUMN_NAMES = { "Partition Tree", "From", "To", "Every",
			"Branch Substitution Model", "Site Rate Model", "Clock Rate Model",
			"Frequency Model" };

	private static final Class<?>[] COLUMN_TYPES = new Class<?>[] {
			JComboBox.class, Integer.class, Integer.class, Integer.class,
			JButton.class, JButton.class, JButton.class, JButton.class };

	private PartitionDataList dataList;

	public PartitionTableModel(PartitionDataList dataList) {
		this.dataList = dataList;
	}// END: Constructor

	@Override
	public int getColumnCount() {
		return COLUMN_NAMES.length;
	}// END: getColumnCount

	@Override
	public int getRowCount() {
		return dataList.size();
	}// END: getRowCount

	@Override
	public Class<?> getColumnClass(int columnIndex) {
		return COLUMN_TYPES[columnIndex];
	}// END: getColumnClass

	@Override
	public Object getValueAt(final int row, final int column) {
		switch (column) {
		case PARTITION_TREE_INDEX:
			
			return dataList.get(row).treeFile == null ? new File("") : dataList.get(row).treeFile;

		case FROM_INDEX:
			return dataList.get(row).from;
		case TO_INDEX:
			return dataList.get(row).to;
		case EVERY_INDEX:
			return dataList.get(row).every;
		case BRANCH_SUBSTITUTION_MODEL_INDEX:

			branchSubstitutionModelEditor = new BranchSubstitutionModelEditor(
					dataList, row);
			final JButton branchSubstModelButton = new JButton(
					COLUMN_NAMES[column]);
			branchSubstModelButton
					.addActionListener(new ListenOpenBranchSubstitutionModelEditor());
			return branchSubstModelButton;

		case SITE_RATE_MODEL_INDEX:

			try {

				siteRateModelEditor = new SiteRateModelEditor(dataList, row);
				final JButton siteRateModelButton = new JButton(
						COLUMN_NAMES[column]);
				siteRateModelButton
						.addActionListener(new ListenOpenSiteRateModelEditor());

				return siteRateModelButton;

			} catch (NumberFormatException e) {
				e.printStackTrace();
			} catch (BadLocationException e) {
				e.printStackTrace();
			}

		case CLOCK_RATE_MODEL_INDEX:

			clockRateModelEditor = new ClockRateModelEditor(dataList, row);
			final JButton clockRateModelButton = new JButton(
					COLUMN_NAMES[column]);
			clockRateModelButton
					.addActionListener(new ListenOpenClockRateModelEditor());
			return clockRateModelButton;

		case FREQUENCY_MODEL_INDEX:

			frequencyModelEditor = new FrequencyModelEditor(dataList, row);
			final JButton frequencyModelButton = new JButton(
					COLUMN_NAMES[column]);
			frequencyModelButton
					.addActionListener(new ListenOpenFrequencyModelEditor());
			return frequencyModelButton;

		default:
			return "Error";
		}

	}// END: getValueAt

	public void setValueAt(Object value, int row, int column) {

		switch (column) {
		case PARTITION_TREE_INDEX:
			dataList.get(row).treeFile = (File) value;
			break;
		case FROM_INDEX:
			dataList.get(row).from = (Integer) value;
			break;
		case TO_INDEX:
			dataList.get(row).to = (Integer) value;
			break;
		case EVERY_INDEX:
			dataList.get(row).every = (Integer) value;
			break;
		case BRANCH_SUBSTITUTION_MODEL_INDEX:
			dataList.get(row).substitutionModel = (Integer) value;
		case SITE_RATE_MODEL_INDEX:
			dataList.get(row).siteModel = (Integer) value;
		case CLOCK_RATE_MODEL_INDEX:
			dataList.get(row).clockModel = (Integer) value;
		case FREQUENCY_MODEL_INDEX:
			dataList.get(row).frequencyModel = (Integer) value;
		default:
			System.out.println("invalid index");
		}

		fireTableCellUpdated(row, column);
		
	}// END: setValueAt

	public void addRow(PartitionData row) {
		dataList.add(row);
		this.fireTableDataChanged();
	}

	public void addDefaultRow() {
		//TODO: this could copy the previous line (new constructor that takes all the elements)
		dataList.add(new PartitionData());
		fireTableRowsInserted(dataList.size() - 1, dataList.size() - 1);
	}

	public void deleteRow(int row) {
		dataList.remove(row);
		this.fireTableDataChanged();
	}

	public String getColumnName(int column) {
		return COLUMN_NAMES[column];
	}// END: getColumnName

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
			return true;
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

	private class ListenOpenBranchSubstitutionModelEditor implements
			ActionListener {
		public void actionPerformed(ActionEvent ev) {

			branchSubstitutionModelEditor.launch();

		}// END: actionPerformed
	}// END: ListenOpenBranchSubstitutionModelEditor

	private class ListenOpenSiteRateModelEditor implements ActionListener {
		public void actionPerformed(ActionEvent ev) {

			siteRateModelEditor.launch();

		}// END: actionPerformed
	}// END: ListenOpenSiteRateModelEditor

	private class ListenOpenClockRateModelEditor implements ActionListener {
		public void actionPerformed(ActionEvent ev) {

			clockRateModelEditor.launch();

		}// END: actionPerformed
	}// END: ListenOpenSiteRateModelEditor

	private class ListenOpenFrequencyModelEditor implements ActionListener {
		public void actionPerformed(ActionEvent ev) {

			frequencyModelEditor.launch();

		}// END: actionPerformed
	}// END: ListenOpenSiteRateModelEditor

}// END: class
