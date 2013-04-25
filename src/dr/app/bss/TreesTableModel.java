package dr.app.bss;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;

import javax.swing.JButton;
import javax.swing.JFileChooser;
import javax.swing.SwingWorker;
import javax.swing.table.AbstractTableModel;

import dr.evolution.util.Taxon;
import dr.evomodel.tree.TreeModel;

@SuppressWarnings("serial")
public class TreesTableModel extends AbstractTableModel {

	private PartitionDataList dataList;
	private MainFrame frame;

	public static String[] COLUMN_NAMES = { "Tree File", "Taxa"};
	private static final Class<?>[] COLUMN_TYPES = new Class<?>[] {
			JButton.class, Integer.class };

	public final static int TREE_FILE_INDEX = 0;
	public final static int TAXA_INDEX = 1;

	public TreesTableModel(PartitionDataList dataList, MainFrame frame) {
		this.dataList = dataList;
		this.frame = frame;
		
		addDefaultRow();
		
	}// END: Constructor

	public void addDefaultRow() {
		dataList.treeFileList.add(new File(""));
		dataList.taxaCounts.add(new Integer(0));
		fireTableDataChanged();
	}//END: addDefaultRow

	public void deleteRow(int row) {
		// remove taxa connected to this row 
		String value = dataList.treeFileList.get(row).getName();
		removeTaxaWithAttributeValue(dataList, Utils.TREE_FILENAME, value);
		dataList.treeFileList.remove(row);
		dataList.taxaCounts.remove(row);
		fireTableDataChanged();
	}//END: deleteRow

	public void setRow(int row, File file, Integer taxaCount) {
		
		// remove taxa connected to this row 
		String value = dataList.treeFileList.get(row).getName();
		removeTaxaWithAttributeValue(dataList, Utils.TREE_FILENAME, value);
		
		dataList.treeFileList.remove(row);
		dataList.treeFileList.add(row, file);
		
		dataList.taxaCounts.remove(row);
		dataList.taxaCounts.add(row, taxaCount);
		
		fireTableDataChanged();
		
	}//END: setRow
	
	@Override
	public int getColumnCount() {
		return COLUMN_NAMES.length;
	}

	@Override
	public int getRowCount() {
		return dataList.treeFileList.size();
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
		default:
			return false;
		}
	}// END: isCellEditable

	public String getColumnName(int column) {
		return COLUMN_NAMES[column];
	}// END: getColumnName

//	//TODO: like in taxa panel
//	private void getTaxaCount() {
//		
//	}
	
	@Override
	public Object getValueAt(int row, int column) {
		switch (column) {
		case TREE_FILE_INDEX:

			JButton treeFileButton = new JButton(Utils.CHOOSE_FILE);
			treeFileButton.addActionListener(new ListenLoadTreeFile(row));
			return treeFileButton;

		case TAXA_INDEX:
			return dataList.taxaCounts.get(row);
			
		default:
			return "Error";
		}
	}// END: getValueAt

	private class ListenLoadTreeFile implements ActionListener {

		private int row;

		public ListenLoadTreeFile(int row) {
			this.row = row;
		}// END: Constructor

		public void actionPerformed(ActionEvent ev) {

			doLoadTreeFile(row);

		}// END: actionPerformed
	}// END: ListenLoadTreeFile

	private void doLoadTreeFile(int row) {

		try {

			JFileChooser chooser = new JFileChooser();
			chooser.setDialogTitle("Select trees file...");
			chooser.setMultiSelectionEnabled(false);
			chooser.setCurrentDirectory(frame.getWorkingDirectory());

			int returnValue = chooser.showOpenDialog(Utils.getActiveFrame());

			if (returnValue == JFileChooser.APPROVE_OPTION) {

				File file = chooser.getSelectedFile();

				if (file != null) {

					// check if file unique, throw exception if not
					if (!isFileInList(file)) {

						loadTreeFile(file, row);

					} else {
						Utils.showDialog("File with this name already loaded. Choose different file.");
					}

					File tmpDir = chooser.getCurrentDirectory();
					if (tmpDir != null) {
						frame.setWorkingDirectory(tmpDir);
					}

				}// END: file opened check
			}// END: dialog cancelled check

		} catch (Exception e) {
			Utils.handleException(e);
		}// END: try-catch block

	}// END: doLoadTreeFile

	private boolean isFileInList(File file) {
		boolean exists = false;

		for (File file2 : dataList.treeFileList) {

			if (file.getName()
					.equalsIgnoreCase(file2.getName())) {
				exists = true;
				break;
			}

		}

		return exists;
	}// END: isFileInList

	private void loadTreeFile(final File file, final int row) {

		frame.setBusy();
		SwingWorker<Void, Void> worker = new SwingWorker<Void, Void>() {

			// Executed in background thread
			public Void doInBackground() {

				try {

					TreeModel tree = Utils.importTreeFromFile(file);
					int taxaCount = 0;
					for (Taxon taxon : tree.asList()) {

						// set attributes to be parsed later in Taxa panel
						if (!Utils.taxonExists(taxon, dataList.taxonList)) {

							double absoluteHeight = Utils
									.getAbsoluteTaxonHeight(taxon, tree);

							taxon.setAttribute(Utils.ABSOLUTE_HEIGHT,
									absoluteHeight);

							taxon.setAttribute(Utils.TREE_FILENAME,
									file.getName());

							dataList.taxonList.addTaxon(taxon);

							taxaCount++;
						}// END: taxon exists check

					}// END: taxon loop

					setRow(row, file, taxaCount);
					
				} catch (Exception e) {
					Utils.handleException(e);
				}// END: try-catch block

				return null;
			}// END: doInBackground()

			// Executed in event dispatch thread
			public void done() {
				frame.setIdle();
				frame.fireTaxaChanged();
			}// END: done
		};

		worker.execute();

	}// END: loadTreeFile

	private void removeTaxaWithAttributeValue(PartitionDataList dataList,
			String attribute, String value) {

		synchronized (dataList.taxonList) {
			for (int i = 0; i < dataList.taxonList.getTaxonCount(); i++) {

				Taxon taxon = dataList.taxonList.getTaxon(i);
				if (taxon.getAttribute(attribute).toString()
						.equalsIgnoreCase(value)) {
					dataList.taxonList.removeTaxon(taxon);
					i--;
				}
			}
		}

	}// END: removeTaxaWithAttributeValue
	
	public void setDataList(PartitionDataList dataList) {
		this.dataList = dataList;
	}

}// END: class
