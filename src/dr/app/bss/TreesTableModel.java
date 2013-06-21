package dr.app.bss;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;

import javax.swing.JButton;
import javax.swing.JFileChooser;
import javax.swing.SwingWorker;
import javax.swing.table.AbstractTableModel;

import dr.evolution.tree.Tree;
import dr.evolution.util.Taxa;
import dr.evolution.util.Taxon;

/**
 * @author Filip Bielejec
 * @version $Id$
 */
@SuppressWarnings("serial")
public class TreesTableModel extends AbstractTableModel {

	private PartitionDataList dataList;
	private MainFrame frame;

	public final static int TREE_FILE_INDEX = 0;
	public final static int TAXA_SET_INDEX = 1;
	public final static int TAXA_COUNT_INDEX = 2;
	
	public static String[] COLUMN_NAMES = { "Tree File", "Taxa Set",  "Taxa count" };
	
	private static final Class<?>[] COLUMN_TYPES = new Class<?>[] {
			JButton.class, // Tree File
			JButton.class, // Taxa Set
			Integer.class // Taxa
	};

	private TaxaEditor taxaEditor;
	
	public TreesTableModel(PartitionDataList dataList, MainFrame frame) {
		this.dataList = dataList;
		this.frame = frame;

		addDefaultRow();

	}// END: Constructor

	public void addDefaultRow() {
		dataList.recordsList.add(new TreesTableRecord());
		fireTableDataChanged();
	}// END: addDefaultRow

	public void deleteRow(int row) {
		
		// remove taxa connected to this row from Taxa panel
		String value = dataList.recordsList.get(row).getName();
		removeTaxaWithAttributeValue(dataList, Utils.TREE_FILENAME, value);
		
		dataList.recordsList.remove(row);
		fireTableDataChanged();
	}// END: deleteRow

	public void setRow(int row, TreesTableRecord record) {

		// remove taxa connected to this row from Taxa panel
		String value = dataList.recordsList.get(row).getName();
		removeTaxaWithAttributeValue(dataList, Utils.TREE_FILENAME, value);

		dataList.recordsList.remove(row);
		dataList.recordsList.add(row, record);

		fireTableDataChanged();

	}// END: setRow

	@Override
	public int getColumnCount() {
		return COLUMN_NAMES.length;
	}

	@Override
	public int getRowCount() {
		return dataList.recordsList.size();
	}

	@Override
	public Class<?> getColumnClass(int columnIndex) {
		return COLUMN_TYPES[columnIndex];
	}// END: getColumnClass

	public boolean isCellEditable(int row, int column) {
		switch (column) {
		case TREE_FILE_INDEX:
			return false;
		case TAXA_COUNT_INDEX:
			return false;
		case TAXA_SET_INDEX:
			return false;
		default:
			return false;
		}
	}// END: isCellEditable

	public String getColumnName(int column) {
		return COLUMN_NAMES[column];
	}// END: getColumnName

	@Override
	public Object getValueAt(int row, int column) {
		switch (column) {
		case TREE_FILE_INDEX:

			JButton treeFileButton = new JButton(Utils.CHOOSE_FILE);
			treeFileButton.addActionListener(new ListenLoadTreeFile(row));
			return treeFileButton;

		case TAXA_COUNT_INDEX:
			return dataList.recordsList.get(row).getTaxaCount();

		case TAXA_SET_INDEX:
			
			JButton taxaEditorButton = new JButton(COLUMN_NAMES[column]);
			taxaEditorButton.addActionListener(new ListenOpenTaxaEditor(row));
			return taxaEditorButton;
			
		default:
			return "Error";
		}
	}// END: getValueAt

	private class ListenOpenTaxaEditor implements ActionListener {
		private int row;

		public ListenOpenTaxaEditor(int row) {
			this.row = row;
		}

		public void actionPerformed(ActionEvent ev) {

			try {
				
				taxaEditor = new TaxaEditor(frame, dataList, row);
				taxaEditor.launch();
			
			} catch (Exception e) {
				Utils.handleException(e);
			}

		}// END: actionPerformed
	}// END: ListenOpenSiteRateModelEditor
	
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

					//TODO
					// check if file unique, throw exception if not
//					if (!isFileInList(file)) {

						loadTreeFile(file, row);

//					} else {
//						Utils.showDialog("File with this name already loaded. Choose different file.");
//					}

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

	//TODO
//	private boolean isFileInList(File file) {
//		boolean exists = false;
//
//		for (File file2 : dataList.treesList) {
//
//			if (file.getName().equalsIgnoreCase(file2.getName())) {
//				exists = true;
//				break;
//			}
//
//		}
//
//		return exists;
//	}// END: isFileInList

	private void loadTreeFile(final File file, final int row) {

		frame.setBusy();
		SwingWorker<Void, Void> worker = new SwingWorker<Void, Void>() {

			// Executed in background thread
			public Void doInBackground() {

				try {

					Tree tree = Utils.importTreeFromFile(file);
					tree.setId(file.getName());
					Taxa taxa = new Taxa();
					for (Taxon taxon : tree.asList()) {

						// set attributes to be parsed later in Taxa panel
						if (!Utils.taxonExists(taxon, dataList.allTaxa)) {

							double absoluteHeight = Utils
									.getAbsoluteTaxonHeight(taxon, tree);

							taxon.setAttribute(Utils.ABSOLUTE_HEIGHT,
									absoluteHeight);

							// for demographic models?
							// taxon.setAttribute("date", new Date(absoluteHeight,
							// Units.Type.YEARS, true));
							
							taxon.setAttribute(Utils.TREE_FILENAME,
									tree.getId());

							taxa.addTaxon(taxon);
							dataList.allTaxa.addTaxon(taxon);
							
						}// END: taxon exists check

					}// END: taxon loop

					setRow(row, new TreesTableRecord(tree.getId(), tree, taxa));

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

		synchronized (dataList.allTaxa) {
			for (int i = 0; i < dataList.allTaxa.getTaxonCount(); i++) {

				Taxon taxon = dataList.allTaxa.getTaxon(i);
				if (taxon.getAttribute(attribute).toString()
						.equalsIgnoreCase(value)) {
					dataList.allTaxa.removeTaxon(taxon);
					i--;
				}
			}
		}

	}// END: removeTaxaWithAttributeValue

	public void setDataList(PartitionDataList dataList) {
		this.dataList = dataList;
	}

}// END: class
