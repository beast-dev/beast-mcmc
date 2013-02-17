package dr.app.bss;

import jam.framework.Exportable;
import jam.panels.OptionsPanel;

import java.awt.BorderLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

import javax.swing.JButton;
import javax.swing.JComponent;
import javax.swing.JFileChooser;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTextField;
import javax.swing.SwingConstants;
import javax.swing.SwingWorker;

import dr.evolution.io.Importer.ImportException;
import dr.evolution.io.NewickImporter;
import dr.evolution.io.NexusImporter;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxon;
import dr.evomodel.tree.TreeModel;

@SuppressWarnings("serial")
public class TreePanel extends JPanel implements Exportable {

	private BeagleSequenceSimulatorFrame frame = null;
	private PartitionDataList dataList = null;
	private OptionsPanel optionPanel;

	private JButton treeFileButton = new JButton("Choose File...");
	private JTextField treeFileNameText = new JTextField("not selected", 16);
	private JLabel treeFileLabel = new JLabel("Input Tree File: ");

	private JButton treesFileButton = new JButton("Choose File...");
	private JTextField treesFileNameText = new JTextField("not selected", 16);
	private JLabel treesFileLabel = new JLabel("Input Trees File: ");
	
	public TreePanel(final BeagleSequenceSimulatorFrame frame,
			final PartitionDataList dataList) {

		super();

		this.frame = frame;
		this.dataList = dataList;
		JPanel tmpPanel;
		
		setOpaque(false);
		setLayout(new BorderLayout());
		optionPanel = new OptionsPanel(12, 12, SwingConstants.CENTER);
		add(optionPanel, BorderLayout.NORTH);

		treeFileNameText.setEditable(false);
		treeFileButton.addActionListener(new ListenTreeFileButton());

		tmpPanel = new JPanel(new BorderLayout(0, 0));
		tmpPanel.setOpaque(false);
		tmpPanel.add(treeFileNameText, BorderLayout.CENTER);
		tmpPanel.add(treeFileButton, BorderLayout.EAST);
		optionPanel.addComponents(treeFileLabel, tmpPanel);

		treesFileNameText.setEditable(false);
		treesFileButton.addActionListener(new ListenTreesFileButton());
		
		tmpPanel = new JPanel(new BorderLayout(0, 0));
		tmpPanel.setOpaque(false);
		tmpPanel.add(treesFileNameText, BorderLayout.CENTER);
		tmpPanel.add(treesFileButton, BorderLayout.EAST);
		optionPanel.addComponents(treesFileLabel, tmpPanel);
		
		
	}// END: Constructor

	// ///////////////////
	// ---IMPORT TREE---//
	// ///////////////////

	private class ListenTreeFileButton implements ActionListener {
		public void actionPerformed(ActionEvent ae) {

			doImportTree();

		}// END: actionPerformed
	}// END: ListenTreeFileButton

	public void doImportTree() {

		try {

			JFileChooser chooser = new JFileChooser();
			chooser.setDialogTitle("Select input tree file...");
			chooser.setMultiSelectionEnabled(false);
			chooser.setCurrentDirectory(frame.getWorkingDirectory());

			int returnValue = chooser.showOpenDialog(Utils.getActiveFrame());

			if (returnValue == JFileChooser.APPROVE_OPTION) {

				File file = chooser.getSelectedFile();

				if (file != null) {

					treeFileNameText.setText(file.getName());

					importTreeFromFile(file);

					File tmpDir = chooser.getCurrentDirectory();
					if (tmpDir != null) {
						frame.setWorkingDirectory(tmpDir);
					}

				}// END: file opened check
			}// END: dialog cancelled check

		} catch (Exception e) {
			Utils.handleException(e);
		}// END: try-catch block

	}// END: doImportTree

	// TODO: this should import from all the different formats
	public void importTreeFromFile(final File file) throws IOException,
			ImportException {

		frame.setBusy();
		SwingWorker<Void, Void> worker = new SwingWorker<Void, Void>() {

			// Executed in background thread
			public Void doInBackground() {

				try {

					BufferedReader reader = new BufferedReader(new FileReader(
							file));

					String line = reader.readLine();

					Tree tree = null;

					if (line.toUpperCase().startsWith("#NEXUS")) {

						NexusImporter importer = new NexusImporter(reader);
						tree = importer.importTree(null);

					} else {

						NewickImporter importer = new NewickImporter(reader);

						// tree = importer.importNextTree();
						tree = importer.importTree(null);

					}

					dataList.forestMap.put(file, new TreeModel(tree));
					for (Taxon taxon : tree.asList()) {

						if (!Utils.taxonExists(taxon, dataList.taxonList)) {
							dataList.taxonList.addTaxon(taxon);
						}// END: taxon exists check

					}

					reader.close();

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

	}// END: importTreeFromFile

	// ////////////////////
	// ---IMPORT TREES---//
	// ////////////////////

	private class ListenTreesFileButton implements ActionListener {
		public void actionPerformed(ActionEvent ae) {

			doImportTrees();

		}// END: actionPerformed
	}// END: ListenTreeFileButton

	public void doImportTrees() {

		try {

			JFileChooser chooser = new JFileChooser();
			chooser.setDialogTitle("Select input trees file...");
			chooser.setMultiSelectionEnabled(false);
			chooser.setCurrentDirectory(frame.getWorkingDirectory());

			int returnValue = chooser.showOpenDialog(Utils.getActiveFrame());

			if (returnValue == JFileChooser.APPROVE_OPTION) {

				File file = chooser.getSelectedFile();

				if (file != null) {

					treesFileNameText.setText(file.getName());
					dataList.treesFilename = file;
					
//					importTreesFromFile(file);

					File tmpDir = chooser.getCurrentDirectory();
					if (tmpDir != null) {
						frame.setWorkingDirectory(tmpDir);
					}

				}// END: file opened check
			}// END: dialog cancelled check

		} catch (Exception e) {
			Utils.handleException(e);
		}// END: try-catch block

	}// END: doImport

//	public void importTreesFromFile(final File file) throws IOException,
//			ImportException {
//
//		frame.setBusy();
//		SwingWorker<Void, Void> worker = new SwingWorker<Void, Void>() {
//
//			// Executed in background thread
//			public Void doInBackground() {
//
//				try {
//
//
//				} catch (Exception e) {
//					Utils.handleException(e);
//				}// END: try-catch block
//
//				return null;
//			}// END: doInBackground()
//
//			// Executed in event dispatch thread
//			public void done() {
//				frame.setIdle();
//			}// END: done
//		};
//
//		worker.execute();
//
//	}// END: importTreesFromFile

	public JComponent getExportableComponent() {
		return this;
	}// END: getExportableComponent

	public void enableTreeFileButton() {
		treeFileLabel.setEnabled(true);
		treeFileNameText.setEnabled(true);
		treeFileButton.setEnabled(true);
	}
	
	public void disableTreeFileButton() {
		treeFileLabel.setEnabled(false);
		treeFileNameText.setEnabled(false);
		treeFileButton.setEnabled(false);
	}
	
	public void enableTreesFileButton() {
		treesFileLabel.setEnabled(true);
		treesFileNameText.setEnabled(true);
		treesFileButton.setEnabled(true);
	}
	
	public void disableTreesFileButton() {
		treesFileLabel.setEnabled(false);
		treesFileNameText.setEnabled(false);
		treesFileButton.setEnabled(false);
	}
	
}// END: class