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
import javax.swing.JPanel;
import javax.swing.JTextField;
import javax.swing.SwingConstants;
import javax.swing.SwingWorker;

import dr.evolution.io.Importer.ImportException;
import dr.evolution.io.NewickImporter;
import dr.evolution.io.NexusImporter;
import dr.evolution.tree.Tree;
import dr.evomodel.tree.TreeModel;

@SuppressWarnings("serial")
public class TreePanel extends JPanel implements Exportable {

	private BeagleSequenceSimulatorFrame frame = null;
	private PartitionDataList dataList = null;
	private OptionsPanel optionPanel;

	private JButton treeFileButton = new JButton("Choose File...");
	private JTextField treeFileNameText = new JTextField("not selected", 16);

	public TreePanel(final BeagleSequenceSimulatorFrame frame,
			final PartitionDataList dataList) {

		super();

		this.frame = frame;
		this.dataList = dataList;

		setOpaque(false);
		setLayout(new BorderLayout());
		optionPanel = new OptionsPanel(12, 12, SwingConstants.CENTER);
		add(optionPanel, BorderLayout.NORTH);

		treeFileNameText.setEditable(false);
		treeFileButton.addActionListener(new ListenTreeFileButton());

		JPanel tmpPanel = new JPanel(new BorderLayout(0, 0));
		tmpPanel.setOpaque(false);
		tmpPanel.add(treeFileNameText, BorderLayout.CENTER);
		tmpPanel.add(treeFileButton, BorderLayout.EAST);
		optionPanel.addComponentWithLabel("Input Tree File: ", tmpPanel);

	}// END: Constructor

	private class ListenTreeFileButton implements ActionListener {
		public void actionPerformed(ActionEvent ae) {

			doImport();

		}// END: actionPerformed
	}// END: ListenTreeFileButton

	public void doImport() {

		try {

			JFileChooser chooser = new JFileChooser();
			chooser.setDialogTitle("Select input tree file...");
			chooser.setMultiSelectionEnabled(false);
			chooser.setCurrentDirectory(frame.getWorkingDirectory());

			int returnValue = chooser.showOpenDialog(Utils.getActiveFrame());

			if (returnValue == JFileChooser.APPROVE_OPTION) {

				File file = chooser.getSelectedFile();

				if (file != null) {

//					dataList.treeFilesList.add(file);
					treeFileNameText.setText(file.getName());

					importFromFile(file);

					File tmpDir = chooser.getCurrentDirectory();
					if (tmpDir != null) {
						frame.setWorkingDirectory(tmpDir);
					}

				}// END: file opened check
			}// END: dialog cancelled check

		} catch (IOException e) {
			e.printStackTrace();
		} catch (ImportException e) {
			e.printStackTrace();
		}// END: try-catch block

	}// END: doImport

	public void importFromFile(final File file) throws IOException,
			ImportException {

		frame.setBusy();
		SwingWorker<Void, Void> worker = new SwingWorker<Void, Void>() {

			// Executed in background thread
			public Void doInBackground() {

				try {

					BufferedReader reader = new BufferedReader(new FileReader(
							file));

					String line = reader.readLine();
					Tree tree;

					if (line.toUpperCase().startsWith("#NEXUS")) {
						NexusImporter importer = new NexusImporter(reader);
						tree = importer.importTree(null);
					} else {
						NewickImporter importer = new NewickImporter(reader);
						tree = importer.importTree(null);
					}

					// TODO Add taxons from a new tree to that list and display
					// them in Taxa panel
					dataList.taxonList = tree;
					dataList.forestMap.put(file, new TreeModel(tree));
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

	}// END: importFromFile

	public JComponent getExportableComponent() {
		return this;
	}// END: getExportableComponent

}// END: class