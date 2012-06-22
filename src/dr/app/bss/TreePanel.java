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

import dr.evolution.io.Importer.ImportException;
import dr.evolution.io.NewickImporter;
import dr.evolution.io.NexusImporter;
import dr.evolution.tree.Tree;
import dr.evomodel.tree.TreeModel;

@SuppressWarnings("serial")
public class TreePanel extends JPanel implements Exportable {

	private BeagleSequenceSimulatorFrame frame = null;
	private BeagleSequenceSimulatorData data = null;
	private OptionsPanel optionPanel;

	private JButton treeFileButton = new JButton("Choose File...");
	private JTextField treeFileNameText = new JTextField("not selected", 16);

	public TreePanel(final BeagleSequenceSimulatorFrame frame,
			final BeagleSequenceSimulatorData data) {

		super();

		this.frame = frame;
		this.data = data;

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

			chooser.showOpenDialog(Utils.getActiveFrame());
			File file = chooser.getSelectedFile();

			if (file != null) {

				data.treeFile = file;
				treeFileNameText.setText(data.treeFile.getName());

				importFromFile(file);

				File tmpDir = chooser.getCurrentDirectory();
				if (tmpDir != null) {
					frame.setWorkingDirectory(tmpDir);
				}

			}// END: file opened check

		} catch (IOException e) {
			e.printStackTrace();
		} catch (ImportException e) {
			e.printStackTrace();
		}// END: try-catch block

	}// END: doImport
	
	public void importFromFile(File file) throws IOException, ImportException {
		
        BufferedReader reader = new BufferedReader(new FileReader(file));

        String line = reader.readLine();
        Tree tree;

        if (line.toUpperCase().startsWith("#NEXUS")) {
            NexusImporter importer = new NexusImporter(reader);
            tree = importer.importTree(null);
        } else {
            NewickImporter importer = new NewickImporter(reader);
            tree = importer.importTree(null);
        }

        data.taxonList = tree;
        data.treeModel = new TreeModel(tree);
        
        reader.close();
        frame.fireTaxaChanged();
		
	}//END: importFromFile
	
	@Override
	public JComponent getExportableComponent() {
		return this;
	}// END: getExportableComponent

}// END: class