package dr.app.bss;

import jam.framework.DocumentFrame;
import jam.framework.Exportable;

import java.awt.BorderLayout;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.event.ActionEvent;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import javax.swing.AbstractAction;
import javax.swing.Action;
import javax.swing.JComponent;
import javax.swing.JFileChooser;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JProgressBar;
import javax.swing.JTabbedPane;
import javax.swing.SwingUtilities;
import javax.swing.plaf.BorderUIResource;

import dr.evolution.io.Importer.ImportException;
import dr.evolution.io.NewickImporter;
import dr.evolution.io.NexusImporter;
import dr.evolution.tree.Tree;

@SuppressWarnings("serial")
public class BeagleSequenceSimulatorFrame extends DocumentFrame {

	private JTabbedPane tabbedPane = new JTabbedPane();
	private TaxaPanel taxaPanel;
	private BeagleSequenceSimulatorData data = null;
	private JLabel statusLabel;
	private File workingDirectory = null;
	
	// TODO:
	// taxa pane
	// tree
	// clock model
	// frequency model
	// branch subst model
	// site model (rate cats)

	public BeagleSequenceSimulatorFrame(String title) {

		super();
		setTitle(title);
		data = new BeagleSequenceSimulatorData();

		taxaPanel = new TaxaPanel(this, data);
		JPanel secondPanel = new JPanel();

		tabbedPane.addTab("Taxa", null, taxaPanel);
		tabbedPane.addTab("second", null, secondPanel);

		statusLabel = new JLabel("No taxa loaded");

		JPanel progressPanel = new JPanel(new BorderLayout(0, 0));
		JLabel progressLabel = new JLabel("");
		JProgressBar progressBar = new JProgressBar();
		progressPanel.add(progressLabel, BorderLayout.NORTH);
		progressPanel.add(progressBar, BorderLayout.CENTER);

		JPanel panel1 = new JPanel(new BorderLayout(0, 0));
		panel1.add(statusLabel, BorderLayout.CENTER);
		panel1.add(progressPanel, BorderLayout.EAST);
		panel1.setBorder(new BorderUIResource.EmptyBorderUIResource(
				new java.awt.Insets(0, 6, 0, 6)));

		JPanel panel = new JPanel(new BorderLayout(0, 0));
		panel.add(tabbedPane, BorderLayout.CENTER);
		panel.add(panel1, BorderLayout.SOUTH);
		panel.setBorder(new BorderUIResource.EmptyBorderUIResource(
				new java.awt.Insets(12, 12, 12, 12)));

		getContentPane().setLayout(new java.awt.BorderLayout(0, 0));
		getContentPane().add(panel, BorderLayout.CENTER);

	}// END: Constructor

	// ///////////////////
	// ---IMPORT TAXA---//
	// ///////////////////

	public boolean useImportAction() {
		return true;
	}

	public Action getImportAction() {
		return importTaxaAction;
	}

	private AbstractAction importTaxaAction = new AbstractAction("Import Taxa...") {
		public void actionPerformed(ActionEvent ae) {
			 doImport();
		}
	};

	public final void doImport() {

		try {

			JFileChooser chooser = new JFileChooser();
			chooser.setDialogTitle("Import Tree or Alignment...");
			chooser.setMultiSelectionEnabled(false);
			chooser.setCurrentDirectory(workingDirectory);

			chooser.showOpenDialog(Utils.getActiveFrame());
			File file = chooser.getSelectedFile();

			if (file != null) {

				importFromFile(file);

				File tmpDir = chooser.getCurrentDirectory();
				if (tmpDir != null) {
					workingDirectory = tmpDir;
				}

			}// END: file opened check

		} catch (IOException e) {
			e.printStackTrace();
		} catch (ImportException e) {
			e.printStackTrace();
		}// END: try catch block

	}// END: doImport
	
	private void importFromFile(File file) throws IOException, ImportException {
		
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
        statusLabel.setText(Integer.toString(data.taxonList.getTaxonCount()) + " taxa loaded.");
        reader.close();

        fireTaxaChanged();
		
	}//END: importFromFile

	// /////////////////////////
	// ---INHERITED METHODS---//
	// /////////////////////////

	@Override
	public JComponent getExportableComponent() {
		JComponent exportable = null;
		Component comp = tabbedPane.getSelectedComponent();

		if (comp instanceof Exportable) {
			exportable = ((Exportable) comp).getExportableComponent();
		} else if (comp instanceof JComponent) {
			exportable = (JComponent) comp;
		}

		return exportable;
	}// END: getExportableComponent

	@Override
	protected void initializeComponents() {

		setSize(new Dimension(800, 600));
		setMinimumSize(new Dimension(260, 100));

		taxaPanel = new TaxaPanel(this, data);

	}// END: initializeComponents

	@Override
	protected boolean readFromFile(File arg0) throws IOException {
		// TODO Auto-generated method stub
		return false;
	}

	@Override
	protected boolean writeToFile(File arg0) throws IOException {
		// TODO Auto-generated method stub
		return false;
	}

	//TODO: does not work
	public void fireTaxaChanged() {
		
		SwingUtilities.invokeLater(new Runnable() {
		    public void run() {

		    tabbedPane.repaint();
		    
		    }
		});

	}// END: fireTaxaChanged

	public void dataSelectionChanged(boolean isSelected) {
		if (isSelected) {
			getDeleteAction().setEnabled(true);
		} else {
			getDeleteAction().setEnabled(false);
		}
	}// END: dataSelectionChanged

	// /////////////////
	// ---DEBUGGING---//
	// /////////////////
	
}// END: class
