package dr.app.bss;

import jam.framework.DocumentFrame;
import jam.framework.Exportable;

import java.awt.BorderLayout;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

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

import dr.app.beagle.tools.BeagleSequenceSimulator;
import dr.evolution.io.Importer.ImportException;

@SuppressWarnings("serial")
public class BeagleSequenceSimulatorFrame extends DocumentFrame {

	private JTabbedPane tabbedPane = new JTabbedPane();
	private TaxaPanel taxaPanel;
	private TreePanel treePanel;
	private BranchSubstitutionModelPanel substModelPanel;
	private ClockRateModelPanel clockPanel;
	private FrequencyModelPanel frequencyPanel;
	private SiteRateModelPanel sitePanel;
	private SimulationPanel simulationPanels;
	
	private BeagleSequenceSimulatorData data = null;
	private JLabel statusLabel;
	private File workingDirectory = null;
	
	public BeagleSequenceSimulatorFrame(String title) {

		super();
		
		setTitle(title);
		data = new BeagleSequenceSimulatorData();

	}// END: Constructor

	@Override
	protected void initializeComponents() {

		setSize(new Dimension(800, 600));
		setMinimumSize(new Dimension(260, 100));

		taxaPanel = new TaxaPanel(this, data);
		treePanel = new TreePanel(this, data);
		substModelPanel = new BranchSubstitutionModelPanel(this, data);
		clockPanel = new ClockRateModelPanel(this, data);
		frequencyPanel = new FrequencyModelPanel(this, data);
		sitePanel = new SiteRateModelPanel(this, data);
		simulationPanels = new SimulationPanel(this, data);
		
		tabbedPane.addTab("Taxa", null, taxaPanel);
		tabbedPane.addTab("Tree", null, treePanel);
		tabbedPane.addTab("Branch Substitution Model", null, substModelPanel);
		tabbedPane.addTab("Clock Rate Model", null, clockPanel);
		tabbedPane.addTab("Frequency Model", null, frequencyPanel);
		tabbedPane.addTab("Site Rate Model", null, sitePanel);
		tabbedPane.addTab("Simulation", null, simulationPanels);
		
		statusLabel = new JLabel("No taxa loaded");

		JPanel progressPanel = new JPanel(new BorderLayout(0, 0));
		JProgressBar progressBar = new JProgressBar();
		progressPanel.add(progressBar, BorderLayout.CENTER);

		JPanel statusPanel = new JPanel(new BorderLayout(0, 0));
		statusPanel.add(statusLabel, BorderLayout.CENTER);
		statusPanel.add(progressPanel, BorderLayout.EAST);
		statusPanel.setBorder(new BorderUIResource.EmptyBorderUIResource(
				new Insets(0, 6, 0, 6)));

		JPanel tabbedPanePanel = new JPanel(new BorderLayout(0, 0));
		tabbedPanePanel.add(tabbedPane, BorderLayout.CENTER);
		tabbedPanePanel.add(statusPanel, BorderLayout.SOUTH);
		tabbedPanePanel.setBorder(new BorderUIResource.EmptyBorderUIResource(
				new Insets(12, 12, 12, 12)));

		getContentPane().setLayout(new java.awt.BorderLayout(0, 0));
		getContentPane().add(tabbedPanePanel, BorderLayout.CENTER);

	}// END: initializeComponents

	// ///////////////////
	// ---IMPORT TAXA---//
	// ///////////////////

//	public Action getImportAction() {
//		return importTaxaAction;
//	}
//
//	private AbstractAction importTaxaAction = new AbstractAction("Import Taxa...") {
//		public void actionPerformed(ActionEvent ae) {
//			 doImport();
//		}
//	};
//
//	public void doImport() {
//
//		try {
//
//			JFileChooser chooser = new JFileChooser();
//			chooser.setDialogTitle("Import Taxa...");
//			chooser.setMultiSelectionEnabled(false);
//			chooser.setCurrentDirectory(workingDirectory);
//
//			chooser.showOpenDialog(Utils.getActiveFrame());
//			File file = chooser.getSelectedFile();
//
//			if (file != null) {
//
//				importFromFile(file);
//
//				File tmpDir = chooser.getCurrentDirectory();
//				if (tmpDir != null) {
//					workingDirectory = tmpDir;
//				}
//
//			}// END: file opened check
//
//		} catch (IOException e) {
//			e.printStackTrace();
//		} catch (ImportException e) {
//			e.printStackTrace();
//		}// END: try catch block
//
//	}// END: doImport
//	
//	public void importFromFile(File file) throws IOException, ImportException {
//		
//        BufferedReader reader = new BufferedReader(new FileReader(file));
//
//        String line = reader.readLine();
//        Tree tree;
//
//        if (line.toUpperCase().startsWith("#NEXUS")) {
//            NexusImporter importer = new NexusImporter(reader);
//            tree = importer.importTree(null);
//        } else {
//            NewickImporter importer = new NewickImporter(reader);
//            tree = importer.importTree(null);
//        }
//
//        data.taxonList = tree;
//        reader.close();
//
//        fireTaxaChanged();
//		
//	}//END: importFromFile

	public void fireTaxaChanged() {
		
		SwingUtilities.invokeLater(new Runnable() {
		    public void run() {

		        taxaPanel.updateUI();
		        statusLabel.setText(Integer.toString(data.taxonList.getTaxonCount()) + " taxa loaded.");
		    
		    }
		});

	}// END: fireTaxaChanged

	// ////////////////
	// ---SIMULATE---//
	// ////////////////

	public Action getExportAction() {
		return simulateAction;
	}

	private AbstractAction simulateAction = new AbstractAction("Simulate...") {
		public void actionPerformed(ActionEvent ae) {
			doExport();
		}
	};
	
	public final void doExport() {

		try {

			JFileChooser chooser = new JFileChooser();
			chooser.setDialogTitle("Simulate...");
			chooser.setMultiSelectionEnabled(false);
			chooser.setCurrentDirectory(workingDirectory);

			chooser.showOpenDialog(Utils.getActiveFrame());
			File file = chooser.getSelectedFile();

			if (file != null) {

				collectAllSettings();
				generateFile(file);

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

	}// END: doExport

	// TODO: simulate
	private void generateFile(File outFile) throws IOException, ImportException {

		PrintWriter writer = new PrintWriter(new FileWriter(outFile));

		BeagleSequenceSimulator beagleSequenceSimulator = new BeagleSequenceSimulator(
				data.treeModel, //
				data.createBranchSubstitutionModel(), //
				data.createSiteRateModel(), //
				data.createBranchRateModel(), //
				data.createFrequencyModel(), //
				data.replicateCount //
		);

		writer.println(beagleSequenceSimulator.simulate().toString());
		writer.close();

	}// END: generateFile
	    
	@Override
	public JComponent getExportableComponent() {
		JComponent exportable = null;
		Component component = tabbedPane.getSelectedComponent();

		if (component instanceof Exportable) {
			exportable = ((Exportable) component).getExportableComponent();
		} else if (component instanceof JComponent) {
			exportable = (JComponent) component;
		}

		return exportable;
	}// END: getExportableComponent

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

	private void collectAllSettings() {
		
		// TODO: collect all settings from panels
		frequencyPanel.collectSettings();
		simulationPanels.collectSettings();
		substModelPanel.collectSettings();
		
	}// END: collectAllSettings
	
	public void dataSelectionChanged(boolean isSelected) {
		if (isSelected) {
			getDeleteAction().setEnabled(true);
		} else {
			getDeleteAction().setEnabled(false);
		}
	}// END: dataSelectionChanged

	public File getWorkingDirectory() {
		return workingDirectory;
	}// END: getWorkingDirectory

	public void setWorkingDirectory(File workingDirectory) {
		this.workingDirectory = workingDirectory;
	}// END: setWorkingDirectory
	
	public void fireModelChanged() {
		substModelPanel.collectSettings();
	}
	
}// END: class
