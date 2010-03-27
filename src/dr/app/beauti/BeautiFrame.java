/*
 * BeautiFrame.java
 *
 * (c) 2002-2005 BEAST Development Core Team
 *
 * This package may be distributed under the
 * Lesser Gnu Public Licence (LGPL)
 */
package dr.app.beauti;

import dr.app.beauti.clockModelsPanel.ClockModelsPanel;
import dr.app.beauti.components.ComponentFactory;
import dr.app.beauti.components.SequenceErrorModelComponentFactory;
import dr.app.beauti.components.TipDateSamplingComponentFactory;
import dr.app.beauti.datapanel.DataPanel;
import dr.app.beauti.enumTypes.TreePriorType;
import dr.app.beauti.generator.BeastGenerator;
import dr.app.beauti.mcmcpanel.MCMCPanel;
import dr.app.beauti.operatorspanel.OperatorsPanel;
import dr.app.beauti.options.*;
import dr.app.beauti.priorsPanel.DefaultPriorDialog;
import dr.app.beauti.priorsPanel.PriorsPanel;
import dr.app.beauti.siteModelsPanel.SiteModelsPanel;
import dr.app.beauti.taxonsetspanel.TaxaPanel;
import dr.app.beauti.tipdatepanel.TipDatesPanel;
import dr.app.beauti.traitspanel.TraitsPanel;
import dr.app.beauti.treespanel.OldTreesPanel;
import dr.app.beauti.treespanel.TreesPanel;
import dr.app.beauti.util.BEAUTiImporter;
import dr.app.java16compat.FileNameExtensionFilter;
import dr.app.util.Arguments;
import dr.app.util.Utils;
import dr.evolution.io.Importer.ImportException;
import dr.evolution.io.NexusImporter.MissingBlockException;
import dr.evolution.util.Taxon;
import org.virion.jam.framework.DocumentFrame;
import org.virion.jam.framework.Exportable;
import org.virion.jam.util.IconUtils;

import javax.swing.*;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import javax.swing.plaf.BorderUIResource;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.EnumSet;
import java.util.List;
import java.util.Map;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: BeautiFrame.java,v 1.22 2006/09/09 16:07:06 rambaut Exp $
 */
public class BeautiFrame extends DocumentFrame {

    private static final long serialVersionUID = 2114148696789612509L;

    private final BeautiOptions beautiOptions;
    private final BeastGenerator generator;

    private final JTabbedPane tabbedPane = new JTabbedPane();
    public final JLabel statusLabel = new JLabel("No data loaded");

    private DataPanel dataPanel;
    private TipDatesPanel tipDatesPanel;
    private TraitsPanel traitsPanel;
    private TaxaPanel taxaPanel;
    private SiteModelsPanel siteModelsPanel;
    private ClockModelsPanel clockModelsPanel;
    private OldTreesPanel oldTreesPanel;
//    private SpeciesTreesPanel speciesTreesPanel;
    private TreesPanel treesPanel;
    private PriorsPanel priorsPanel;
    private OperatorsPanel operatorsPanel;
    private MCMCPanel mcmcPanel;

    private BeautiPanel currentPanel;

    private JFileChooser importChooser; // make JFileChooser chooser remember previous path
    private JFileChooser exportChooser; // make JFileChooser chooser remember previous path

    final Icon gearIcon = IconUtils.getIcon(this.getClass(), "images/gear.png");

    public BeautiFrame(String title) {
        super();

        setTitle(title);

        // Prevent the application to close in requestClose()
        // after a user cancel or a failure in beast file generation
        setDefaultCloseOperation(DO_NOTHING_ON_CLOSE);

        getOpenAction().setEnabled(false);
        getSaveAction().setEnabled(false);

        getFindAction().setEnabled(false);

        getZoomWindowAction().setEnabled(false);

        ComponentFactory[] components = {
                SequenceErrorModelComponentFactory.INSTANCE,
                TipDateSamplingComponentFactory.INSTANCE,
//                DiscreteTraitsComponentFactory.INSTANCE
        };

        beautiOptions = new BeautiOptions(components);
        generator = new BeastGenerator(beautiOptions, components);
    }

    public void initializeComponents() {

        dataPanel = new DataPanel(this, getImportAction(), getDeleteAction(), getImportTraitsAction());
        tipDatesPanel = new TipDatesPanel(this);
        traitsPanel = new TraitsPanel(this, getImportTraitsAction());
        taxaPanel = new TaxaPanel(this);
        siteModelsPanel = new SiteModelsPanel(this, getDeleteAction());
        clockModelsPanel = new ClockModelsPanel(this);
        oldTreesPanel = new OldTreesPanel(this);
        treesPanel = new TreesPanel(this, getDeleteAction());
//        speciesTreesPanel = new SpeciesTreesPanel(this);
        priorsPanel = new PriorsPanel(this, false);
        operatorsPanel = new OperatorsPanel(this);
        mcmcPanel = new MCMCPanel(this);

        tabbedPane.addTab("Data Partitions", dataPanel);
        tabbedPane.addTab("Taxon Sets", taxaPanel);
        tabbedPane.addTab("Tip Dates", tipDatesPanel);
        tabbedPane.addTab("Traits", traitsPanel);
        tabbedPane.addTab("Site Models", siteModelsPanel);
        tabbedPane.addTab("Clock Models", clockModelsPanel);
        if (DataPanel.ALLOW_UNLINKED_TREES) {
            tabbedPane.addTab("Trees", treesPanel);
        } else {
            tabbedPane.addTab("Trees", oldTreesPanel);
        }
        tabbedPane.addTab("Priors", priorsPanel);
        tabbedPane.addTab("Operators", operatorsPanel);
        tabbedPane.addTab("MCMC", mcmcPanel);
        currentPanel = (BeautiPanel) tabbedPane.getSelectedComponent();

        tabbedPane.addChangeListener(new ChangeListener() {
            public void stateChanged(ChangeEvent e) {
                BeautiPanel selectedPanel = (BeautiPanel) tabbedPane.getSelectedComponent();
                if (selectedPanel == dataPanel) {
                    dataPanel.selectionChanged();
                } else {
                    getDeleteAction().setEnabled(false);
                }
                currentPanel.getOptions(beautiOptions);
                setAllOptions();
                currentPanel = selectedPanel;
            }
        });

        JPanel panel = new JPanel(new BorderLayout(6, 6));
        panel.setBorder(new BorderUIResource.EmptyBorderUIResource(new java.awt.Insets(12, 12, 12, 12)));
        panel.add(tabbedPane, BorderLayout.CENTER);
        panel.setPreferredSize(new java.awt.Dimension(800, 600));

        getExportAction().setEnabled(false);
        JButton generateButton = new JButton(getExportAction());
        generateButton.putClientProperty("JButton.buttonType", "roundRect");

        JPanel panel2 = new JPanel(new BorderLayout(6, 6));
        panel2.add(statusLabel, BorderLayout.CENTER);
        panel2.add(generateButton, BorderLayout.EAST);

        panel.add(panel2, BorderLayout.SOUTH);

        JScrollPane scrollPane = new JScrollPane(panel);
        scrollPane.setOpaque(false);

        getContentPane().setLayout(new java.awt.BorderLayout(0, 0));
        getContentPane().add(scrollPane, BorderLayout.CENTER);

        setAllOptions();

        setSize(new java.awt.Dimension(1024, 768));

        // make JFileChooser chooser remember previous path
        exportChooser = new JFileChooser(Utils.getCWD());
        exportChooser.setFileFilter(new FileNameExtensionFilter("BEAST XML File", "xml", "beast"));
        exportChooser.setDialogTitle("Generate BEAST XML File...");


        importChooser = new JFileChooser(Utils.getCWD());

        importChooser.setMultiSelectionEnabled(true);
        importChooser.setFileFilter(new FileNameExtensionFilter(
                "NEXUS (*.nex) & BEAST (*.xml) Files", "nex", "nexus", "nx", "xml", "beast", "fa", "fasta", "afa"));
        importChooser.setDialogTitle("Import Aligment...");
    }

    /**
     * set all the options for all panels
     */
    public void setAllOptions() {
        dataPanel.setOptions(beautiOptions);
        tipDatesPanel.setOptions(beautiOptions);
        traitsPanel.setOptions(beautiOptions);
        taxaPanel.setOptions(beautiOptions);
        siteModelsPanel.setOptions(beautiOptions);
        clockModelsPanel.setOptions(beautiOptions);
//        if (beautiOptions.isSpeciesAnalysis()) {
//            speciesTreesPanel.setOptions(beautiOptions);
//        } else
        if (DataPanel.ALLOW_UNLINKED_TREES) {
            treesPanel.setOptions(beautiOptions);
        } else {
            oldTreesPanel.setOptions(beautiOptions);
        }
        priorsPanel.setOptions(beautiOptions);
        operatorsPanel.setOptions(beautiOptions);
        mcmcPanel.setOptions(beautiOptions);

        setStatusMessage();
    }

    /**
     * get all the options for all panels
     */
    private void getAllOptions() {
        dataPanel.getOptions(beautiOptions);
        tipDatesPanel.getOptions(beautiOptions);
        traitsPanel.getOptions(beautiOptions);
        taxaPanel.getOptions(beautiOptions);
        siteModelsPanel.getOptions(beautiOptions);
        clockModelsPanel.getOptions(beautiOptions);
//        if (beautiOptions.isSpeciesAnalysis()) {
//            speciesTreesPanel.getOptions(beautiOptions);
//        } else
        if (DataPanel.ALLOW_UNLINKED_TREES) {
            treesPanel.getOptions(beautiOptions);
        } else {
            oldTreesPanel.getOptions(beautiOptions);
        }
        priorsPanel.getOptions(beautiOptions);
        operatorsPanel.getOptions(beautiOptions);
        mcmcPanel.getOptions(beautiOptions);
    }

    public void doSelectAll() {
        if (currentPanel == dataPanel) {
            dataPanel.selectAll();
        }
    }

    public final void dataSelectionChanged(boolean isSelected) {
        getDeleteAction().setEnabled(isSelected);
    }

    public final void modelSelectionChanged(boolean isSelected) {
        getDeleteAction().setEnabled(isSelected);
    }

    public void doDelete() {
        if (tabbedPane.getSelectedComponent() == dataPanel) {
            dataPanel.removeSelection();
//        } else if (tabbedPane.getSelectedComponent() == modelsPanel) {
//            modelsPanel.removeSelection();
//        } else if (tabbedPane.getSelectedComponent() == treesPanel) {
//        	treesPanel.removeSelection();
        } else {
            throw new RuntimeException("Delete should only be accessable from the Data and Models panels");
        }

        setStatusMessage();
    }

    public boolean requestClose() {
        if (isDirty() && beautiOptions.hasData()) {
            int option = JOptionPane.showConfirmDialog(this,
                    "You have made changes but have not generated\n" +
                            "a BEAST XML file. Do you wish to generate\n" +
                            "before closing this window?",
                    "Unused changes",
                    JOptionPane.YES_NO_CANCEL_OPTION,
                    JOptionPane.WARNING_MESSAGE);

            if (option == JOptionPane.YES_OPTION) {
                return !doGenerate();
            } else if (option == JOptionPane.CANCEL_OPTION || option == JOptionPane.DEFAULT_OPTION) {
                return false;
            }
            return true;
        }
        return true;
    }

    public void doApplyTemplate() {
        FileDialog dialog = new FileDialog(this,
                "Apply Template",
                FileDialog.LOAD);
        dialog.setVisible(true);
        if (dialog.getFile() != null) {
            File file = new File(dialog.getDirectory(), dialog.getFile());
            try {
                readFromFile(file);
            } catch (FileNotFoundException fnfe) {
                JOptionPane.showMessageDialog(this, "Unable to open template file: File not found",
                        "Unable to open file",
                        JOptionPane.ERROR_MESSAGE);
            } catch (IOException ioe) {
                JOptionPane.showMessageDialog(this, "Unable to read template file: " + ioe.getMessage(),
                        "Unable to read file",
                        JOptionPane.ERROR_MESSAGE);
            }
        }
    }

    protected boolean readFromFile(File file) throws IOException {
        return false;
    }

    public String getDefaultFileName() {
        return beautiOptions.fileNameStem + ".beauti";
    }

    protected boolean writeToFile(File file) throws IOException {
        return false;
    }

    public final void doImport() {
        int returnVal = importChooser.showOpenDialog(this);
        if (returnVal == JFileChooser.APPROVE_OPTION) {
            File[] files = importChooser.getSelectedFiles();
            for (File file : files) {
                if (file == null || file.getName().equals("")) {
                    JOptionPane.showMessageDialog(this, "Invalid file name",
                            "Invalid file name", JOptionPane.ERROR_MESSAGE);
                } else {
                    try {
                        BEAUTiImporter beautiImporter = new BEAUTiImporter(beautiOptions);
                        beautiImporter.importFromFile(this, file);

                        setDirty();
//                    } catch (FileNotFoundException fnfe) {
//                        JOptionPane.showMessageDialog(this, "Unable to open file: File not found",
//                                "Unable to open file", JOptionPane.ERROR_MESSAGE);
                    } catch (IOException ioe) {
                        JOptionPane.showMessageDialog(this, "File I/O Error unable to read file: " + ioe.getMessage(),
                                "Unable to read file", JOptionPane.ERROR_MESSAGE);
                        ioe.printStackTrace();
                        return;

                    } catch (MissingBlockException ex) {
                        JOptionPane.showMessageDialog(this, "TAXON, DATA or CHARACTERS block is missing in Nexus file: " + ex,
                                "Missing Block in Nexus File",
                                JOptionPane.ERROR_MESSAGE);
                        ex.printStackTrace();

                    } catch (ImportException ime) {
                        JOptionPane.showMessageDialog(this, "Error parsing imported file: " + ime,
                                "Error reading file",
                                JOptionPane.ERROR_MESSAGE);
                        ime.printStackTrace();
                        return;

                    } catch (Exception ex) {
                        JOptionPane.showMessageDialog(this, "Fatal exception: " + ex,
                                "Error reading file",
                                JOptionPane.ERROR_MESSAGE);
                        ex.printStackTrace();
                        return;
                    }
                }
            }

            if (beautiOptions.allowDifferentTaxa) {
                setAllOptions();
                dataPanel.selectAll();
                dataPanel.unlinkTrees();
            }

            setStatusMessage();
            setAllOptions();

//          // @Todo templates are not implemented yet...
////        getOpenAction().setEnabled(true);
////        getSaveAction().setEnabled(true);
            getExportAction().setEnabled(true);
        }
    }

    public int allowDifferentTaxaJOptionPane() {
        // AR - Yes and No are perfectly good answers to this question
        return JOptionPane.showOptionDialog(this, "This file contains different taxa from the previously loaded\n"
                + "data partitions. This may be because the taxa are mislabelled\n" + "and need correcting before reloading.\n\n"
                + "Would you like to allow different taxa for each partition?\n", "Validation of Non-matching Taxon Name(s)",
                JOptionPane.YES_NO_OPTION, JOptionPane.WARNING_MESSAGE, null, new String[] { "Yes", "No" }, "No");
    }

    public final void doImportTraits() {
        if (beautiOptions.taxonList != null) { // validation of check empty taxonList
            FileDialog dialog = new FileDialog(this,
                    "Import Traits File...",
                    FileDialog.LOAD);

            dialog.setVisible(true);
            if (dialog.getFile() != null) {
                final File file = new File(dialog.getDirectory(), dialog.getFile());

                importTraitsFromFile(file);

                traitsPanel.fireTraitsChanged();
            }
        } else {
            JOptionPane.showMessageDialog(this, "No taxa loaded yet, please import Alignment file!",
                    "No taxa loaded", JOptionPane.ERROR_MESSAGE);
        }

    }

    protected void importTraitsFromFile(final File file) {

        try {
            importMultiTraits(file);
        } catch (FileNotFoundException fnfe) {
            JOptionPane.showMessageDialog(this, "Unable to open file: File not found",
                    "Unable to open file",
                    JOptionPane.ERROR_MESSAGE);
        } catch (IOException ioe) {
            JOptionPane.showMessageDialog(this, "Unable to read file: " + ioe.getMessage(),
                    "Unable to read file",
                    JOptionPane.ERROR_MESSAGE);
        }

        if (beautiOptions.starBEASTOptions.isSpeciesAnalysis()) { // species
            setupSpeciesAnalysis();
        }

        setAllOptions();
    }

    private void importMultiTraits(final File file) throws IOException {

        try {
            Map<String, List<String[]>> traits = Utils.importTraitsFromFile(file, "\t");

            for (Map.Entry<String, List<String[]>> e : traits.entrySet()) {
                final Class c = Utils.detectTYpe(e.getValue().get(0)[1]);
                final String traitName = e.getKey();

                Boolean warningGiven = false;
                for (String[] v : e.getValue()) {
                    final Class c1 = Utils.detectTYpe(v[1]);
                    if (c != c1 && !warningGiven) {
                        JOptionPane.showMessageDialog(this, "Not all values of same type in column" + traitName,
                                "Incompatible values", JOptionPane.WARNING_MESSAGE);
                        warningGiven = true;
                        // TODO Error - not all values of same type
                    }
                }

                TraitData.TraitType t = (c == Boolean.class || c == String.class) ? TraitData.TraitType.DISCRETE :
                        (c == Integer.class) ? TraitData.TraitType.INTEGER : TraitData.TraitType.CONTINUOUS;

                TraitData newTrait = new TraitData(traitName, file.getName(), t);

                if (validateTraitName(traitName))
                    traitsPanel.addTrait(newTrait, traitsPanel.traitsTable);

                for (final String[] v : e.getValue()) {
                    final int index = beautiOptions.taxonList.getTaxonIndex(v[0]);
                    if (index >= 0) {
                        // if the taxon isn't in the list then ignore it.
                        // TODO provide a warning of unmatched taxa
                        final Taxon taxon = beautiOptions.taxonList.getTaxon(index);
                        taxon.setAttribute(traitName, Utils.constructFromString(c, v[1]));
                    }
                }
            }
        } catch (Arguments.ArgumentException e) {
            JOptionPane.showMessageDialog(this, "Error in loading traits file: " + e.getMessage(),
                    "Error Loading file", JOptionPane.ERROR_MESSAGE);
        }
    }

    public boolean validateTraitName(String traitName) {
        // check that the name is valid
        if (traitName.trim().length() == 0) {
            Toolkit.getDefaultToolkit().beep();
            return false;
        }

        // disallow a trait called 'date'
        if (traitName.equalsIgnoreCase("date")) {
            JOptionPane.showMessageDialog(this,
                    "This trait name has a special meaning. Use the 'Tip Date' panel\n" +
                            " to set dates for taxa.",
                    "Reserved trait name",
                    JOptionPane.WARNING_MESSAGE);

            return false;
        }

        // check that the trait name doesn't exist
        if (BeautiOptions.containTrait(traitName)) {
            int option = JOptionPane.showConfirmDialog(this,
                    "A trait of this name already exists. Do you wish to replace\n" +
                            "it with this new trait? This may result in the loss or change\n" +
                            "in trait values for the taxa.",
                    "Overwrite trait?",
                    JOptionPane.YES_NO_OPTION,
                    JOptionPane.WARNING_MESSAGE);

            if (option == JOptionPane.NO_OPTION) {
                return false;
            }
        }

        return true;
    }

    public void setupSpeciesAnalysis() {
        dataPanel.selectAll();
        dataPanel.unlinkAll();
//        dataPanel.unlinkModels();
//        dataPanel.unlinkTrees();

//        if (beautiOptions.getPartitionClockModels().size() > 1) {
//        	dataPanel.linkClocks();
//        }

//        beautiOptions.rateOptionClockModel = FixRateType.FIX_FIRST_PARTITION;
//        beautiOptions.activedSameTreePrior.setNodeHeightPrior(TreePriorType.SPECIES_YULE);
//
//        int i = tabbedPane.indexOfTab("Trees");
//        tabbedPane.removeTabAt(i);
//        tabbedPane.insertTab("Trees", null, speciesTreesPanel, "", i);
//        speciesTreesPanel.getOptions(beautiOptions);

        treesPanel.updatePriorPanelForSpeciesAnalysis();

        beautiOptions.starBEASTOptions = new STARBEASTOptions(beautiOptions);
        beautiOptions.fileNameStem = "LogStem";

        setStatusMessage();
    }

    public void removeSepciesAnalysis() {
//        beautiOptions.activedSameTreePrior.setNodeHeightPrior(TreePriorType.CONSTANT);
//
//        int i = tabbedPane.indexOfTab("Trees");
//        tabbedPane.removeTabAt(i);
//        if (DataPanel.ALLOW_UNLINKED_TREES) {
//            tabbedPane.insertTab("Trees", null, treesPanel, "", i);
//        } else {
//            tabbedPane.insertTab("Trees", null, oldTreesPanel, "", i);
//        }

    	treesPanel.updatePriorPanelForSpeciesAnalysis();

        setStatusMessage();
    }

    public void updateDiscreteTraitAnalysis() {
    	setStatusMessage();
    }

    public void setupEBSP() {
    	dataPanel.selectAll();
    	dataPanel.unlinkAll();

    	setAllOptions();
    }

    public PartitionTreePrior getCurrentPartitionTreePrior() {
    	return treesPanel.currentTreeModel.getPartitionTreePrior();
    }

    public void removeSpecifiedTreePrior(boolean isChecked) { // TipDatesPanel usingTipDates
    	//TODO: wait for new implementation in BEAST
	    if (DataPanel.ALLOW_UNLINKED_TREES) {
			treesPanel.setCheckedTipDate(isChecked);
		} else {
			if (isChecked) {
				oldTreesPanel.treePriorCombo.removeItem(TreePriorType.YULE);
				oldTreesPanel.treePriorCombo.removeItem(TreePriorType.BIRTH_DEATH);
			} else {
				oldTreesPanel.treePriorCombo = new JComboBox(EnumSet.range(TreePriorType.CONSTANT, TreePriorType.BIRTH_DEATH).toArray());
			}
		}
    }

    public void setStatusMessage() {
        GUIValidate();

        statusLabel.setText(beautiOptions.statusMessage());
    }

     public void GUIValidate() {
         if (beautiOptions.starBEASTOptions.isSpeciesAnalysis()) {
             if (beautiOptions.starBEASTOptions.getSpeciesList() == null) {
                 JOptionPane.showMessageDialog(this, "Species value is empty."
                         + "\nPlease go to Traits panel, either Import Traits,"
                         + "\nor Guess trait values", "*BEAST Error Message",
                    JOptionPane.ERROR_MESSAGE);
             }
         }
    }

    public final boolean doGenerate() {

        try {
            generator.checkOptions();
        } catch (IllegalArgumentException iae) {
            JOptionPane.showMessageDialog(this, iae.getMessage(), "Unable to generate file",
                    JOptionPane.ERROR_MESSAGE);
            return false;
        }

        DefaultPriorDialog defaultPriorDialog = new DefaultPriorDialog(this);
        if (!defaultPriorDialog.showDialog(beautiOptions)) {
           return false;
        }

        // offer stem as default
        exportChooser.setSelectedFile(new File(beautiOptions.fileNameStem + ".xml"));

        final int returnVal = exportChooser.showSaveDialog(this);
        if( returnVal == JFileChooser.APPROVE_OPTION ) {
            File file = exportChooser.getSelectedFile();

            int n = JOptionPane.YES_OPTION;

            if (file.exists()) {
            	n = JOptionPane.showConfirmDialog(this, file.getName(),
            		     "Overwrite the existing file?", JOptionPane.YES_NO_OPTION);
            }

            if (n == JOptionPane.YES_OPTION) {
	            try {
	                getAllOptions();
                    generator.generateXML(file);

	            } catch (IOException ioe) {
	                JOptionPane.showMessageDialog(this, "Unable to generate file: " + ioe.getMessage(),
	                        "Unable to generate file", JOptionPane.ERROR_MESSAGE);
	                return false;
	            } catch (Exception e) {
	                JOptionPane.showMessageDialog(this, "Unable to generate file: " + e.getMessage(),
	                        "Unable to generate file", JOptionPane.ERROR_MESSAGE);
	                return false;
	            }
            } else {
            	doGenerate();
            }
        }

        clearDirty();
        return true;
    }

    public JComponent getExportableComponent() {

        JComponent exportable = null;
        Component comp = tabbedPane.getSelectedComponent();

        if (comp instanceof Exportable) {
            exportable = ((Exportable) comp).getExportableComponent();
        } else if (comp instanceof JComponent) {
            exportable = (JComponent) comp;
        }

        return exportable;
    }

    public boolean doSave() {
        return doSaveAs();
    }

    public boolean doSaveAs() {
        FileDialog dialog = new FileDialog(this,
                "Save Template As...",
                FileDialog.SAVE);

        dialog.setVisible(true);
        if (dialog.getFile() == null) {
            // the dialog was cancelled...
            return false;
        }

        File file = new File(dialog.getDirectory(), dialog.getFile());

        try {
            if (writeToFile(file)) {

                clearDirty();
            }
        } catch (IOException ioe) {
            JOptionPane.showMessageDialog(this, "Unable to save file: " + ioe,
                    "Unable to save file",
                    JOptionPane.ERROR_MESSAGE);
        }

        return true;
    }

    public Action getOpenAction() {
        return openTemplateAction;
    }

    private final AbstractAction openTemplateAction = new AbstractAction("Apply Template...") {
        private static final long serialVersionUID = 2450459627280385426L;

        public void actionPerformed(ActionEvent ae) {
            doApplyTemplate();
        }
    };

    public Action getSaveAction() {
        return saveAsAction;
    }

    public Action getSaveAsAction() {
        return saveAsAction;
    }

    private final AbstractAction saveAsAction = new AbstractAction("Save Template As...") {
        private static final long serialVersionUID = 2424923366448459342L;

        public void actionPerformed(ActionEvent ae) {
            doSaveAs();
        }
    };

    public Action getImportAction() {
        return importAlignmentAction;
    }

    protected AbstractAction importAlignmentAction = new AbstractAction("Import Alignment(s)") {
        private static final long serialVersionUID = 3217702096314745005L;

        public void actionPerformed(java.awt.event.ActionEvent ae) {
            doImport();
        }
    };

    public Action getImportTraitsAction() {
        return importTraitsAction;
    }

    protected AbstractAction importTraitsAction = new AbstractAction("Import Trait(s)") {
        private static final long serialVersionUID = 3217702096314745005L;

        public void actionPerformed(java.awt.event.ActionEvent ae) {
            doImportTraits();
        }
    };

    public Action getExportAction() {
        return generateAction;
    }

    protected AbstractAction generateAction = new AbstractAction("Generate BEAST File...", gearIcon) {
        private static final long serialVersionUID = -5329102618630268783L;

        public void actionPerformed(java.awt.event.ActionEvent ae) {
            doGenerate();
        }
    };

}
