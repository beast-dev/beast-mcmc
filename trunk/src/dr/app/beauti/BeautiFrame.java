/*
 * BeautiFrame.java
 *
 * (c) 2002-2005 BEAST Development Core Team
 *
 * This package may be distributed under the
 * Lesser Gnu Public Licence (LGPL)
 */
package dr.app.beauti;

import dr.app.beauti.clockModelsPanel.OldClockModelsPanel;
import dr.app.beauti.components.ComponentFactory;
import dr.app.beauti.components.SequenceErrorModelComponentFactory;
import dr.app.beauti.components.TipDateSamplingComponentFactory;
import dr.app.beauti.datapanel.DataPanel;
import dr.app.beauti.generator.BeastGenerator;
import dr.app.beauti.generator.Generator;
import dr.app.beauti.mcmcpanel.MCMCPanel;
import dr.app.beauti.operatorspanel.OperatorsPanel;
import dr.app.beauti.options.BeautiOptions;
import dr.app.beauti.options.PartitionTreePrior;
import dr.app.beauti.options.STARBEASTOptions;
import dr.app.beauti.options.TraitData;
import dr.app.beauti.priorsPanel.DefaultPriorDialog;
import dr.app.beauti.priorsPanel.PriorsPanel;
import dr.app.beauti.siteModelsPanel.SiteModelsPanel;
import dr.app.beauti.taxonsetspanel.TaxaPanel;
import dr.app.beauti.tipdatepanel.TipDatesPanel;
import dr.app.beauti.traitspanel.TraitsPanel;
import dr.app.beauti.treespanel.TreesPanel;
import dr.app.beauti.util.BEAUTiImporter;
import dr.app.beauti.util.TextUtil;
import dr.app.gui.FileDrop;
import dr.app.util.OSType;
import dr.app.util.Utils;
import dr.evolution.io.Importer.ImportException;
import dr.evolution.io.NexusImporter.MissingBlockException;
import jam.framework.DocumentFrame;
import jam.framework.Exportable;
import jam.util.IconUtils;

import javax.swing.*;
import javax.swing.border.Border;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import javax.swing.filechooser.FileNameExtensionFilter;
import javax.swing.plaf.BorderUIResource;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.HierarchyBoundsListener;
import java.awt.event.HierarchyEvent;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: BeautiFrame.java,v 1.22 2006/09/09 16:07:06 rambaut Exp $
 */
public class BeautiFrame extends DocumentFrame {

    private static final long serialVersionUID = 2114148696789612509L;

    private final BeautiOptions options;
    private final BeastGenerator generator;

    public final JTabbedPane tabbedPane = new JTabbedPane();
    public final JLabel statusLabel = new JLabel();

    private DataPanel dataPanel;
    private TipDatesPanel tipDatesPanel;
    private TraitsPanel traitsPanel;
    private TaxaPanel taxaPanel;
    private SiteModelsPanel siteModelsPanel;
    private OldClockModelsPanel clockModelsPanel;
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

        options = new BeautiOptions(components);
        generator = new BeastGenerator(options, components);

        this.getContentPane().addHierarchyBoundsListener(new HierarchyBoundsListener() {
            public void ancestorMoved(HierarchyEvent e) {
            }

            public void ancestorResized(HierarchyEvent e) {
                setStatusMessage();
            }
        });
    }

    public void initializeComponents() {

        dataPanel = new DataPanel(this, getImportAction(), getDeleteAction()/*, getImportTraitsAction()*/);
        tipDatesPanel = new TipDatesPanel(this);
        traitsPanel = new TraitsPanel(this, getImportTraitsAction());
        taxaPanel = new TaxaPanel(this);
        siteModelsPanel = new SiteModelsPanel(this, getDeleteAction());
        clockModelsPanel = new OldClockModelsPanel(this);
//        oldTreesPanel = new OldTreesPanel(this);
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
        tabbedPane.addTab("Trees", treesPanel);
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
                currentPanel.getOptions(options);
                setAllOptions();
                currentPanel = selectedPanel;
            }
        });

        JPanel basePanel = new JPanel(new BorderLayout(6, 6));
        basePanel.setBorder(new BorderUIResource.EmptyBorderUIResource(new java.awt.Insets(12, 12, 12, 12)));
//        basePanel.setPreferredSize(new java.awt.Dimension(800, 600));

        getExportAction().setEnabled(false);
        JButton generateButton = new JButton(getExportAction());
        generateButton.putClientProperty("JButton.buttonType", "roundRect");

        JPanel panel2 = new JPanel(new BorderLayout(6, 6));
        panel2.add(statusLabel, BorderLayout.WEST);
        panel2.add(generateButton, BorderLayout.EAST);
        panel2.setMinimumSize(new java.awt.Dimension(10, 10));

        basePanel.add(tabbedPane, BorderLayout.CENTER);
        basePanel.add(panel2, BorderLayout.SOUTH);

        add(basePanel, BorderLayout.CENTER);

        Toolkit tk = Toolkit.getDefaultToolkit();
        Dimension d = tk.getScreenSize();
        System.out.println("Screen width = " + d.width);
        System.out.println("Screen height = " + d.height);

        if (d.width < 1000 || d.height < 700) {
            setSize(new java.awt.Dimension(700, 550));
        } else {
            setSize(new java.awt.Dimension(1024, 768));
        }

        // todo Andrew:
        // It is really nasty on the Mac not to have a minimum window size (probably
        // other platforms too?). We surely don't require a working interface less
        // than this...
        if (OSType.isMac()) {
            setMinimumSize(new java.awt.Dimension(640, 480));
        }

        setAllOptions();

        // make JFileChooser chooser remember previous path
        exportChooser = new JFileChooser(Utils.getCWD());
        exportChooser.setFileFilter(new FileNameExtensionFilter("BEAST XML File", "xml", "beast"));
        exportChooser.setDialogTitle("Generate BEAST XML File...");


        importChooser = new JFileChooser(Utils.getCWD());

        importChooser.setMultiSelectionEnabled(true);
        importChooser.setFileFilter(new FileNameExtensionFilter(
                        "Microsatellite (tab-delimited *.txt) Files", "txt"));
        importChooser.setFileFilter(new FileNameExtensionFilter(
                "NEXUS (*.nex) & BEAST (*.xml) Files", "nex", "nexus", "nx", "xml", "beast", "fa", "fasta", "afa"));
        importChooser.setDialogTitle("Import Aligment...");

        Color focusColor = UIManager.getColor("Focus.color");
        Border focusBorder = BorderFactory.createMatteBorder(2, 2, 2, 2, focusColor);
        dataPanel.setBorder(BorderFactory.createEmptyBorder(2, 2, 2, 2));
        new FileDrop(null, dataPanel, focusBorder, new FileDrop.Listener() {
            public void filesDropped(java.io.File[] files) {
                importFiles(files);
            }   // end filesDropped
        }); // end FileDrop.Listener


    }

    /**
     * set all the options for all panels
     */
    public void setAllOptions() {
        try {
            dataPanel.setOptions(options);
            tipDatesPanel.setOptions(options);
            traitsPanel.setOptions(options);
            taxaPanel.setOptions(options);
            siteModelsPanel.setOptions(options);
            clockModelsPanel.setOptions(options);
            treesPanel.setOptions(options);
            priorsPanel.setOptions(options);
            operatorsPanel.setOptions(options);
            mcmcPanel.setOptions(options);

            setStatusMessage();
        } catch (IllegalArgumentException illegEx) {
            JOptionPane.showMessageDialog(this, illegEx.getMessage(),
                    "Illegal Argument Exception", JOptionPane.ERROR_MESSAGE);
        }
    }

    /**
     * get all the options for all panels
     */
    private void getAllOptions() {
        try {
            dataPanel.getOptions(options);
            tipDatesPanel.getOptions(options);
            traitsPanel.getOptions(options);
            taxaPanel.getOptions(options);
            siteModelsPanel.getOptions(options);
            clockModelsPanel.getOptions(options);
            treesPanel.getOptions(options);
            priorsPanel.getOptions(options);
            operatorsPanel.getOptions(options);
            mcmcPanel.getOptions(options);
        } catch (IllegalArgumentException illegEx) {
            JOptionPane.showMessageDialog(this, illegEx.getMessage(),
                    "Illegal Argument Exception", JOptionPane.ERROR_MESSAGE);
        }
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
//            modelsPanel.delete();
//        } else if (tabbedPane.getSelectedComponent() == treesPanel) {
//        	treesPanel.delete();
        } else {
            throw new RuntimeException("Delete should only be accessable from the Data and Models panels");
        }

        setStatusMessage();
    }

    public boolean requestClose() {
        if (isDirty() && options.hasData()) {
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
        return options.fileNameStem + ".beauti";
    }

    protected boolean writeToFile(File file) throws IOException {
        return false;
    }

    public final void doImport() {
        int returnVal = importChooser.showOpenDialog(this);
        if (returnVal == JFileChooser.APPROVE_OPTION) {
            File[] files = importChooser.getSelectedFiles();
            importFiles(files);
        }
    }

    private void importFiles(File[] files) {
        for (File file : files) {
            if (file == null || file.getName().equals("")) {
                JOptionPane.showMessageDialog(this, "Invalid file name",
                        "Invalid file name", JOptionPane.ERROR_MESSAGE);
            } else {
                try {
                    BEAUTiImporter beautiImporter = new BEAUTiImporter(this, options);
                    beautiImporter.importFromFile(file);

                    setDirty();
//                    } catch (FileNotFoundException fnfe) {
//                        JOptionPane.showMessageDialog(this, "Unable to open file: File not found",
//                                "Unable to open file", JOptionPane.ERROR_MESSAGE);
                } catch (IOException ioe) {
                    JOptionPane.showMessageDialog(this, "File I/O Error unable to read file: " + ioe.getMessage(),
                            "Unable to read file", JOptionPane.ERROR_MESSAGE);
                    ioe.printStackTrace();
                    // there may be other files in the list so don't return
//                    return;

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
                    // there may be other files in the list so don't return
//                    return;
                } catch (IllegalArgumentException illegEx) {
                    JOptionPane.showMessageDialog(this, illegEx.getMessage(),
                            "Illegal Argument Exception", JOptionPane.ERROR_MESSAGE);

                } catch (Exception ex) {
                    JOptionPane.showMessageDialog(this, "Fatal exception: " + ex,
                            "Error reading file",
                            JOptionPane.ERROR_MESSAGE);
                    ex.printStackTrace();
                    return;
                }
            }
        }

        if (options.allowDifferentTaxa) {
            setAllOptions();
            dataPanel.selectAll();
            dataPanel.unlinkTrees();
        }

        setAllOptions();

//          // @Todo templates are not implemented yet...
////        getOpenAction().setEnabled(true);
////        getSaveAction().setEnabled(true);
        getExportAction().setEnabled(true);
    }

    public final boolean doImportTraits() {
        if (options.taxonList != null) { // validation of check empty taxonList
            FileDialog dialog = new FileDialog(this,
                    "Import Traits File...",
                    FileDialog.LOAD);

            dialog.setVisible(true);
            if (dialog.getFile() != null) {
                final File file = new File(dialog.getDirectory(), dialog.getFile());

                try {
                    BEAUTiImporter beautiImporter = new BEAUTiImporter(this, options);
                    beautiImporter.importTraits(file);
                } catch (FileNotFoundException fnfe) {
                    JOptionPane.showMessageDialog(this, "Unable to open file: File not found",
                            "Unable to open file",
                            JOptionPane.ERROR_MESSAGE);
                    return false;
                } catch (IOException ioe) {
                    JOptionPane.showMessageDialog(this, "Unable to read file: " + ioe.getMessage(),
                            "Unable to read file",
                            JOptionPane.ERROR_MESSAGE);
                    return false;
                } catch (Exception ex) {
                    JOptionPane.showMessageDialog(this, "Fatal exception: " + ex,
                            "Error reading file",
                            JOptionPane.ERROR_MESSAGE);
                    ex.printStackTrace();
                    return false;
                }
            } else {
                return false;
            }

            traitsPanel.fireTraitsChanged();
            setAllOptions();
            return true;

        } else {
            JOptionPane.showMessageDialog(this, "No taxa loaded yet, please import Alignment file.",
                    "No taxa loaded", JOptionPane.ERROR_MESSAGE);
            return false;
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
        if (options.traitExists(traitName)) {
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

    public void setupStarBEAST(boolean useStarBEAST) {
        if (useStarBEAST) {
            dataPanel.selectAll();
            dataPanel.unlinkAll();

            options.starBEASTOptions = new STARBEASTOptions(options);
            options.fileNameStem = "StarBEASTLog";

            if (!options.traitExists(TraitData.TRAIT_SPECIES)) {
                if (!traitsPanel.addTrait(TraitData.TRAIT_SPECIES)) {
                    dataPanel.useStarBEASTCheck.setSelected(false); // go back to unchecked
                    useStarBEAST = false;
                }

                // why delete this? The user may want to use it again
//        } else { // remove species
//            options.removeTrait(TraitData.TRAIT_SPECIES);
            }
        }

        options.useStarBEAST = useStarBEAST;

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
        treesPanel.setCheckedTipDate(isChecked);
    }

    public void setStatusMessage() {
        int width = this.getWidth() - 260; // minus generate button size
        if (width < 100) width = 100; // prevent too narrow
        String tw = TextUtil.wrapText(options.statusMessage(), statusLabel, width);
//        System.out.println(this.getWidth() + "   " + tw);
        statusLabel.setText(tw);
    }

    public final boolean doGenerate() {

        try {
            generator.checkOptions();
        } catch (IllegalArgumentException iae) {
            JOptionPane.showMessageDialog(this, iae.getMessage(), "Invalid BEAUti setting : ",
                    JOptionPane.ERROR_MESSAGE);
            return false;
        }

        DefaultPriorDialog defaultPriorDialog = new DefaultPriorDialog(this);
        if (!defaultPriorDialog.showDialog(options)) {
            return false;
        }

        // offer stem as default
        exportChooser.setSelectedFile(new File(options.fileNameStem + ".xml"));

        final int returnVal = exportChooser.showSaveDialog(this);
        if (returnVal == JFileChooser.APPROVE_OPTION) {
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
                    JOptionPane.showMessageDialog(this, "Unable to generate file due to I/O issue: " + ioe.getMessage(),
                            "Unable to generate file", JOptionPane.ERROR_MESSAGE);
                    return false;
                } catch (Generator.GeneratorException e) {
                    JOptionPane.showMessageDialog(this, "The BEAST XML is incomplete because :\n" + e.getMessage(),
                            "The BEAST XML is incomplete", JOptionPane.ERROR_MESSAGE);
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

    protected AbstractAction importAlignmentAction = new AbstractAction("Import Data...") {
        private static final long serialVersionUID = 3217702096314745005L;

        public void actionPerformed(java.awt.event.ActionEvent ae) {
            doImport();
        }
    };

    public Action getImportTraitsAction() {
        return importTraitsAction;
    }

    protected AbstractAction importTraitsAction = new AbstractAction("Import Traits") {
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
