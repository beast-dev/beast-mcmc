package dr.app.phylogeography.spread;

import dr.app.gui.DeleteActionResponder;
import dr.app.gui.SelectAllActionResponder;
import dr.app.java16compat.FileNameExtensionFilter;
import dr.app.phylogeography.generator.Generator;
import dr.app.phylogeography.generator.KMLGenerator;
import dr.app.phylogeography.spread.layerspanel.LayersPanel;
import dr.app.phylogeography.spread.inputpanel.InputPanel;
import dr.app.util.Utils;
import dr.evolution.io.Importer;
import dr.evolution.io.NexusImporter;
import dr.evolution.tree.Tree;
import dr.evolution.util.TaxonList;
import jam.framework.DocumentFrame;
import jam.framework.Exportable;
import jam.util.IconUtils;

import javax.swing.*;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import javax.swing.plaf.BorderUIResource;
import java.awt.*;
import java.io.*;
import java.util.ArrayList;

/**
 * @author Andrew Rambaut
 * @version $Id$
 */
public class SpreadFrame extends DocumentFrame {
    private static final long serialVersionUID = 2114148696789612509L;

    private final SpreadDocument document = new SpreadDocument();
//    private final BeastGenerator generator;

    private final JTabbedPane tabbedPane = new JTabbedPane();
    private final JLabel statusLabel = new JLabel("No data loaded");
    private final InputPanel inputPanel;
    private final LayersPanel layersPanel = new LayersPanel(this, document);

    private JFileChooser importChooser; // make JFileChooser chooser remember previous path
    private JFileChooser exportChooser; // make JFileChooser chooser remember previous path

    private final Icon gearIcon = IconUtils.getIcon(this.getClass(), "images/gear.png");

    private final java.util.List<Generator> generators;

    public SpreadFrame(String title) {
        super();

        generators = new ArrayList<Generator>();
        generators.add(new KMLGenerator());

        setTitle(title);

        // Prevent the application to close in requestClose()
        // after a user cancel or a failure in beast file generation
        setDefaultCloseOperation(DO_NOTHING_ON_CLOSE);

//        getOpenAction().setEnabled(false);
//        getSaveAction().setEnabled(false);

        getFindAction().setEnabled(false);

        getZoomWindowAction().setEnabled(false);

        inputPanel = new InputPanel(this, document, getImportAction());

        tabbedPane.addChangeListener(new ChangeListener() {
            public void stateChanged(final ChangeEvent e) {
                tabbedPaneChanged();
            }
        });

        document.addListener(new SpreadDocument.Listener() {
            public void dataChanged() {
                setDirty();
            }

            public void settingsChanged() {
                setDirty();
            }
        });
    }

    public void initializeComponents() {

        final TimelinePanel timeLinePanel = new TimelinePanel(this, document);
        final OutputPanel outputPanel = new OutputPanel(this, document, generators);

        tabbedPane.addTab("Input", inputPanel);
        tabbedPane.addTab("Layers", layersPanel);
        tabbedPane.addTab("Timeline", timeLinePanel);
        tabbedPane.addTab("Output", outputPanel);

        JPanel panel = new JPanel(new BorderLayout(6, 6));
        panel.setBorder(new BorderUIResource.EmptyBorderUIResource(new java.awt.Insets(12, 12, 12, 12)));
        panel.add(tabbedPane, BorderLayout.CENTER);
        panel.setPreferredSize(new java.awt.Dimension(800, 600));

//        getExportAction().setEnabled(false);
//        JButton generateButton = new JButton(getExportAction());
//        generateButton.putClientProperty("JButton.buttonType", "roundRect");

        JPanel panel2 = new JPanel(new BorderLayout(6, 6));
        panel2.add(statusLabel, BorderLayout.CENTER);
//        panel2.add(generateButton, BorderLayout.EAST);

        panel.add(panel2, BorderLayout.SOUTH);

        getContentPane().setLayout(new java.awt.BorderLayout(0, 0));
        getContentPane().add(panel2, BorderLayout.CENTER);

//        setAllOptions();

        setSize(new java.awt.Dimension(1024, 768));
        setMinimumSize(new java.awt.Dimension(800, 600));

        // make JFileChooser chooser remember previous path
        exportChooser = new JFileChooser(Utils.getCWD());
        exportChooser.setDialogTitle("Generate Map File...");


        importChooser = new JFileChooser(Utils.getCWD());

        importChooser.setMultiSelectionEnabled(true);
        importChooser.setFileFilter(new FileNameExtensionFilter(
                "NEXUS (*.nex) & BEAST Tree (*.trees) Files", "nex", "nexus", "nx", "tree", "tre", "trees"));
        importChooser.setDialogTitle("Import Data Files...");
    }

    public final void tabbedPaneChanged() {
        if (tabbedPane.getSelectedComponent() instanceof DeleteActionResponder) {
            getDeleteAction().setEnabled(
                    ((DeleteActionResponder)(tabbedPane.getSelectedComponent()))
                            .getDeleteAction()
                            .isEnabled()
            );
        }
    }

    public void doDelete() {
        if (tabbedPane.getSelectedComponent() instanceof DeleteActionResponder) {
            ((DeleteActionResponder)(tabbedPane.getSelectedComponent())).delete();
        }
    }

    public void doSelectAll() {
        if (tabbedPane.getSelectedComponent() instanceof SelectAllActionResponder) {
            ((SelectAllActionResponder)(tabbedPane.getSelectedComponent())).selectAll();
        }
    }


    public void setRemoveActionEnabled(final DeleteActionResponder responder, final boolean isEnabled) {
        if (responder == tabbedPane.getSelectedComponent()) {
            getDeleteAction().setEnabled(isEnabled);
        }
    }


    public boolean requestClose() {
        if (isDirty()) {
            int option = JOptionPane.showConfirmDialog(this,
                    "You have made changes but have not generated\n" +
                            "an output file. Do you wish to generate\n" +
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

    protected boolean readFromFile(File file) throws IOException {
        return false;
    }

    public String getDefaultFileName() {
        return "untitled";
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
                        importDataFile(file);

                        setDirty();
//                    } catch (FileNotFoundException fnfe) {
//                        JOptionPane.showMessageDialog(this, "Unable to open file: File not found",
//                                "Unable to open file", JOptionPane.ERROR_MESSAGE);
                    } catch (IOException ioe) {
                        JOptionPane.showMessageDialog(this, "File I/O Error unable to read file: " + ioe.getMessage(),
                                "Unable to read file", JOptionPane.ERROR_MESSAGE);
                        return;

                    } catch (NexusImporter.MissingBlockException ex) {
                        JOptionPane.showMessageDialog(this, "TAXON, DATA or CHARACTERS block is missing in Nexus file: " + ex,
                                "Missing Block in Nexus File",
                                JOptionPane.ERROR_MESSAGE);

                    } catch (Importer.ImportException ime) {
                        JOptionPane.showMessageDialog(this, "Error parsing imported file: " + ime,
                                "Error reading file",
                                JOptionPane.ERROR_MESSAGE);
                        return;

                    } catch (Exception ex) {
                        JOptionPane.showMessageDialog(this, "Fatal exception: " + ex,
                                "Error reading file",
                                JOptionPane.ERROR_MESSAGE);
                        return;
                    }
                }
            }

//            setAllOptions();
        }
    }

    public void importDataFile(File file) throws IOException, Importer.ImportException {
        Reader reader = new FileReader(file);

        BufferedReader bufferedReader = new BufferedReader(reader);
        String line = bufferedReader.readLine();
        while (line != null && line.length() == 0) {
            line = bufferedReader.readLine();
        }

        if ((line != null && line.toUpperCase().contains("#NEXUS"))) {
            int treeCount = 0;

            do {
                line = bufferedReader.readLine();
            } while (line != null && !line.startsWith("tree STATE"));

            while (line != null) {
                if (line.startsWith("tree STATE")) {
                    treeCount++;
                }
                line = bufferedReader.readLine();
            }

            // is a NEXUS file
            Tree tree = importFirstTree(file);
            if (tree != null) {
                if (treeCount > 0) {
                    InputFile inputFile = new InputFile(file, tree, treeCount);
                    document.addTreeFile(inputFile);
                } else {
                    InputFile inputFile = new InputFile(file, tree);
                    document.addTreeFile(inputFile);
                }
            }
        } else {
            JOptionPane.showMessageDialog(this, "Error parsing imported file. This may not be a NEXUS file",
                    "Error reading file",
                    JOptionPane.ERROR_MESSAGE);
        }
    }

    // nexus
    private Tree importFirstTree(File file) throws IOException, Importer.ImportException {
        TaxonList taxa = null;
        Tree tree = null;

        FileReader reader = new FileReader(file);

        NexusImporter importer = new NexusImporter(reader);

        tree = importer.importTree(taxa);


        return tree;
    }


    public void setStatusMessage(String text) {
        statusLabel.setText(text);
    }

    public final boolean doGenerate() {

        try {
//            generator.checkOptions();
        } catch (IllegalArgumentException iae) {
            JOptionPane.showMessageDialog(this, iae.getMessage(), "Unable to generate file",
                    JOptionPane.ERROR_MESSAGE);
            return false;
        }

        // offer stem as default
//        exportChooser.setSelectedFile(new File(beautiOptions.fileNameStem + ".xml"));

        final int returnVal = exportChooser.showSaveDialog(this);
        if( returnVal == JFileChooser.APPROVE_OPTION ) {
            File file = exportChooser.getSelectedFile();

            int n = JOptionPane.YES_OPTION;

            if (file.exists()) {
                n = JOptionPane.showConfirmDialog(this, file.getName(),
                        "Overwrite the exsting file?", JOptionPane.YES_NO_OPTION);
            }

            if (n == JOptionPane.YES_OPTION) {
                try {
                    generate(file);

                } catch (IOException ioe) {
                    JOptionPane.showMessageDialog(this, "Unable to generate file: " + ioe.getMessage(),
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

    protected void generate(File file) throws IOException {
//        getAllOptions();

        FileWriter fw = new FileWriter(file);
//        generator.generateXML(fw);
        fw.close();
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

    public Action getImportAction() {
        return importAlignmentAction;
    }

    protected AbstractAction importAlignmentAction = new AbstractAction("Import Alignment...") {
        private static final long serialVersionUID = 3217702096314745005L;

        public void actionPerformed(java.awt.event.ActionEvent ae) {
            doImport();
        }
    };

//    public Action getImportTraitsAction() {
//        return importTraitsAction;
//    }
//
//    protected AbstractAction importTraitsAction = new AbstractAction("Import Traits...") {
//        private static final long serialVersionUID = 3217702096314745005L;
//
//        public void actionPerformed(java.awt.event.ActionEvent ae) {
//            doImportTraits();
//        }
//    };

    public Action getExportAction() {
        return generateAction;
    }

    protected AbstractAction generateAction = new AbstractAction("Generate Map File...", gearIcon) {
        private static final long serialVersionUID = -5329102618630268783L;

        public void actionPerformed(java.awt.event.ActionEvent ae) {
            doGenerate();
        }
    };



}
