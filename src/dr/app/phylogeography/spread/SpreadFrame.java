package dr.app.phylogeography.spread;

import dr.app.util.Arguments;
import dr.app.util.Utils;
import dr.evolution.io.NexusImporter.MissingBlockException;
import dr.evolution.io.Importer.ImportException;
import dr.evolution.util.Taxon;
import dr.app.java16compat.FileNameExtensionFilter;
import jam.framework.DocumentFrame;
import jam.framework.Exportable;
import jam.util.IconUtils;

import javax.swing.*;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import javax.swing.plaf.BorderUIResource;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.io.*;
import java.util.EnumSet;
import java.util.List;
import java.util.Map;

/**
 * @author Andrew Rambaut
 * @version $Id: BeautiFrame.java,v 1.22 2006/09/09 16:07:06 rambaut Exp $
 */
public class SpreadFrame extends DocumentFrame {

    private static final long serialVersionUID = 2114148696789612509L;

//    private final BeautiOptions beautiOptions;
//    private final BeastGenerator generator;

    private final JTabbedPane tabbedPane = new JTabbedPane();
    public final JLabel statusLabel = new JLabel("No data loaded");

    private JFileChooser importChooser; // make JFileChooser chooser remember previous path
    private JFileChooser exportChooser; // make JFileChooser chooser remember previous path

    final Icon gearIcon = IconUtils.getIcon(this.getClass(), "images/gear.png");

    public SpreadFrame(String title) {
        super();

        setTitle(title);

        // Prevent the application to close in requestClose()
        // after a user cancel or a failure in beast file generation
        setDefaultCloseOperation(DO_NOTHING_ON_CLOSE);

//        getOpenAction().setEnabled(false);
//        getSaveAction().setEnabled(false);

        getFindAction().setEnabled(false);

        getZoomWindowAction().setEnabled(false);
    }

    public void initializeComponents() {

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

//        setAllOptions();

        setSize(new java.awt.Dimension(1024, 768));

        // make JFileChooser chooser remember previous path
        exportChooser = new JFileChooser(Utils.getCWD());
        exportChooser.setFileFilter(new FileNameExtensionFilter("BEAST XML File", "xml", "beast"));
        exportChooser.setDialogTitle("Generate BEAST XML File...");


        importChooser = new JFileChooser(Utils.getCWD());

        importChooser.setMultiSelectionEnabled(true);
        importChooser.setFileFilter(new FileNameExtensionFilter(
                "NEXUS (*.nex) & BEAST (*.xml) Files", "nex", "nexus", "nx", "xml", "beast"));
        importChooser.setDialogTitle("Import Aligment...");
    }

    public final void dataSelectionChanged(boolean isSelected) {
        getDeleteAction().setEnabled(isSelected);
    }

    public final void modelSelectionChanged(boolean isSelected) {
        getDeleteAction().setEnabled(isSelected);
    }

    public void doDelete() {
        setStatusMessage();
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
//        int returnVal = importChooser.showOpenDialog(this);
//        if (returnVal == JFileChooser.APPROVE_OPTION) {
//            File[] files = importChooser.getSelectedFiles();
//            for (File file : files) {
//                if (file == null || file.getName().equals("")) {
//                    JOptionPane.showMessageDialog(this, "Invalid file name",
//                            "Invalid file name", JOptionPane.ERROR_MESSAGE);
//                } else {
//                    try {
////                        beautiOptions.beautiImporter.importFromFile(this, file);
//
//                        setDirty();
////                    } catch (FileNotFoundException fnfe) {
////                        JOptionPane.showMessageDialog(this, "Unable to open file: File not found",
////                                "Unable to open file", JOptionPane.ERROR_MESSAGE);
//                    } catch (IOException ioe) {
//                        JOptionPane.showMessageDialog(this, "File I/O Error unable to read file: " + ioe.getMessage(),
//                                "Unable to read file", JOptionPane.ERROR_MESSAGE);
//                        ioe.printStackTrace();
//                        return;
//
//                    } catch (MissingBlockException ex) {
//                        JOptionPane.showMessageDialog(this, "TAXON, DATA or CHARACTERS block is missing in Nexus file: " + ex,
//                                "Missing Block in Nexus File",
//                                JOptionPane.ERROR_MESSAGE);
//                        ex.printStackTrace();
//
//                    } catch (ImportException ime) {
//                        JOptionPane.showMessageDialog(this, "Error parsing imported file: " + ime,
//                                "Error reading file",
//                                JOptionPane.ERROR_MESSAGE);
//                        ime.printStackTrace();
//                        return;
//
//                    } catch (Exception ex) {
//                        JOptionPane.showMessageDialog(this, "Fatal exception: " + ex,
//                                "Error reading file",
//                                JOptionPane.ERROR_MESSAGE);
//                        ex.printStackTrace();
//                        return;
//                    }
//                }
//            }
//
////            setAllOptions();
//        }
    }

    public final void doImportTraits() {
        FileDialog dialog = new FileDialog(this,
                "Import Traits File...",
                FileDialog.LOAD);

        dialog.setVisible(true);
        if (dialog.getFile() != null) {
            final File file = new File(dialog.getDirectory(), dialog.getFile());

//            importTraitsFromFile(file);
        }
    }

    public void setStatusMessage() {
//        statusLabel.setText(beautiOptions.statusMessage());
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

//    public Action getImportAction() {
//        return importAlignmentAction;
//    }
//
//    protected AbstractAction importAlignmentAction = new AbstractAction("Import Alignment...") {
//        private static final long serialVersionUID = 3217702096314745005L;
//
//        public void actionPerformed(java.awt.event.ActionEvent ae) {
//            doImport();
//        }
//    };
//
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
