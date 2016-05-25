/*
 * TreeSpaceFrame.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
 *
 * This file is part of BEAST.
 * See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * BEAST is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 *  BEAST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package dr.app.treespace;

import dr.app.gui.DeleteActionResponder;
import dr.app.gui.SelectAllActionResponder;
import dr.app.treespace.inputpanel.InputPanel;
import dr.app.util.Utils;
import dr.evolution.io.Importer;
import dr.evolution.io.TreeImporter;
import dr.evolution.tree.Tree;
import jam.framework.DocumentFrame;
import jam.framework.Exportable;
import jam.util.IconUtils;
import jebl.evolution.io.ImportException;
import jebl.evolution.io.NexusImporter;
import jebl.evolution.trees.RootedTree;

import javax.swing.*;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import javax.swing.filechooser.FileNameExtensionFilter;
import javax.swing.plaf.BorderUIResource;
import java.awt.*;
import java.io.*;

/**
 * @author Andrew Rambaut
 * @version $Id$
 */
public class TreeSpaceFrame extends DocumentFrame {
    private static final long serialVersionUID = 2114148696789612509L;

    private final TreeSpaceDocument document = new TreeSpaceDocument();

    private final JTabbedPane tabbedPane = new JTabbedPane();
    private final JLabel statusLabel = new JLabel("No data loaded");
    private final InputPanel inputPanel;

    private final CladePanel cladePanel;
    private final CladePlotter cladePlotter;

    private final TreePanel treePanel;

    private final TreePlotter treePlotter;
    private final JPanel plotterPanel;

    private JFileChooser importChooser; // make JFileChooser chooser remember previous path
    private JFileChooser exportChooser; // make JFileChooser chooser remember previous path

    private final Icon gearIcon = IconUtils.getIcon(this.getClass(), "images/gear.png");

    private CladeSystem cladeSystem = new CladeSystem();

    public TreeSpaceFrame(String title) {
        super();

        setTitle(title);

        // Prevent the application to close in requestClose()
        // after a user cancel or a failure in beast file generation
        setDefaultCloseOperation(DO_NOTHING_ON_CLOSE);

//        getOpenAction().setEnabled(false);
//        getSaveAction().setEnabled(false);

        getFindAction().setEnabled(false);

        getZoomWindowAction().setEnabled(false);

        inputPanel = new InputPanel(this, document, getImportAction());
        cladePanel = new CladePanel(this, document);
        treePanel = new TreePanel(this, document);
        cladePlotter = new CladePlotter(cladeSystem);

        treePlotter = new TreePlotter();
        plotterPanel = new JPanel(new BorderLayout());
        plotterPanel.add(new JScrollPane(treePlotter), BorderLayout.CENTER);

        tabbedPane.addChangeListener(new ChangeListener() {
            public void stateChanged(final ChangeEvent e) {
                tabbedPaneChanged();
            }
        });

        document.addListener(new TreeSpaceDocument.Listener() {
            public void dataChanged() {
                setDirty();
            }

            public void settingsChanged() {
                setDirty();
            }
        });
    }

    public void initializeComponents() {

        tabbedPane.addTab("Input", inputPanel);
        tabbedPane.addTab("Trees", treePanel);
        tabbedPane.addTab("Clades", cladePanel);
        tabbedPane.addTab("Graph", cladePlotter);
        tabbedPane.addTab("Plot", plotterPanel);

        JPanel panel = new JPanel(new BorderLayout(6, 6));
        panel.setBorder(new BorderUIResource.EmptyBorderUIResource(new java.awt.Insets(12, 12, 12, 12)));
        panel.add(tabbedPane, BorderLayout.CENTER);
        panel.setPreferredSize(new java.awt.Dimension(800, 600));

//        getExportAction().setEnabled(false);
//        JButton generateButton = new JButton(getExportAction());
//        generateButton.putClientProperty("JButton.buttonType", "roundRect");

        JPanel panel2 = new JPanel(new BorderLayout(6, 6));
        panel2.add(statusLabel, BorderLayout.NORTH);
//        panel2.add(generateButton, BorderLayout.EAST);

        panel.add(panel2, BorderLayout.SOUTH);

        getContentPane().setLayout(new java.awt.BorderLayout(0, 0));
        getContentPane().add(panel, BorderLayout.CENTER);

//        setAllOptions();

        setSize(new java.awt.Dimension(1024, 768));
        setMinimumSize(new java.awt.Dimension(800, 600));

        // make JFileChooser chooser remember previous path
        exportChooser = new JFileChooser(Utils.getCWD());
        exportChooser.setDialogTitle("Generate File...");


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
//        if (isDirty()) {
//            int option = JOptionPane.showConfirmDialog(this,
//                    "You have made changes but have not generated\n" +
//                            "an output file. Do you wish to generate\n" +
//                            "before closing this window?",
//                    "Unused changes",
//                    JOptionPane.YES_NO_CANCEL_OPTION,
//                    JOptionPane.WARNING_MESSAGE);
//
//            if (option == JOptionPane.YES_OPTION) {
//                return !doGenerate();
//            } else if (option == JOptionPane.CANCEL_OPTION || option == JOptionPane.DEFAULT_OPTION) {
//                return false;
//            }
//            return true;
//        }
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

                    } catch (ImportException ime) {
                        JOptionPane.showMessageDialog(this, "Error parsing imported file: " + ime,
                                "Error reading file",
                                JOptionPane.ERROR_MESSAGE);
                        return;

//                    } catch (Exception ex) {
//                        JOptionPane.showMessageDialog(this, "Fatal exception: " + ex,
//                                "Error reading file",
//                                JOptionPane.ERROR_MESSAGE);
//                        return;
                    }
                }
            }

//            setAllOptions();
        }
    }

    public void importDataFile(File file) throws IOException, ImportException {
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
            } while (line != null && !line.contains("("));

            while (line != null) {
                if (line.contains("(")) {
                    treeCount++;
                    if (treeCount > 1) {
                        break;
                    }
                }
                line = bufferedReader.readLine();
            }

            // is a NEXUS file
            RootedTree tree = importFirstTree(file);
            if (tree != null) {
                if (treeCount > 1) {
                    InputFile inputFile = new InputFile(file, tree, -1);
                    document.addTreeFile(inputFile);

                    treePlotter.setTrees(loadTrees(inputFile));

                    inputFile.setTreeCount(treeCount);

                    document.fireDataChanged();
                } else {
                    InputFile inputFile = new InputFile(file, tree);
                    document.addTreeFile(inputFile);

                    treePlotter.setTrees(loadTrees(inputFile));

                    inputFile.setTreeCount(treeCount);

                    document.fireDataChanged();
                }
            }
        } else {
            JOptionPane.showMessageDialog(this, "Error parsing imported file. This may not be a NEXUS file",
                    "Error reading file",
                    JOptionPane.ERROR_MESSAGE);
        }
    }

    // nexus
    private RootedTree importFirstTree(File file) throws IOException, ImportException {
        RootedTree tree = null;

        FileReader reader = new FileReader(file);

        NexusImporter importer = new NexusImporter(reader);

        tree = (RootedTree)importer.importNextTree();

        return tree;
    }

    private TreeLineages loadTrees(InputFile inputFile) throws IOException {

        PrintStream progressStream = System.out;

        int totalTrees = 10000;
        int totalTreesUsed = 0;

        progressStream.println("Reading trees (bar assumes 10,000 trees)...");
        progressStream.println("0              25             50             75            100");
        progressStream.println("|--------------|--------------|--------------|--------------|");

        int stepSize = totalTrees / 60;
        if (stepSize < 1) stepSize = 1;

        TreeLineages treeLineages = new TreeLineages();

        FileReader fileReader = new FileReader(inputFile.getFile());
        jebl.evolution.io.NexusImporter importer = new NexusImporter(fileReader);
        try {
            totalTrees = 0;
            while (importer.hasTree()) {
                RootedTree tree = (RootedTree)importer.importNextTree();

                if (totalTrees >= inputFile.getBurnin()) {
                    treeLineages.addTree(tree);
                    totalTreesUsed += 1;
                }

                if (totalTrees > 0 && totalTrees % stepSize == 0) {
                    progressStream.print("*");
                    progressStream.flush();
                }
                totalTrees++;
            }

        } catch (ImportException e) {
            System.err.println("Error Parsing Input Tree: " + e.getMessage());
            return null;
        }
        fileReader.close();
        progressStream.println();
        progressStream.println();

        if (totalTrees < 1) {
            System.err.println("No trees");
            return null;
        }
        if (totalTreesUsed <= 1) {
            if (inputFile.getBurnin() > 0) {
                System.err.println("No trees to use: burnin too high");
                return null;
            }
        }

        progressStream.println("Total trees read: " + totalTrees);
        if (inputFile.getBurnin() > 0) {
            progressStream.println("Ignoring first " + inputFile.getBurnin() + " trees.");
        }

        treeLineages.setupTrees();

        return treeLineages;
    }

    private int processTrees1(InputFile inputFile) throws IOException {

        PrintStream progressStream = System.out;

        int totalTrees = 10000;
        int totalTreesUsed = 0;

        progressStream.println("Reading trees (bar assumes 10,000 trees)...");
        progressStream.println("0              25             50             75            100");
        progressStream.println("|--------------|--------------|--------------|--------------|");

        int stepSize = totalTrees / 60;
        if (stepSize < 1) stepSize = 1;

        CladeSystem cladeSystem = document.getCladeSystem();

        FileReader fileReader = new FileReader(inputFile.getFile());
        TreeImporter importer = new dr.evolution.io.NexusImporter(fileReader);
        try {
            totalTrees = 0;
            while (importer.hasTree()) {
                Tree tree = importer.importNextTree();

                if (totalTrees >= inputFile.getBurnin()) {
                    cladeSystem.add(tree, true);
                    totalTreesUsed += 1;
                }

                if (totalTrees > 0 && totalTrees % stepSize == 0) {
                    progressStream.print("*");
                    progressStream.flush();
                }
                totalTrees++;
            }

        } catch (Importer.ImportException e) {
            System.err.println("Error Parsing Input Tree: " + e.getMessage());
            return 0;
        }
        fileReader.close();
        progressStream.println();
        progressStream.println();

        if (totalTrees < 1) {
            System.err.println("No trees");
            return 0;
        }
        if (totalTreesUsed <= 1) {
            if (inputFile.getBurnin() > 0) {
                System.err.println("No trees to use: burnin too high");
                return 0;
            }
        }
        cladeSystem.normalizeClades(totalTreesUsed);

        progressStream.println("Total trees read: " + totalTrees);
        if (inputFile.getBurnin() > 0) {
            progressStream.println("Ignoring first " + inputFile.getBurnin() + " trees.");
        }

        int cladeCount = cladeSystem.getCladeMap().keySet().size();

        progressStream.println("Total unique clades: " + cladeCount);
        progressStream.println();

        progressStream.println("Processing trees for correlated clades:");
        fileReader = new FileReader(inputFile.getFile());
        importer = new dr.evolution.io.NexusImporter(fileReader);
        try {
            totalTrees = 0;
            while (importer.hasTree()) {
                Tree tree = importer.importNextTree();

                if (totalTrees >= inputFile.getBurnin()) {
                    cladeSystem.addCooccurances(tree);
                }

                if (totalTrees > 0 && totalTrees % stepSize == 0) {
                    progressStream.print("*");
                    progressStream.flush();
                }
                totalTrees++;
            }

        } catch (Importer.ImportException e) {
            System.err.println("Error Parsing Input Tree: " + e.getMessage());
            return 0;
        }
        fileReader.close();
        progressStream.println();
        progressStream.println();

        double THRESHOLD = 0.05;

        PrintWriter writer = new PrintWriter("clade_co-occurance.txt");

        writer.println("source\tsize\ttarget\tco-occurence");
        java.util.List<CladeSystem.Clade> allClades = cladeSystem.getClades();

        for (CladeSystem.Clade clade1 : allClades) {
            String name1;
            int card1 = clade1.bits.cardinality();
            if (card1 == 1) {
                name1 = clade1.label;
            } else {
                name1 = "clade" + (clade1.index + 1);
            }

            if (clade1.parents != null) {
                for (CladeSystem.Clade clade2 : clade1.parents.keySet()) {
                    String name2;
                    int card2 = clade2.bits.cardinality();
                    name2 = "clade" + (clade2.index+1);

                    double value = clade1.parents.get(clade2);
                    value /= totalTreesUsed;
                    if (value > THRESHOLD)  {
                        if (card1 > card2) {
                            writer.println(name1 + "_" + card1 + "\t" + card1 + "\t" + name2 + "_" + card2 + "\t" + value);
                        } else {
                            writer.println(name2 + "_" + card2 + "\t" + card2 + "\t" + name1 + "_" + card1 + "\t" + value);
                        }
                    }
                }
            }
        }
        writer.close();

        writer = new PrintWriter("clade_frequency.txt");

        writer.println("source\tsize\tfrequency");

        for (CladeSystem.Clade clade1 : allClades) {
            String name1;
            int card1 = clade1.bits.cardinality();
            if (card1 == 1) {
                name1 = clade1.label;
            } else {
                name1 = "clade" + (clade1.index + 1);
            }

            double value = clade1.count;
            value /= totalTreesUsed;
            if (value > THRESHOLD)  {
                writer.println(name1 + "_" + card1 + "\t" + card1 + "\t" + value);
            }
        }
        writer.close();
        progressStream.println();

        cladePlotter.setCladeSystem(cladeSystem);

        return totalTreesUsed;
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
