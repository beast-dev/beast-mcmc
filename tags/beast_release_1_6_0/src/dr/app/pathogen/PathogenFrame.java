/*
 * BeautiFrame.java
 *
 * (c) 2002-2005 BEAST Development Core Team
 *
 * This package may be distributed under the
 * Lesser Gnu Public Licence (LGPL)
 */
package dr.app.pathogen;

import dr.evolution.io.Importer;
import dr.evolution.io.NexusImporter;
import dr.evolution.tree.Tree;
import dr.evolution.tree.FlexibleTree;
import dr.evolution.util.TaxonList;
import dr.app.tools.NexusExporter;
import dr.app.pathogen.TemporalRooting;
import dr.util.NumberFormatter;
import jam.framework.DocumentFrame;
import jam.framework.Exportable;

import javax.swing.*;
import javax.swing.plaf.BorderUIResource;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.io.*;
import java.util.ArrayList;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: BeautiFrame.java,v 1.22 2006/09/09 16:07:06 rambaut Exp $
 */
public class PathogenFrame extends DocumentFrame {

    /**
     *
     */
    private static final long serialVersionUID = 2114148696789612509L;

    private JTabbedPane tabbedPane = new JTabbedPane();
    private JLabel statusLabel = new JLabel("No data loaded");

    private SamplesPanel samplesPanel;
    private TreesPanel treesPanel;

    TaxonList taxa = null;
    java.util.List<Tree> trees = new ArrayList<Tree>();

    public PathogenFrame(String title) {
        super();

        setTitle(title);

        // Prevent the application to close in requestClose()
        // after a user cancel or a failure in beast file generation
        setDefaultCloseOperation(DO_NOTHING_ON_CLOSE);

        getOpenAction().setEnabled(true);
        getSaveAction().setEnabled(false);
        getSaveAsAction().setEnabled(false);

        getFindAction().setEnabled(false);

        getZoomWindowAction().setEnabled(false);
    }

    public void initializeComponents() {

        samplesPanel = new SamplesPanel(this, taxa);
        treesPanel = new TreesPanel(this, trees.get(0));

        tabbedPane.addTab("Sample Dates", samplesPanel);
        tabbedPane.addTab("Analysis", treesPanel);

        JPanel panel = new JPanel(new BorderLayout(6, 6));
        panel.setBorder(new BorderUIResource.EmptyBorderUIResource(new Insets(12, 12, 12, 12)));
        panel.add(tabbedPane, BorderLayout.CENTER);

        panel.add(statusLabel, BorderLayout.SOUTH);

        getContentPane().setLayout(new BorderLayout(0, 0));
        getContentPane().add(panel, BorderLayout.CENTER);

        setSize(new Dimension(1024, 768));

        setStatusMessage();
    }

    public void timeScaleChanged() {
        treesPanel.timeScaleChanged();
        setStatusMessage();
    }

    protected boolean readFromFile(File file) throws IOException {

        try {
            FileReader reader = new FileReader(file);

            NexusImporter importer = new NexusImporter(reader);

            boolean done = false;

            int count = 0;

            while (!done) {
                try {

                    NexusImporter.NexusBlock block = importer.findNextBlock();

                    if (block == NexusImporter.TAXA_BLOCK) {

                        if (taxa != null) {
                            throw new NexusImporter.MissingBlockException("TAXA block already defined");
                        }

                        taxa = importer.parseTaxaBlock();

                    } else if (block == NexusImporter.TREES_BLOCK) {

                        Tree[] treeArray = importer.parseTreesBlock(taxa);
                        for (Tree tree : treeArray) {
                            // get a rooted version of the tree to clone
                            FlexibleTree binaryTree = new FlexibleTree(tree, true);
                            binaryTree.resolveTree();

                            String id = binaryTree.getId();
                            if (id == null || id.trim().length() == 0) {
                                count++;
                                binaryTree.setId("tree_" + (count + 1));
                            }
                            trees.add(binaryTree);
                        }

                        if (taxa == null && trees.size() > 0) {
                            taxa = trees.get(0);
                        }


                    } else {
                        // Ignore the block..
                    }

                } catch (EOFException ex) {
                    done = true;
                }
            }

        } catch (Importer.ImportException ime) {
            JOptionPane.showMessageDialog(this, "Error parsing imported file: " + ime,
                    "Error reading file",
                    JOptionPane.ERROR_MESSAGE);
            ime.printStackTrace();
            return false;
        } catch (IOException ioex) {
            JOptionPane.showMessageDialog(this, "File I/O Error: " + ioex,
                    "File I/O Error",
                    JOptionPane.ERROR_MESSAGE);
            ioex.printStackTrace();
            return false;
        } catch (Exception ex) {
            JOptionPane.showMessageDialog(this, "Fatal exception: " + ex,
                    "Error reading file",
                    JOptionPane.ERROR_MESSAGE);
            ex.printStackTrace();
            return false;
        }

        getExportTreeAction().setEnabled(true);
        getExportDataAction().setEnabled(true);

        return true;
    }

    protected boolean writeToFile(File file) throws IOException {
        return false;  //To change body of implemented methods use File | Settings | File Templates.
    }

    protected void doExportTree() {
        FileDialog dialog = new FileDialog(this,
                "Export Tree File...",
                FileDialog.SAVE);

        dialog.setVisible(true);
        if (dialog.getFile() != null) {
            File file = new File(dialog.getDirectory(), dialog.getFile());

            PrintStream ps = null;
            try {
                ps = new PrintStream(file);
                writeTreeFile(ps, false);
                ps.close();
            } catch (IOException ioe) {
                JOptionPane.showMessageDialog(this, "Error writing tree file: " + ioe.getMessage(),
                        "Export Error",
                        JOptionPane.ERROR_MESSAGE);
            }

        }
    }

    protected void writeTreeFile(PrintStream ps, boolean newickFormat) throws IOException {

        Tree tree = treesPanel.getTreeAsViewed();

//        if (newickFormat) {
//            NewickExporter newickExporter = new NewickExporter(ps);
//            newickExporter.exportTree(tree);
//        } else {
            NexusExporter nexusExporter = new NexusExporter(new PrintStream(ps));
            nexusExporter.exportTree(tree);
//        }
    }

//    protected void doExportGraphic() {
//        ExportDialog export = new ExportDialog();
//        export.showExportDialog( this, "Export view as ...", treeViewer.getContentPane(), "export" );
//    }

    protected void doExportData() {
        FileDialog dialog = new FileDialog(this,
                "Export Data File...",
                FileDialog.SAVE);

        dialog.setVisible(true);
        if (dialog.getFile() != null) {
            File file = new File(dialog.getDirectory(), dialog.getFile());

            Writer writer = null;
            try {
                writer = new PrintWriter(file);
                treesPanel.writeDataFile(writer);
                writer.close();
            } catch (IOException ioe) {
                JOptionPane.showMessageDialog(this, "Error writing data file: " + ioe.getMessage(),
                        "Export Error",
                        JOptionPane.ERROR_MESSAGE);
            }

        }
    }

    private void setStatusMessage() {
        Tree tree = treesPanel.getTree();
        if (tree != null) {
            String message = "";
            message += "Tree loaded, " + tree.getTaxonCount() + " taxa";

            TemporalRooting tr = treesPanel.getTemporalRooting();
            if (tr.isContemporaneous()) {
                message += ", contemporaneous tips";
            } else {
                NumberFormatter nf = new NumberFormatter(3);
                message += ", dated tips with range " + nf.format(tr.getDateRange());
            }
            statusLabel.setText(message);
        }
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

    public Action getExportTreeAction() {
        return exportTreeAction;
    }

//    public Action getExportGraphicAction() {
//        return exportGraphicAction;
//    }

    public Action getExportDataAction() {
        return exportDataAction;
    }

    protected AbstractAction exportTreeAction = new AbstractAction("Export Tree...") {
        public void actionPerformed(ActionEvent ae) {
            doExportTree();
        }
    };

//    protected AbstractAction exportGraphicAction = new AbstractAction("Export Graphic...") {
//        public void actionPerformed(ActionEvent ae) {
//            doExportGraphic();
//        }
//    };

    protected AbstractAction exportDataAction = new AbstractAction("Export Data...") {
        public void actionPerformed(ActionEvent ae) {
            doExportData();
        }
    };
}