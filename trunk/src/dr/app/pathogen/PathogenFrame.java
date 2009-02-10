/*
 * BeautiFrame.java
 *
 * (c) 2002-2005 BEAST Development Core Team
 *
 * This package may be distributed under the
 * Lesser Gnu Public Licence (LGPL)
 */
package dr.app.pathogen;

import dr.evolution.alignment.SimpleAlignment;
import dr.evolution.alignment.Alignment;
import dr.evolution.io.Importer;
import dr.evolution.io.NexusImporter;
import dr.evolution.tree.Tree;
import dr.evolution.util.TaxonList;
import dr.evolution.util.Units;
import org.virion.jam.framework.DocumentFrame;
import org.virion.jam.framework.Exportable;
import org.virion.jam.util.IconUtils;

import javax.swing.*;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import javax.swing.plaf.BorderUIResource;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;

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

        getFindAction().setEnabled(false);

        getZoomWindowAction().setEnabled(false);
    }

    public void initializeComponents() {

        samplesPanel = new SamplesPanel(this, taxa);
        treesPanel = new TreesPanel(this, trees);

        tabbedPane.addTab("Sample Dates", samplesPanel);
        tabbedPane.addTab("Trees", treesPanel);

        JPanel panel = new JPanel(new BorderLayout(6, 6));
        panel.setBorder(new BorderUIResource.EmptyBorderUIResource(new Insets(12, 12, 12, 12)));
        panel.add(tabbedPane, BorderLayout.CENTER);

        panel.add(statusLabel, BorderLayout.SOUTH);

        getContentPane().setLayout(new BorderLayout(0, 0));
        getContentPane().add(panel, BorderLayout.CENTER);

        setSize(new Dimension(1024, 768));
    }

    public void timeScaleChanged() {
        treesPanel.timeScaleChanged();
    }

    protected boolean readFromFile(File file) throws IOException {

        try {
            FileReader reader = new FileReader(file);

            NexusImporter importer = new NexusImporter(reader);

            boolean done = false;

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
                        trees.addAll(Arrays.asList(treeArray));

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

        int count = 0;
        if (trees.size() > 0) {
            for (Tree tree : trees) {
                String id = tree.getId();
                if (id == null || id.trim().length() == 0) {
                    count++;
                    tree.setId("tree_" + (count + 1));
                }
            }
        }

        setStatusMessage();

        getExportAction().setEnabled(true);

        return true;
    }

    protected boolean writeToFile(File file) throws IOException {
        return false;  //To change body of implemented methods use File | Settings | File Templates.
    }

    private void setStatusMessage() {
//        if (trees.size() > 0) {
//            message += "Trees only : " + beautiOptions.trees.size() +
//                    (beautiOptions.trees.size() > 1 ? " trees, " : " tree, ") +
//                    beautiOptions.taxonList.getTaxonCount() + " taxa";
//        } else {
//            message += "Taxa only: " + beautiOptions.taxonList.getTaxonCount() + " taxa";
//        }
//        statusLabel.setText(message);
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

    public Action getExportAction() {
        return generateAction;
    }

    protected AbstractAction generateAction = new AbstractAction("Export...") {
        public void actionPerformed(ActionEvent ae) {
//            doGenerate();
        }
    };

}