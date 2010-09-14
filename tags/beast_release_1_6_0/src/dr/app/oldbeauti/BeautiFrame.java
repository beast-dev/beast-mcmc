/*
 * BeautiFrame.java
 *
 * (c) 2002-2005 BEAST Development Core Team
 *
 * This package may be distributed under the
 * Lesser Gnu Public Licence (LGPL)
 */
package dr.app.oldbeauti;

import dr.evolution.alignment.Patterns;
import dr.evolution.alignment.SimpleAlignment;
import dr.evolution.distance.DistanceMatrix;
import dr.evolution.distance.JukesCantorDistanceMatrix;
import dr.evolution.io.Importer;
import dr.evolution.io.NexusImporter;
import dr.evolution.tree.Tree;
import dr.evolution.util.Units;
import org.jdom.Document;
import org.jdom.JDOMException;
import org.jdom.input.SAXBuilder;
import org.jdom.output.XMLOutputter;
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

/**
 * @author			Andrew Rambaut
 * @author			Alexei Drummond
 * @version			$Id: BeautiFrame.java,v 1.22 2006/09/09 16:07:06 rambaut Exp $
 */
public class BeautiFrame extends DocumentFrame {

    /**
     *
     */
    private static final long serialVersionUID = 2114148696789612509L;

    private final BeastGenerator beautiOptions = new BeastGenerator();

    private final JTabbedPane tabbedPane = new JTabbedPane();
    private final JLabel statusLabel = new JLabel("No data loaded");

    private DataPanel dataPanel;
    private TaxaPanel taxaPanel;
    private ModelPanel modelPanel;
    private PriorsPanel priorsPanel;
    private OperatorsPanel operatorsPanel;
    private MCMCPanel mcmcPanel;

//	final Icon dataIcon = new ImageIcon(Utils.getImage(this, "/images/data-icon.gif"));
//	final Icon modelIcon = new ImageIcon(Utils.getImage(this, "/images/model-icon.gif"));
//	final Icon mcmcIcon = new ImageIcon(Utils.getImage(this, "/images/mcmc-icon.gif"));

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
    }

    public void initializeComponents() {

        dataPanel = new DataPanel(this);
        taxaPanel = new TaxaPanel(this);
        modelPanel = new ModelPanel(this);
        priorsPanel = new PriorsPanel(this);
        operatorsPanel = new OperatorsPanel(this);
        mcmcPanel = new MCMCPanel(this);

        tabbedPane.addTab("Data", dataPanel);
        tabbedPane.addTab("Taxa", taxaPanel);
        tabbedPane.addTab("Model", modelPanel);
        tabbedPane.addTab("Priors", priorsPanel);
        tabbedPane.addTab("Operators", operatorsPanel);
        tabbedPane.addTab("MCMC", mcmcPanel);
        tabbedPane.addChangeListener(new ChangeListener() {
            public void stateChanged(ChangeEvent e) {
                if (tabbedPane.getSelectedComponent() == dataPanel) {
                    dataPanel.selectionChanged();
                } else {
                    getDeleteAction().setEnabled(false);
                }
            }
        });

        JPanel panel = new JPanel(new BorderLayout(6, 6));
        panel.setBorder(new BorderUIResource.EmptyBorderUIResource(new java.awt.Insets(12, 12, 12, 12)));
        panel.add(tabbedPane, BorderLayout.CENTER);

        getExportAction().setEnabled(false);
        JButton generateButton = new JButton(getExportAction());
	    generateButton.putClientProperty("JButton.buttonType", "roundRect");

        JPanel panel2 = new JPanel(new BorderLayout(6, 6));
        panel2.add(statusLabel, BorderLayout.CENTER);
        panel2.add(generateButton, BorderLayout.EAST);

        panel.add(panel2, BorderLayout.SOUTH);

        getContentPane().setLayout(new java.awt.BorderLayout(0, 0));
        getContentPane().add(panel, BorderLayout.CENTER);

        dataPanel.setOptions(beautiOptions);
        taxaPanel.setOptions(beautiOptions);
        modelPanel.setOptions(beautiOptions);
        priorsPanel.setOptions(beautiOptions);
        operatorsPanel.setOptions(beautiOptions);
        mcmcPanel.setOptions(beautiOptions);

        setSize(new java.awt.Dimension(800, 600));
    }

    public final void dataChanged() {
        taxaPanel.setOptions(beautiOptions);
        modelPanel.setOptions(beautiOptions);
        priorsPanel.setOptions(beautiOptions);
        operatorsPanel.setOptions(beautiOptions);
        setDirty();
    }

    public final void dataSelectionChanged(boolean isSelected) {
        if (isSelected) {
            getDeleteAction().setEnabled(true);
        } else {
            getDeleteAction().setEnabled(false);
        }
    }

    public void taxonSetsChanged() {
        priorsPanel.setOptions(beautiOptions);
        setDirty();
    }


    public void doDelete() {
        if (tabbedPane.getSelectedComponent() == dataPanel) {
            dataPanel.deleteSelection();
        } else {
            throw new RuntimeException("Delete should only be accessable from the Data panel");
        }
    }

    public final void modelChanged() {
        modelPanel.getOptions(beautiOptions);

        priorsPanel.setOptions(beautiOptions);
        operatorsPanel.setOptions(beautiOptions);
        setDirty();
    }

    public final void operatorsChanged() {
        setDirty();
    }

    public final void priorsChanged() {
        priorsPanel.getOptions(beautiOptions);

        operatorsPanel.setOptions(beautiOptions);

        priorsPanel.setOptions(beautiOptions);

        setDirty();
    }

    public final void mcmcChanged() {
        setDirty();
    }

    public boolean requestClose() {
        if (isDirty()) {
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
                JOptionPane.showMessageDialog(this, "Unable to read template file: " + ioe,
                        "Unable to read file",
                        JOptionPane.ERROR_MESSAGE);
            }
        }
    }

    protected boolean readFromFile(File file) throws IOException {
        try {
            SAXBuilder parser = new SAXBuilder();
            Document doc = parser.build(file);
            beautiOptions.parse(doc);

            if (beautiOptions.guessDates) {
                beautiOptions.guessDates();
            }

            dataPanel.setOptions(beautiOptions);
            taxaPanel.setOptions(beautiOptions);
            modelPanel.setOptions(beautiOptions);
            priorsPanel.setOptions(beautiOptions);
            operatorsPanel.setOptions(beautiOptions);
            mcmcPanel.setOptions(beautiOptions);

            getExportAction().setEnabled(beautiOptions.alignment != null);
            getSaveAction().setEnabled(beautiOptions.alignment != null);
            getSaveAsAction().setEnabled(beautiOptions.alignment != null);

        } catch (dr.xml.XMLParseException xpe) {
            JOptionPane.showMessageDialog(this, "Error reading file: This may not be a BEAUti file",
                    "Error reading file",
                    JOptionPane.ERROR_MESSAGE);
            return false;
        } catch (JDOMException e) {
            JOptionPane.showMessageDialog(this, "Unable to open file: This may not be a BEAUti file",
                    "Unable to open file",
                    JOptionPane.ERROR_MESSAGE);
            return false;
        }
        return true;
    }

    public String getDefaultFileName() { return beautiOptions.fileNameStem+".beauti"; }

    protected boolean writeToFile(File file) throws IOException {
        dataPanel.getOptions(beautiOptions);
        taxaPanel.getOptions(beautiOptions);
        modelPanel.getOptions(beautiOptions);
        priorsPanel.getOptions(beautiOptions);
        operatorsPanel.getOptions(beautiOptions);
        mcmcPanel.getOptions(beautiOptions);

        Document doc = beautiOptions.create(false, true);

        FileWriter fw = new FileWriter(file);

        XMLOutputter outputter = new XMLOutputter(org.jdom.output.Format.getPrettyFormat());

        outputter.output(doc, fw);

        fw.close();
        return true;
    }

    public final void doImport() {

        FileDialog dialog = new FileDialog(this,
                "Import NEXUS File...",
                FileDialog.LOAD);

        dialog.setVisible(true);
        if (dialog.getFile() != null) {
            File file = new File(dialog.getDirectory(), dialog.getFile());

            try {
                importFromFile(file);

                setDirty();
            } catch (FileNotFoundException fnfe) {
                JOptionPane.showMessageDialog(this, "Unable to open file: File not found",
                        "Unable to open file",
                        JOptionPane.ERROR_MESSAGE);
            } catch (IOException ioe) {
                JOptionPane.showMessageDialog(this, "Unable to read file: " + ioe,
                        "Unable to read file",
                        JOptionPane.ERROR_MESSAGE);
            }
        }

    }

    protected void importFromFile(File file) throws IOException {

        try {
            FileReader reader = new FileReader(file);

            NexusApplicationImporter importer = new NexusApplicationImporter(reader);

            boolean done = false;

            beautiOptions.originalAlignment = null;
            beautiOptions.alignment = null;
            beautiOptions.tree = null;
            beautiOptions.taxonList = null;

            while (!done) {
                try {

                    NexusImporter.NexusBlock block = importer.findNextBlock();

                    if (block == NexusImporter.TAXA_BLOCK) {

                        if (beautiOptions.taxonList != null) {
                            throw new NexusImporter.MissingBlockException("TAXA block already defined");
                        }

                        beautiOptions.taxonList = importer.parseTaxaBlock();

                    } else if (block == NexusImporter.CALIBRATION_BLOCK) {
                        if (beautiOptions.taxonList == null) {
                            throw new NexusImporter.MissingBlockException("TAXA or DATA block must be defined before a CALIBRATION block");
                        }

                        importer.parseCalibrationBlock(beautiOptions.taxonList);

                    } else if (block == NexusImporter.CHARACTERS_BLOCK) {

                        if (beautiOptions.taxonList == null) {
                            throw new NexusImporter.MissingBlockException("TAXA block must be defined before a CHARACTERS block");
                        }

                        if (beautiOptions.originalAlignment != null) {
                            throw new NexusImporter.MissingBlockException("CHARACTERS or DATA block already defined");
                        }

                        beautiOptions.originalAlignment = (SimpleAlignment)importer.parseCharactersBlock(beautiOptions.taxonList);

                    } else if (block == NexusImporter.DATA_BLOCK) {

                        if (beautiOptions.originalAlignment != null) {
                            throw new NexusImporter.MissingBlockException("CHARACTERS or DATA block already defined");
                        }

                        // A data block doesn't need a taxon block before it
                        // but if one exists then it will use it.
                        beautiOptions.originalAlignment = (SimpleAlignment)importer.parseDataBlock(beautiOptions.taxonList);
                        if (beautiOptions.taxonList == null) {
                            beautiOptions.taxonList = beautiOptions.originalAlignment;
                        }

                    } else if (block == NexusImporter.TREES_BLOCK) {

                        if (beautiOptions.taxonList == null) {
                            throw new NexusImporter.MissingBlockException("TAXA or DATA block must be defined before a TREES block");
                        }

                        if (beautiOptions.tree != null) {
                            throw new NexusImporter.MissingBlockException("TREES block already defined");
                        }

                        Tree[] trees = importer.parseTreesBlock(beautiOptions.taxonList);
                        if (trees.length > 0) {
                            beautiOptions.tree = trees[0];
                        }

/*					} else if (block == NexusApplicationImporter.PAUP_BLOCK) {

						importer.parsePAUPBlock(beautiOptions);

					} else if (block == NexusApplicationImporter.MRBAYES_BLOCK) {

						importer.parseMrBayesBlock(beautiOptions);

					} else if (block == NexusApplicationImporter.RHINO_BLOCK) {

						importer.parseRhinoBlock(beautiOptions);
*/
                    } else {
                        // Ignore the block..
                    }

                } catch (EOFException ex) {
                    done = true;
                }
            }

            // Allow the user to load taxa only (perhaps from a tree file) so that they can sample from a prior...
            if (beautiOptions.originalAlignment == null && beautiOptions.taxonList == null) {
                throw new NexusImporter.MissingBlockException("TAXON, DATA or CHARACTERS block is missing");
            }

        } catch (Importer.ImportException ime) {
            JOptionPane.showMessageDialog(this, "Error parsing imported file: " + ime,
                    "Error reading file",
                    JOptionPane.ERROR_MESSAGE);
            return;
        } catch (IOException ioex) {
            JOptionPane.showMessageDialog(this, "File I/O Error: " + ioex,
                    "File I/O Error",
                    JOptionPane.ERROR_MESSAGE);
            return;
        } catch (Exception ex) {
            JOptionPane.showMessageDialog(this, "Fatal exception: " + ex,
                    "Error reading file",
                    JOptionPane.ERROR_MESSAGE);
            return;
        }

        // check the taxon names for invalid characters
        boolean foundAmp = false;
        for (int i = 0; i < beautiOptions.taxonList.getTaxonCount(); i++) {
            String name = beautiOptions.taxonList.getTaxon(i).getId();
            if (name.indexOf('&') >= 0) {
                foundAmp = true;
            }
        }
        if (foundAmp) {
            JOptionPane.showMessageDialog(this, "One or more taxon names include an illegal character ('&').\n" +
                    "These characters will prevent BEAST from reading the resulting XML file.\n\n" +
                    "Please edit the taxon name(s) before generating the BEAST file.",
                    "Illegal Taxon Name(s)",
                    JOptionPane.WARNING_MESSAGE);
        }


        // make sure they all have dates...
        for (int i = 0; i < beautiOptions.taxonList.getTaxonCount(); i++) {
            if (beautiOptions.taxonList.getTaxonAttribute(i, "date") == null) {
                java.util.Date origin = new java.util.Date(0);

                dr.evolution.util.Date date = dr.evolution.util.Date.createTimeSinceOrigin(0.0, Units.Type.YEARS, origin);
                beautiOptions.taxonList.getTaxon(i).setAttribute("date", date);
            }
        }

        beautiOptions.fileNameStem = dr.app.util.Utils.trimExtensions(file.getName(),
                new String[] {"nex", "NEX", "tre", "TRE", "nexus", "NEXUS"});

        beautiOptions.alignment = beautiOptions.originalAlignment;
        beautiOptions.alignmentReset = true;
        if (beautiOptions.alignment != null) {
            Patterns patterns = new Patterns(beautiOptions.alignment);
            DistanceMatrix distances = new JukesCantorDistanceMatrix(patterns);
            beautiOptions.meanDistance = distances.getMeanDistance();

            statusLabel.setText("Alignment: " + beautiOptions.alignment.getTaxonCount() + " taxa, " +
                    beautiOptions.alignment.getSiteCount() + " sites");
            beautiOptions.dataType = beautiOptions.alignment.getDataType().getType();
        } else {
            statusLabel.setText("Taxa only: " + beautiOptions.taxonList.getTaxonCount() + " taxa");
            beautiOptions.meanDistance = 0.0;
        }

        dataPanel.setOptions(beautiOptions);
        taxaPanel.setOptions(beautiOptions);
        modelPanel.setOptions(beautiOptions);
        priorsPanel.setOptions(beautiOptions);
        operatorsPanel.setOptions(beautiOptions);
        mcmcPanel.setOptions(beautiOptions);


        getOpenAction().setEnabled(true);
        getSaveAction().setEnabled(true);
        getExportAction().setEnabled(true);
    }

    public final boolean doGenerate() {

        try {
            beautiOptions.checkOptions();
        } catch(IllegalArgumentException iae) {
            JOptionPane.showMessageDialog(this, iae.getMessage(),
                    "Unable to generate file",
                    JOptionPane.ERROR_MESSAGE);
            return false;
        }

        FileDialog dialog = new FileDialog(this,
                "Generate BEAST File...",
                FileDialog.SAVE);

        dialog.setVisible(true);
        dialog.setFile(beautiOptions.fileNameStem + ".xml");
        if (dialog.getFile() != null) {
            File file = new File(dialog.getDirectory(), dialog.getFile());

            try {
                generate(file);

            } catch (IOException ioe) {
                JOptionPane.showMessageDialog(this, "Unable to generate file: " + ioe.getMessage(),
                        "Unable to generate file",
                        JOptionPane.ERROR_MESSAGE);
                return false;
            }
        }

        clearDirty();
        return true;
    }

    protected void generate(File file) throws IOException {
        dataPanel.getOptions(beautiOptions);
        taxaPanel.getOptions(beautiOptions);
        modelPanel.getOptions(beautiOptions);
        priorsPanel.getOptions(beautiOptions);
        operatorsPanel.getOptions(beautiOptions);
        mcmcPanel.getOptions(beautiOptions);

        FileWriter fw = new FileWriter(file);
        beautiOptions.generateXML(fw);
        fw.close();
    }

    public JComponent getExportableComponent() {

        JComponent exportable = null;
        Component comp = tabbedPane.getSelectedComponent();

        if (comp instanceof Exportable) {
            exportable = ((Exportable)comp).getExportableComponent();
        } else if (comp instanceof JComponent) {
            exportable = (JComponent)comp;
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

    public Action getImportAction() { return importNexusAction; }

    protected AbstractAction importNexusAction = new AbstractAction("Import NEXUS...") {
        private static final long serialVersionUID = 3217702096314745005L;

        public void actionPerformed(java.awt.event.ActionEvent ae) {
            doImport();
        }
    };

    public Action getExportAction() { return generateAction; }

    protected AbstractAction generateAction = new AbstractAction("Generate BEAST File...", gearIcon) {
        private static final long serialVersionUID = -5329102618630268783L;

        public void actionPerformed(java.awt.event.ActionEvent ae) {
            doGenerate();
        }
    };

}
