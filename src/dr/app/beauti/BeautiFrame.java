/*
 * BeautiFrame.java
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

/*
 * BeautiFrame.java
 *
 * (c) 2002-2005 BEAST Development Core Team
 *
 * This package may be distributed under the
 * Lesser Gnu Public Licence (LGPL)
 */
package dr.app.beauti;

import dr.app.beauti.ancestralstatespanel.AncestralStatesPanel;
import dr.app.beauti.clockmodelspanel.ClockModelsPanel;
import dr.app.beauti.components.ComponentFactory;
import dr.app.beauti.components.ancestralstates.AncestralStatesComponentFactory;
import dr.app.beauti.components.continuous.ContinuousComponentFactory;
import dr.app.beauti.components.discrete.DiscreteTraitsComponentFactory;
import dr.app.beauti.components.dollo.DolloComponentFactory;
import dr.app.beauti.components.hpm.HierarchicalModelComponentFactory;
import dr.app.beauti.components.linkedparameters.LinkedParameterComponentFactory;
import dr.app.beauti.components.marginalLikelihoodEstimation.MarginalLikelihoodEstimationComponentFactory;
import dr.app.beauti.components.sequenceerror.SequenceErrorModelComponentFactory;
import dr.app.beauti.components.tipdatesampling.TipDateSamplingComponentFactory;
import dr.app.beauti.datapanel.DataPanel;
import dr.app.beauti.generator.BeastGenerator;
import dr.app.beauti.generator.Generator;
import dr.app.beauti.mcmcpanel.MCMCPanel;
import dr.app.beauti.operatorspanel.OperatorsPanel;
import dr.app.beauti.options.BeautiOptions;
import dr.app.beauti.options.PartitionTreePrior;
import dr.app.beauti.options.STARBEASTOptions;
import dr.app.beauti.options.TraitData;
import dr.app.beauti.priorspanel.DefaultPriorTableDialog;
import dr.app.beauti.priorspanel.PriorsPanel;
import dr.app.beauti.sitemodelspanel.SiteModelsPanel;
import dr.app.beauti.taxonsetspanel.SpeciesSetPanel;
import dr.app.beauti.taxonsetspanel.TaxonSetPanel;
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
import org.jdom.JDOMException;

import javax.swing.*;
import javax.swing.border.Border;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import javax.swing.filechooser.FileNameExtensionFilter;
import javax.swing.plaf.BorderUIResource;
import java.awt.*;
import java.awt.event.HierarchyBoundsListener;
import java.awt.event.HierarchyEvent;
import java.io.*;
import java.util.HashMap;
import java.util.Map;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: BeautiFrame.java,v 1.22 2006/09/09 16:07:06 rambaut Exp $
 */
public class BeautiFrame extends DocumentFrame {

    private static final long serialVersionUID = 2114148696789612509L;

    public final static String DATA_PARTITIONS = "Partitions";
    public final static String TAXON_SETS = "Taxa";
    public final static String TIP_DATES = "Tips";
    public final static String TRAITS = "Traits";
    public final static String SITE_MODELS = "Sites";
    public final static String CLOCK_MODELS = "Clocks";
    public final static String TREES = "Trees";
    public final static String ANCESTRAL_STATES = "States";
    public final static String PRIORS = "Priors";
    public final static String OPERATORS = "Operators";
    public final static String MCMC = "MCMC";

    private BeautiOptions options;
    private BeastGenerator generator;
    private final ComponentFactory[] components;

    public final JTabbedPane tabbedPane = new JTabbedPane();
    public final JLabel statusLabel = new JLabel();

    private DataPanel dataPanel;
    private TipDatesPanel tipDatesPanel;
    private TraitsPanel traitsPanel;
    private TaxonSetPanel taxonSetPanel;
    private SpeciesSetPanel speciesSetPanel;
    private SiteModelsPanel siteModelsPanel;
    private AncestralStatesPanel ancestralStatesPanel;
    private ClockModelsPanel clockModelsPanel;
    private TreesPanel treesPanel;
    private PriorsPanel priorsPanel;
    private OperatorsPanel operatorsPanel;
    private MCMCPanel mcmcPanel;

    private BeautiPanel currentPanel;

    private Map<String, FileDialog> fileDialogs = new HashMap<String, FileDialog>();
    private Map<String, JFileChooser> fileChoosers = new HashMap<String, JFileChooser>();

    final Icon gearIcon = IconUtils.getIcon(this.getClass(), "images/gear.png");

    public BeautiFrame(String title) {
        super();

        setTitle(title);

        // Prevent the application to close in requestClose()
        // after a user cancel or a failure in beast file generation
        setDefaultCloseOperation(DO_NOTHING_ON_CLOSE);

//        getOpenAction().setEnabled(false);
//        getSaveAction().setEnabled(false);

        getFindAction().setEnabled(false);

        getZoomWindowAction().setEnabled(false);

        components = new ComponentFactory[] {
                AncestralStatesComponentFactory.INSTANCE,
                ContinuousComponentFactory.INSTANCE,
                DiscreteTraitsComponentFactory.INSTANCE,
//                DnDsComponentFactory.INSTANCE,
                DolloComponentFactory.INSTANCE,
                LinkedParameterComponentFactory.INSTANCE,
                HierarchicalModelComponentFactory.INSTANCE,
                MarginalLikelihoodEstimationComponentFactory.INSTANCE,
                SequenceErrorModelComponentFactory.INSTANCE,
                TipDateSamplingComponentFactory.INSTANCE
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
        traitsPanel = new TraitsPanel(this, dataPanel, getImportTraitsAction());
        taxonSetPanel = new TaxonSetPanel(this);
        speciesSetPanel = new SpeciesSetPanel(this);
        siteModelsPanel = new SiteModelsPanel(this, getDeleteAction());
        ancestralStatesPanel = new AncestralStatesPanel(this);
        clockModelsPanel = new ClockModelsPanel(this);
//        clockmodelspanel = new OldClockModelsPanel(this);
//        oldTreesPanel = new OldTreesPanel(this);
        treesPanel = new TreesPanel(this, getDeleteAction());
//        speciesTreesPanel = new SpeciesTreesPanel(this);
        priorsPanel = new PriorsPanel(this, false);
        operatorsPanel = new OperatorsPanel(this);
        mcmcPanel = new MCMCPanel(this);

        int index = 0;
        tabbedPane.addTab(DATA_PARTITIONS, dataPanel);
        tabbedPane.setToolTipTextAt(index++, "<html>" +
                "Import sequence alignments, organize data partitions,<br>" +
                "link models between partitions and select *BEAST</html>");

        tabbedPane.addTab(TAXON_SETS, taxonSetPanel);
        tabbedPane.setToolTipTextAt(index++, "<html>" +
                "Create and edit sets of taxa which can be used to <br>" +
                "define times of most recent common ancestors and <br>" +
                "to keep groups monophyletic.</html>");
//        tabbedPane.addTab("Species Sets", speciesSetPanel);
        tabbedPane.addTab(TIP_DATES, tipDatesPanel);
        tabbedPane.setToolTipTextAt(index++, "<html>" +
                "Specify sampling dates of tips for use in temporal <br>" +
                "analyses of measurably evolving populations.</html>");
        tabbedPane.addTab(TRAITS, traitsPanel);
        tabbedPane.setToolTipTextAt(index++, "<html>" +
                "Import and organize continuous and discrete traits <br>" +
                "for taxa, convert them into data partitions for evolutionary<br>" +
                "analysis.</html>");
        tabbedPane.addTab(SITE_MODELS, siteModelsPanel);
        tabbedPane.setToolTipTextAt(index++, "<html>" +
                "Select evolutionary models to be used for each data <br>" +
                "partition including substitution models, codon partitioning<br>" +
                "and trait evolution models.</html>");
        tabbedPane.addTab(CLOCK_MODELS, clockModelsPanel);
        tabbedPane.setToolTipTextAt(index++, "<html>" +
                "Select relaxed molecular clock models to be used across <br>" +
                "the tree. Specify sampling of rates.</html>");
        tabbedPane.addTab(TREES, treesPanel);
        tabbedPane.setToolTipTextAt(index++, "<html>" +
                "Select the priors on trees including coalescent models<br>" +
                "birth-death speciation models and the *BEAST gene tree,<br>" +
                "species tree options.</html>");

        tabbedPane.addTab(ANCESTRAL_STATES, ancestralStatesPanel);
        tabbedPane.setToolTipTextAt(index++, "<html>" +
                "Select options for sampling ancestral states at specific<br>" +
                "or all common ancestors, models of counting state changes<br>" +
                "and models of sequencing error for data partitions.</html>");

        tabbedPane.addTab(PRIORS, priorsPanel);
        tabbedPane.setToolTipTextAt(index++, "<html>" +
                "Specify prior probability distributions on each and every<br>" +
                "parameter of the current model.</html>");
        tabbedPane.addTab(OPERATORS, operatorsPanel);
        tabbedPane.setToolTipTextAt(index++, "<html>" +
                "Select and adjust the menu of operators that will be used<br>" +
                "to propose changes to the parameters. Switch off operators<br>" +
                "on certain parameters to fix them to initial values.</html>");
        tabbedPane.addTab(MCMC, mcmcPanel);
        tabbedPane.setToolTipTextAt(index++, "<html>" +
                "Specify the details of MCMC sampling. This includes chain<br>" +
                "length, sampling frequencies, log file names and more.</html>");

        for (int i = 1; i < tabbedPane.getTabCount(); i++) {
            tabbedPane.setEnabledAt(i, false);
        }

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
//        System.out.println("Screen width = " + d.width);
//        System.out.println("Screen height = " + d.height);

        if (d.width < 1000 || d.height < 700) {
            setSize(new java.awt.Dimension(700, 550));
        } else {
            setSize(new java.awt.Dimension(1024, 768));
        }

        if (OSType.isMac()) {
            setMinimumSize(new java.awt.Dimension(640, 480));
        }

        setAllOptions();

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
            if (options.useStarBEAST) {
                speciesSetPanel.setOptions(options);
            } else {
                taxonSetPanel.setOptions(options);
            }
            siteModelsPanel.setOptions(options);
            clockModelsPanel.setOptions(options);
            treesPanel.setOptions(options);
            ancestralStatesPanel.setOptions(options);
            priorsPanel.setOptions(options);
            operatorsPanel.setOptions(options);
            mcmcPanel.setOptions(options);

            setStatusMessage();
        } catch (IllegalArgumentException iae) {
            iae.printStackTrace(System.err);
            JOptionPane.showMessageDialog(this, iae.getMessage(),
                    "Illegal Argument Exception", JOptionPane.ERROR_MESSAGE);
        }

        // enable/disable the other tabs and generate option depending on whether any
        // data has been loaded.
        boolean enabled = options.getDataPartitions().size() > 0;
        for (int i = 1; i < tabbedPane.getTabCount(); i++) {
            tabbedPane.setEnabledAt(i, enabled);
        }
        getExportAction().setEnabled(enabled);
    }

    /**
     * get all the options for all panels
     */
    private void getAllOptions() {
        try {
            dataPanel.getOptions(options);
            tipDatesPanel.getOptions(options);
            traitsPanel.getOptions(options);
            if (options.useStarBEAST) {
                speciesSetPanel.getOptions(options);
            } else {
                taxonSetPanel.getOptions(options);
            }
            siteModelsPanel.getOptions(options);
            clockModelsPanel.getOptions(options);
            treesPanel.getOptions(options);
            ancestralStatesPanel.getOptions(options);
            priorsPanel.getOptions(options);
            operatorsPanel.getOptions(options);
            mcmcPanel.getOptions(options);
        } catch (IllegalArgumentException iae) {
            iae.printStackTrace(System.err);
            JOptionPane.showMessageDialog(this, iae.getMessage(),
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
        if (isDirty() && options.hasData() && isVisible()) {
            int option = JOptionPane.showConfirmDialog(this,
                    "You have made changes but have not generated\n" +
                            "a BEAST XML file. Do you wish to generate\n" +
                            "before closing this window?",
                    "Unused changes",
                    JOptionPane.YES_NO_CANCEL_OPTION,
                    JOptionPane.WARNING_MESSAGE);

            if (option == JOptionPane.YES_OPTION) {
                return doGenerate();
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
        FileInputStream fileIn =
                new FileInputStream(file);
        try {
            ObjectInputStream in = new ObjectInputStream(fileIn);
            options = (BeautiOptions) in.readObject();
            in.close();
        } catch (IOException ioe) {
            JOptionPane.showMessageDialog(this,
                    "Unable to read BEAUti file. BEAUti can only read files\n" +
                            "created by 'Saving' within BEAUti. It cannot read BEAST\n" +
                            "XML files. To read data within BEAST XML files, use\n" +
                            "the 'Import' option.",
                    "Unable to read file",
                    JOptionPane.ERROR_MESSAGE);
        } catch (ClassNotFoundException cnfe) {
            JOptionPane.showMessageDialog(this, "Unable to read BEAUti file: " + cnfe.getMessage(),
                    "Unable to read file",
                    JOptionPane.ERROR_MESSAGE);
            cnfe.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
            return false;
        }
        fileIn.close();

        options.registerComponents(components);
        generator = new BeastGenerator(options, components);

        return true;
    }

    public String getDefaultFileName() {
        return options.fileNameStem + ".beauti";
    }

    protected boolean writeToFile(File file) throws IOException {
        OutputStream fileOut = new FileOutputStream(file);
        ObjectOutputStream out =
                new ObjectOutputStream(fileOut);
        out.writeObject(options);
        out.close();
        fileOut.close();
        return true;
    }

    public final void doImport() {
        File[] files = selectImportFiles("Import Aligment...", true, new FileNameExtensionFilter[] {
                new FileNameExtensionFilter( "Microsatellite (tab-delimited *.txt) Files", "txt"),
                new FileNameExtensionFilter(
                        "NEXUS, BEAST or FASTA Files", "nex", "nexus", "nx", "xml", "beast", "fa", "fasta", "afa")});
        // new FileNameExtensionFilter( "Microsatellite (tab-delimited *.txt) Files", "txt");
        if (files != null && files.length != 0) {
            importFiles(files);
            tabbedPane.setSelectedComponent(dataPanel);
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
                } catch (JDOMException jde) {
                    JOptionPane.showMessageDialog(this, "Error parsing imported file: " + jde,
                            "Error reading file",
                            JOptionPane.ERROR_MESSAGE);
                    jde.printStackTrace();
                }
            }
        }

        if (!options.hasIdenticalTaxa()) {
            setAllOptions(); // need this to refresh panels otherwise it will throw exception
            dataPanel.selectAll();
            dataPanel.unlinkTreeModels();
        }

        setAllOptions();

    }

    public final boolean doImportTraits() {
        if (options.taxonList != null) { // validation of check empty taxonList
            File[] files = selectImportFiles("Import Traits File...", false, new FileNameExtensionFilter[] {
                    new FileNameExtensionFilter("Tab-delimited text files", "txt", "tab", "dat") });

            if (files != null && files.length != 0) {
                try {
                    BEAUTiImporter beautiImporter = new BEAUTiImporter(this, options);
                    beautiImporter.importTraits(files[0]);
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
                    ex.printStackTrace(System.err);
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

            tabbedPane.setSelectedComponent(traitsPanel);
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

        if (options.useStarBEAST && traitName.equalsIgnoreCase(TraitData.TRAIT_SPECIES)) {
            JOptionPane.showMessageDialog(this,
                    "This trait name is already in used to denote species\n" +
                            "for *BEAST. Please select a different name.",
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

    /**
     * Attempts to set up starBEAST - returns true if successful
     * @param useStarBEAST
     * @return
     */
    public boolean setupStarBEAST(boolean useStarBEAST) {
        if (useStarBEAST) {
            if (!options.traitExists(TraitData.TRAIT_SPECIES)) {
                if (!traitsPanel.addTrait(
                        "<html><p>" +
                                "StarBEAST requires a trait to give species designations<br>" +
                                "for each taxon. Create or import a discrete trait<br>" +
                                "labelled 'species'.</p></html>",
                        TraitData.TRAIT_SPECIES,
                        true /* isSpeciesTrait */
                )) {
                    return false;
                }
            } else if (options.getTraitPartitions(options.getTrait(TraitData.TRAIT_SPECIES)).size() > 0) {
                int option = JOptionPane.showConfirmDialog(this,
                        "The trait named '" + TraitData.TRAIT_SPECIES + "', used to denote species in *BEAST, is\n" +
                                "already in use as a data partition. Do you wish to continue?",
                        "Species trait already in use",
                        JOptionPane.YES_NO_OPTION,
                        JOptionPane.WARNING_MESSAGE);
                if (option == JOptionPane.NO_OPTION) {
                    return false;
                }

            }

            dataPanel.selectAll();
            dataPanel.unlinkAll();

            options.starBEASTOptions = new STARBEASTOptions(options);
            options.fileNameStem = "StarBEASTLog";

            tabbedPane.removeTabAt(1);
            tabbedPane.insertTab("Species Sets", null, speciesSetPanel, null, 1);

        } else { // remove species
            options.fileNameStem = MCMCPanel.DEFAULT_FILE_NAME_STEM;

            tabbedPane.removeTabAt(1);
            tabbedPane.insertTab("Taxon Sets", null, taxonSetPanel, null, 1);
        }

        options.useStarBEAST = useStarBEAST;

        treesPanel.updatePriorPanelForSpeciesAnalysis();

        setStatusMessage();

        return true;
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
        treesPanel.setOptions(options); // need this to refresh the currentTreeModel
        return treesPanel.currentTreeModel.getPartitionTreePrior();
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
        } catch (Generator.GeneratorException ge) {
            ge.printStackTrace(System.err);
            JOptionPane.showMessageDialog(this, ge.getMessage(), "Invalid BEAUti setting : ",
                    JOptionPane.ERROR_MESSAGE);
            if (ge.getSwitchToPanel() != null) {
                switchToPanel(ge.getSwitchToPanel());
            }
            return false;
        }

        DefaultPriorTableDialog defaultPriorDialog = new DefaultPriorTableDialog(this);
        if (!defaultPriorDialog.showDialog(options)) {
            return false;
        }

        File file = selectExportFile("Generate BEAST XML File...", new FileNameExtensionFilter("BEAST XML File", "xml", "beast"));

        if (file != null) {
            try {
                getAllOptions();
                generator.generateXML(file);

            } catch (IOException ioe) {
                ioe.printStackTrace(System.err);
                JOptionPane.showMessageDialog(this, "Unable to generate file due to I/O issue: " + ioe.getMessage(),
                        "Unable to generate file", JOptionPane.ERROR_MESSAGE);
                return false;
            } catch (Generator.GeneratorException e) {
                e.printStackTrace(System.err);
                JOptionPane.showMessageDialog(this, "The BEAST XML is incomplete because :\n" + e.getMessage(),
                        "The BEAST XML is incomplete", JOptionPane.ERROR_MESSAGE);
                return false;
            } catch (Exception e) {
                e.printStackTrace(System.err);
                JOptionPane.showMessageDialog(this, "Unable to generate file: " + e.getMessage(),
                        "Unable to generate file", JOptionPane.ERROR_MESSAGE);
                return false;
            }

            clearDirty();
            return true;
        }

        return false;
    }

    /**
     * Use the native file dialog on the Mac because the Swing one is bad. On linux, the native
     * one is bad. No preference on Windows.
     * @param title
     * @return
     */
    public File[] selectImportFiles(final String title, boolean multipleSelection, FileNameExtensionFilter[] fileNameExtensionFilters) {
        if (Boolean.parseBoolean(System.getProperty("use.native.choosers", Boolean.toString(OSType.isMac())))) {
            FileDialog importDialog = fileDialogs.get(title);
            if (importDialog == null) {
                importDialog = new FileDialog(this, title, FileDialog.LOAD);
                fileDialogs.put(title, importDialog);
            }

            importDialog.setVisible(true);
            if (importDialog.getFile() != null) {
                return new File[] { new File(importDialog.getDirectory(), importDialog.getFile()) };
            }
        } else {
            JFileChooser importChooser = fileChoosers.get(title);
            if (importChooser == null) {
                importChooser = new JFileChooser(Utils.getCWD());

                importChooser.setMultiSelectionEnabled(multipleSelection);
                for (FileNameExtensionFilter fileNameExtensionFilter : fileNameExtensionFilters) {
                    importChooser.setFileFilter(fileNameExtensionFilter);
                }
                importChooser.setDialogTitle(title);

                fileChoosers.put(title, importChooser);
            }

            int returnVal = importChooser.showOpenDialog(this);
            if (returnVal == JFileChooser.APPROVE_OPTION) {
                if (importChooser.isMultiSelectionEnabled()) {
                    return importChooser.getSelectedFiles();
                } else {
                    return new File[] { importChooser.getSelectedFile() };
                }
            }
        }

        return null;
    }

    /**
     * Use the native file dialog on the Mac because the Swing one is bad. On linux, the native
     * one is bad. No preference on Windows.
     * @param title
     * @return
     */
    private File selectExportFile(final String title, FileNameExtensionFilter fileNameExtensionFilter) {
        if (Boolean.parseBoolean(System.getProperty("use.native.choosers", Boolean.toString(OSType.isMac())))) {
            FileDialog exportDialog = fileDialogs.get(title);
            if (exportDialog == null) {
                exportDialog = new FileDialog(this, title, FileDialog.SAVE);
                fileDialogs.put(title, exportDialog);
            }

            exportDialog.setFile(options.fileNameStem + ".xml");

            exportDialog.setVisible(true);

            // Mac dialog box will already have asked about overwriting file...
            if (exportDialog.getFile() != null) {
                return new File(exportDialog.getDirectory(), exportDialog.getFile());
            }
        } else {
            JFileChooser exportChooser = fileChoosers.get(title);
            if (exportChooser == null) {
                exportChooser = new JFileChooser(Utils.getCWD());

                // make JFileChooser chooser remember previous path
                exportChooser = new JFileChooser(Utils.getCWD());
                exportChooser.setFileFilter(fileNameExtensionFilter);
                exportChooser.setDialogTitle(title);

                fileChoosers.put(title, exportChooser);
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
                    return file;
                }
            }
        }

        return null;
    }

    public void switchToPanel(String panelName) {
        for (int i = 0; i < tabbedPane.getTabCount(); i++) {
            if (tabbedPane.getTitleAt(i).equals(panelName)) {
                tabbedPane.setSelectedIndex(i);
                break;
            }
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

    protected AbstractAction importTraitsAction = new AbstractAction("Import Traits...") {
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
