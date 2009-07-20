/*
 * BeautiFrame.java
 *
 * (c) 2002-2005 BEAST Development Core Team
 *
 * This package may be distributed under the
 * Lesser Gnu Public Licence (LGPL)
 */
package dr.app.beauti;

import dr.app.beauti.components.ComponentFactory;
import dr.app.beauti.components.SequenceErrorModelComponentFactory;
import dr.app.beauti.components.TipDateSamplingComponentFactory;
import dr.app.beauti.datapanel.DataPanel;
import dr.app.beauti.generator.BeastGenerator;
import dr.app.beauti.mcmcpanel.MCMCPanel;
import dr.app.beauti.modelsPanel.ModelsPanel;
import dr.app.beauti.operatorspanel.OperatorsPanel;
import dr.app.beauti.options.*;
import dr.app.beauti.priorsPanel.PriorsPanel;
import dr.app.beauti.taxonsetspanel.TaxaPanel;
import dr.app.beauti.tipdatepanel.TipDatesPanel;
import dr.app.beauti.traitspanel.TraitsPanel;
import dr.app.beauti.treespanel.OldTreesPanel;
import dr.app.beauti.treespanel.SpeciesTreesPanel;
import dr.app.beauti.treespanel.TreesPanel;
import dr.app.beauti.util.NexusApplicationImporter;
import dr.app.util.Arguments;
import dr.app.util.Utils;
import dr.evolution.alignment.Alignment;
import dr.evolution.alignment.SimpleAlignment;
import dr.evolution.io.Importer;
import dr.evolution.io.NexusImporter;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxa;
import dr.evolution.util.Taxon;
import dr.evolution.util.TaxonList;
import dr.evolution.util.Units;
import dr.java16compat.FileNameExtensionFilter;
import org.jdom.JDOMException;
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
import java.util.EnumSet;
import java.util.List;
import java.util.Map;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: BeautiFrame.java,v 1.22 2006/09/09 16:07:06 rambaut Exp $
 */
public class BeautiFrame extends DocumentFrame {

    /**
     *
     */
    private static final long serialVersionUID = 2114148696789612509L;

    private final BeautiOptions beautiOptions;
    private final BeastGenerator generator;

    private final JTabbedPane tabbedPane = new JTabbedPane();
    public final JLabel statusLabel = new JLabel("No data loaded");

    private DataPanel dataPanel;
    private TipDatesPanel tipDatesPanel;
    private TraitsPanel traitsPanel;
    private TaxaPanel taxaPanel;
    private ModelsPanel modelsPanel;
    private OldTreesPanel oldTreesPanel;
    private SpeciesTreesPanel speciesTreesPanel;
    private TreesPanel treesPanel;
    private PriorsPanel priorsPanel;
    private OperatorsPanel operatorsPanel;
    private MCMCPanel mcmcPanel;

    private BeautiPanel currentPanel;

    private JFileChooser chooser; // make JFileChooser chooser remember previous path

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
                TipDateSamplingComponentFactory.INSTANCE
        };

        beautiOptions = new BeautiOptions(components);
        generator = new BeastGenerator(beautiOptions, components);
    }

    public void initializeComponents() {

        dataPanel = new DataPanel(this, getImportAction(), getDeleteAction());
        tipDatesPanel = new TipDatesPanel(this);
        traitsPanel = new TraitsPanel(this, getImportTraitsAction());
        taxaPanel = new TaxaPanel(this);
        modelsPanel = new ModelsPanel(this, getDeleteAction());
        oldTreesPanel = new OldTreesPanel(this);
        treesPanel = new TreesPanel(this, getDeleteAction());
        speciesTreesPanel = new SpeciesTreesPanel(this);
        priorsPanel = new PriorsPanel(this);
        operatorsPanel = new OperatorsPanel(this);
        mcmcPanel = new MCMCPanel(this);

        tabbedPane.addTab("Data Partitions", dataPanel);
        tabbedPane.addTab("Taxon Sets", taxaPanel);
        tabbedPane.addTab("Tip Dates", tipDatesPanel);
        tabbedPane.addTab("Traits", traitsPanel);
        tabbedPane.addTab("Models", modelsPanel);
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

        getExportAction().setEnabled(false);
        JButton generateButton = new JButton(getExportAction());
        generateButton.putClientProperty("JButton.buttonType", "roundRect");

        JPanel panel2 = new JPanel(new BorderLayout(6, 6));
        panel2.add(statusLabel, BorderLayout.CENTER);
        panel2.add(generateButton, BorderLayout.EAST);

        panel.add(panel2, BorderLayout.SOUTH);

        getContentPane().setLayout(new java.awt.BorderLayout(0, 0));
        getContentPane().add(panel, BorderLayout.CENTER);

        setAllOptions();

        setSize(new java.awt.Dimension(1024, 768));

        // make JFileChooser chooser remember previous path
        chooser = new JFileChooser(Utils.getCWD());
    }

    /**
     * set all the options for all panels
     */
    private void setAllOptions() {
        dataPanel.setOptions(beautiOptions);
        tipDatesPanel.setOptions(beautiOptions);
        traitsPanel.setOptions(beautiOptions);
        taxaPanel.setOptions(beautiOptions);
        modelsPanel.setOptions(beautiOptions);
        if (beautiOptions.isSpeciesAnalysis()) {
            speciesTreesPanel.setOptions(beautiOptions);
        } else if (DataPanel.ALLOW_UNLINKED_TREES) {
            treesPanel.setOptions(beautiOptions);
        } else {
            oldTreesPanel.setOptions(beautiOptions);
        }
        priorsPanel.setOptions(beautiOptions);
        operatorsPanel.setOptions(beautiOptions);
        mcmcPanel.setOptions(beautiOptions);
    }

    /**
     * get all the options for all panels
     */
    private void getAllOptions() {
        dataPanel.getOptions(beautiOptions);
        tipDatesPanel.getOptions(beautiOptions);
        traitsPanel.getOptions(beautiOptions);
        taxaPanel.getOptions(beautiOptions);
        modelsPanel.getOptions(beautiOptions);
        if (beautiOptions.isSpeciesAnalysis()) {
            speciesTreesPanel.getOptions(beautiOptions);
        } else if (DataPanel.ALLOW_UNLINKED_TREES) {
            treesPanel.getOptions(beautiOptions);
        } else {
            oldTreesPanel.getOptions(beautiOptions);
        }
        priorsPanel.getOptions(beautiOptions);
        operatorsPanel.getOptions(beautiOptions);
        mcmcPanel.getOptions(beautiOptions);
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
        chooser.setMultiSelectionEnabled(true);
        FileNameExtensionFilter filter = new FileNameExtensionFilter(
                "NEXUS (*.nex) & BEAST (*.xml) Files", "nex", "nexus", "nx", "xml", "beast");
        chooser.setFileFilter(filter);

        int returnVal = chooser.showDialog(this, "Import Aligment...");
        if (returnVal == JFileChooser.APPROVE_OPTION) {
            File[] files = chooser.getSelectedFiles();
            for (File file : files) {
                if (file == null || file.getName().equals("")) {
                    JOptionPane.showMessageDialog(this, "Invalid file name",
                            "Invalid file name", JOptionPane.ERROR_MESSAGE);
                } else {
                    try {
                        importFromFile(file);

                        setDirty();
                    } catch (FileNotFoundException fnfe) {
                        JOptionPane.showMessageDialog(this, "Unable to open file: File not found",
                                "Unable to open file", JOptionPane.ERROR_MESSAGE);
                    } catch (IOException ioe) {
                        JOptionPane.showMessageDialog(this, "Unable to read file: " + ioe.getMessage(),
                                "Unable to read file", JOptionPane.ERROR_MESSAGE);
                    }
                }
            }
        }
    }

    protected void importFromFile(File file) throws IOException {

        Reader reader = new FileReader(file);

        BufferedReader bufferedReader = new BufferedReader(reader);
        String line = bufferedReader.readLine();
        while (line != null && line.length() == 0) {
            line = bufferedReader.readLine();
        }

        if ((line != null && line.toUpperCase().contains("#NEXUS"))) {
            // is a NEXUS file
            importNexusFile(file);
        } else {
            // assume it is a BEAST XML file and see if that works...
            importBEASTFile(file);
        }

        setStatusMessage();

        setAllOptions();

        // @Todo templates are not implemented yet...
//        getOpenAction().setEnabled(true);
//        getSaveAction().setEnabled(true);
        getExportAction().setEnabled(true);
    }

    protected void importBEASTFile(File file) throws IOException {

        try {
            FileReader reader = new FileReader(file);

            BeastImporter importer = new BeastImporter(reader);

            java.util.List<TaxonList> taxonLists = new ArrayList<TaxonList>();
            java.util.List<Alignment> alignments = new ArrayList<Alignment>();

            importer.importBEAST(taxonLists, alignments);

            TaxonList taxa = taxonLists.get(0);

            for (Alignment alignment : alignments) {
                setData(taxa, alignment, null, null, null, file.getName());
            }
        } catch (JDOMException e) {
            JOptionPane.showMessageDialog(this,
                    "Unable to import file: " + e.getMessage(),
                    "Unable to import file",
                    JOptionPane.WARNING_MESSAGE);
        } catch (Importer.ImportException e) {
            JOptionPane.showMessageDialog(this,
                    "Unable to import file: " + e.getMessage(),
                    "Unable to import file",
                    JOptionPane.WARNING_MESSAGE);
        }

    }

    protected void importNexusFile(File file) throws IOException {

        TaxonList taxa = null;
        SimpleAlignment alignment = null;
        java.util.List<Tree> trees = new ArrayList<Tree>();
        PartitionSubstitutionModel model = null;
        java.util.List<NexusApplicationImporter.CharSet> charSets = new ArrayList<NexusApplicationImporter.CharSet>();

        try {
            FileReader reader = new FileReader(file);


            NexusApplicationImporter importer = new NexusApplicationImporter(reader);

            boolean done = false;

            while (!done) {
                try {

                    NexusImporter.NexusBlock block = importer.findNextBlock();

                    if (block == NexusImporter.TAXA_BLOCK) {

                        if (taxa != null) {
                            throw new NexusImporter.MissingBlockException("TAXA block already defined");
                        }

                        taxa = importer.parseTaxaBlock();

                    } else if (block == NexusImporter.CALIBRATION_BLOCK) {
                        if (taxa == null) {
                            throw new NexusImporter.MissingBlockException("TAXA or DATA block must be defined before a CALIBRATION block");
                        }

                        importer.parseCalibrationBlock(taxa);

                    } else if (block == NexusImporter.CHARACTERS_BLOCK) {

                        if (taxa == null) {
                            throw new NexusImporter.MissingBlockException("TAXA block must be defined before a CHARACTERS block");
                        }

                        if (alignment != null) {
                            throw new NexusImporter.MissingBlockException("CHARACTERS or DATA block already defined");
                        }

                        alignment = (SimpleAlignment) importer.parseCharactersBlock(taxa);

                    } else if (block == NexusImporter.DATA_BLOCK) {

                        if (alignment != null) {
                            throw new NexusImporter.MissingBlockException("CHARACTERS or DATA block already defined");
                        }

                        // A data block doesn't need a taxon block before it
                        // but if one exists then it will use it.
                        alignment = (SimpleAlignment) importer.parseDataBlock(taxa);
                        if (taxa == null) {
                            taxa = alignment;
                        }

                    } else if (block == NexusImporter.TREES_BLOCK) {

                        // I guess there is no reason not to allow multiple trees blocks
//                        if (trees.size() > 0) {
//                            throw new NexusImporter.MissingBlockException("TREES block already defined");
//                        }

                        Tree[] treeArray = importer.parseTreesBlock(taxa);
                        trees.addAll(Arrays.asList(treeArray));

                        if (taxa == null && trees.size() > 0) {
                            taxa = trees.get(0);
                        }


                    } else if (block == NexusApplicationImporter.PAUP_BLOCK) {

                        model = importer.parsePAUPBlock(beautiOptions, charSets);

                    } else if (block == NexusApplicationImporter.MRBAYES_BLOCK) {

                        model = importer.parseMrBayesBlock(beautiOptions, charSets);

                    } else if (block == NexusApplicationImporter.ASSUMPTIONS_BLOCK) {

                        importer.parseAssumptionsBlock(charSets);

                    } else {
                        // Ignore the block..
                    }

                } catch (EOFException ex) {
                    done = true;
                }
            }

            // Allow the user to load taxa only (perhaps from a tree file) so that they can sample from a prior...
            if (alignment == null && taxa == null) {
                throw new NexusImporter.MissingBlockException("TAXON, DATA or CHARACTERS block is missing");
            }

        } catch (Importer.ImportException ime) {
            JOptionPane.showMessageDialog(this, "Error parsing imported file: " + ime,
                    "Error reading file",
                    JOptionPane.ERROR_MESSAGE);
            ime.printStackTrace();
            return;
        } catch (IOException ioex) {
            JOptionPane.showMessageDialog(this, "File I/O Error: " + ioex,
                    "File I/O Error",
                    JOptionPane.ERROR_MESSAGE);
            ioex.printStackTrace();
            return;
        } catch (Exception ex) {
            JOptionPane.showMessageDialog(this, "Fatal exception: " + ex,
                    "Error reading file",
                    JOptionPane.ERROR_MESSAGE);
            ex.printStackTrace();
            return;
        }

        setData(taxa, alignment, trees, model, charSets, file.getName());
    }

    private void setData(
            TaxonList taxa,
            Alignment alignment,
            java.util.List<Tree> trees,
            PartitionSubstitutionModel model,
            java.util.List<NexusApplicationImporter.CharSet> charSets,
            String fileName) {
        String fileNameStem = dr.app.util.Utils.trimExtensions(fileName,
                new String[]{"NEX", "NEXUS", "TRE", "TREE", "XML"});

        if (beautiOptions.taxonList == null) {
            // This is the first partition to be loaded...

            beautiOptions.taxonList = new Taxa(taxa);

            // check the taxon names for invalid characters
            boolean foundAmp = false;
            for (Taxon taxon : taxa) {
                String name = taxon.getId();
                if (name.indexOf('&') >= 0) {
                    foundAmp = true;
                }
            }
            if (foundAmp) {
                JOptionPane.showMessageDialog(this, "One or more taxon names include an illegal character ('&').\n" +
                        "These characters will prevent BEAST from reading the resulting XML file.\n\n" +
                        "Please edit the taxon name(s) before reloading the data file.",
                        "Illegal Taxon Name(s)",
                        JOptionPane.WARNING_MESSAGE);
                return;
            }

            // make sure they all have dates...
            for (int i = 0; i < taxa.getTaxonCount(); i++) {
                if (taxa.getTaxonAttribute(i, "date") == null) {
                    java.util.Date origin = new java.util.Date(0);

                    dr.evolution.util.Date date = dr.evolution.util.Date.createTimeSinceOrigin(0.0, Units.Type.YEARS, origin);
                    taxa.getTaxon(i).setAttribute("date", date);
                }
            }


            beautiOptions.fileNameStem = fileNameStem;
        } else {
            // This is an additional partition so check it uses the same taxa
            if (!beautiOptions.allowDifferentTaxa) { // not allow Different Taxa
                java.util.List<String> oldTaxa = new ArrayList<String>();
                for (int i = 0; i < beautiOptions.taxonList.getTaxonCount(); i++) {
                    oldTaxa.add(beautiOptions.taxonList.getTaxon(i).getId());
                }
                java.util.List<String> newTaxa = new ArrayList<String>();
                for (int i = 0; i < taxa.getTaxonCount(); i++) {
                    newTaxa.add(taxa.getTaxon(i).getId());
                }

                if (!(oldTaxa.containsAll(newTaxa) && oldTaxa.size() == newTaxa.size())) {
                    // AR - Yes and No are perfectly good answers to this question
                    int adt = JOptionPane.showOptionDialog(this,
                            "This file contains different taxa from the previously loaded\n" +
                                    "data partitions. This may be because the taxa are mislabelled\n" +
                                    "and need correcting before reloading.\n\n" +
                                    "Would you like to allow different taxa for each partition?\n",
                            "Validation of Non-matching Taxon Name(s)",
                            JOptionPane.YES_NO_OPTION,
                            JOptionPane.WARNING_MESSAGE,
                            null,
                            new String[]{"Yes", "No"},
                            "No"
                    ); // default button title

                    if (adt == JOptionPane.YES_OPTION) {
                        // set to Allow Different Taxa
                        beautiOptions.allowDifferentTaxa = true;
                        //changeTabs();// can be added, if required in future

                        java.util.List<String> prevTaxa = new ArrayList<String>();
                        for (int i = 0; i < beautiOptions.taxonList.getTaxonCount(); i++) {
                            prevTaxa.add(beautiOptions.taxonList.getTaxon(i).getId());
                        }
                        for (int i = 0; i < taxa.getTaxonCount(); i++) {
                            if (!prevTaxa.contains(taxa.getTaxon(i).getId())) {
                                beautiOptions.taxonList.addTaxon(taxa.getTaxon(i));
                            }
                        }

                    } else {
                        return;
                    }
                }
            } else { // allow Different Taxa
                // AR - it will be much simpler just to consider beautiOptions.taxonList
                // to be the union set of all taxa. Each data partition has an alignment
                // which is a taxon list containing the taxa specific to that partition

//                for (Taxon taxon : taxa) { // not working
//                    if (!beautiOptions.taxonList.contains(taxon)) {
//                        beautiOptions.taxonList.addTaxon(taxon);
//                    }
//                }

                // add the new diff taxa
                java.util.List<String> prevTaxa = new ArrayList<String>();
                for (int i = 0; i < beautiOptions.taxonList.getTaxonCount(); i++) {
                    prevTaxa.add(beautiOptions.taxonList.getTaxon(i).getId());
                }
                for (int i = 0; i < taxa.getTaxonCount(); i++) {
                    if (!prevTaxa.contains(taxa.getTaxon(i).getId())) {
                        beautiOptions.taxonList.addTaxon(taxa.getTaxon(i));
                    }
                }

            }
        }

        addAlignment(alignment, charSets, model, fileName, fileNameStem);

        addTrees(trees);
    }

    private void addAlignment(Alignment alignment, java.util.List<NexusApplicationImporter.CharSet> charSets,
                              PartitionSubstitutionModel model,
                              String fileName, String fileNameStem) {
        if (alignment != null) {
            java.util.List<PartitionData> partitions = new ArrayList<PartitionData>();
            if (charSets != null && charSets.size() > 0) {
                for (NexusApplicationImporter.CharSet charSet : charSets) {
                    partitions.add(new PartitionData(charSet.getName(), fileName,
                            alignment, charSet.getFromSite(), charSet.getToSite(), charSet.getEvery()));
                }
            } else {
                partitions.add(new PartitionData(fileNameStem, fileName, alignment));
            }
            for (PartitionData partition : partitions) {
                beautiOptions.dataPartitions.add(partition);
                //TODO Cannot load Substitution Model and Tree Model from BEAST file yet
                if (model != null) {
                    partition.setPartitionSubstitutionModel(model);
//                    beautiOptions.addPartitionSubstitutionModel(model);
                    
                    // use same tree model and same tree prior in beginning
                    for (PartitionTreeModel ptm : beautiOptions.getPartitionTreeModels()) {                        
                    	partition.setPartitionTreeModel(ptm);
                    	
                    	if (ptm.getPartitionTreePrior() == null || 
                    			!(ptm.getPartitionTreePrior().getName().equalsIgnoreCase(beautiOptions.activedSameTreePrior.getName()))) {
                    		PartitionTreePrior ptp = new PartitionTreePrior(beautiOptions, ptm);
                            ptm.setPartitionTreePrior(ptp);
                    	}
                    }
                    if (partition.getPartitionTreeModel() == null) {
                    	// PartitionTreeModel based on PartitionData
                    	PartitionTreeModel ptm = new PartitionTreeModel(beautiOptions, partition);
                        partition.setPartitionTreeModel(ptm);
                        
                        // PartitionTreePrior always based on PartitionTreeModel
                        PartitionTreePrior ptp = new PartitionTreePrior(beautiOptions, ptm);
                        ptm.setPartitionTreePrior(ptp);
                        
//                        beautiOptions.addPartitionTreeModel(ptm);
                        beautiOptions.shareSameTreePrior = true;
                        beautiOptions.activedSameTreePrior = ptp;                        
                    }

                	// use same clock model in beginning, have to create after partition.setPartitionTreeModel(ptm);
                	for (PartitionClockModel pcm : beautiOptions.getPartitionClockModels()) {                        
                        partition.setPartitionClockModel(pcm);                        
                	}
                	if (partition.getPartitionClockModel() == null) {
                        // PartitionClockModel based on PartitionData
                		PartitionClockModel pcm = new PartitionClockModel(beautiOptions, partition);
                        partition.setPartitionClockModel(pcm);
//                        beautiOptions.addPartitionClockModel(pcm);
                    }
                } else {// only this works                    
                	for (PartitionSubstitutionModel psm : beautiOptions.getPartitionSubstitutionModels()) {
                        if (psm.getDataType() == alignment.getDataType()) { // use same substitution model in beginning
                            partition.setPartitionSubstitutionModel(psm);
                        }
                	}
                	if (partition.getPartitionSubstitutionModel() == null) {
                        // PartitionSubstitutionModel based on PartitionData
                    	PartitionSubstitutionModel psm = new PartitionSubstitutionModel(beautiOptions, partition);
                        partition.setPartitionSubstitutionModel(psm);
//                        beautiOptions.addPartitionSubstitutionModel(psm);
                    }
                	
                	// use same tree model and same tree prior in beginning
                    for (PartitionTreeModel ptm : beautiOptions.getPartitionTreeModels()) {                        
                    	partition.setPartitionTreeModel(ptm);
                    	
                    	if (ptm.getPartitionTreePrior() == null || 
                    			!(ptm.getPartitionTreePrior().getName().equalsIgnoreCase(beautiOptions.activedSameTreePrior.getName()))) {
                    		PartitionTreePrior ptp = new PartitionTreePrior(beautiOptions, ptm);
                            ptm.setPartitionTreePrior(ptp);
                    	}
                    }
                    if (partition.getPartitionTreeModel() == null) {
                    	// PartitionTreeModel based on PartitionData
                    	PartitionTreeModel ptm = new PartitionTreeModel(beautiOptions, partition);
                        partition.setPartitionTreeModel(ptm);
                        
                        // PartitionTreePrior always based on PartitionTreeModel
                        PartitionTreePrior ptp = new PartitionTreePrior(beautiOptions, ptm);
                        ptm.setPartitionTreePrior(ptp);
                        
//                        beautiOptions.addPartitionTreeModel(ptm);
                        beautiOptions.shareSameTreePrior = true;
                        beautiOptions.activedSameTreePrior = ptp;                        
                    }

                	// use same clock model in beginning, have to create after partition.setPartitionTreeModel(ptm);
                	for (PartitionClockModel pcm : beautiOptions.getPartitionClockModels()) {                        
                        partition.setPartitionClockModel(pcm);                        
                	}
                	if (partition.getPartitionClockModel() == null) {
                        // PartitionClockModel based on PartitionData
                		PartitionClockModel pcm = new PartitionClockModel(beautiOptions, partition);
                        partition.setPartitionClockModel(pcm);
//                        beautiOptions.addPartitionClockModel(pcm);
                    }
                    
                }
            }
            beautiOptions.updatePartitionClockTreeLinks();
        }
    }

    private void addTrees(java.util.List<Tree> trees) {
        if (trees != null && trees.size() > 0) {
            for (Tree tree : trees) {
                String id = tree.getId();
                if (id == null || id.trim().length() == 0) {
                    tree.setId("tree_" + (beautiOptions.userTrees.size() + 1));
                } else {
                    String newId = id;
                    int count = 1;
                    for (Tree tree1 : beautiOptions.userTrees) {
                        if (tree1.getId().equals(newId)) {
                            newId = id + "_" + count;
                            count++;
                        }
                    }
                    tree.setId(newId);
                }
                beautiOptions.userTrees.add(tree);
            }
        }
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

        if (beautiOptions.isSpeciesAnalysis()) { // species
            setupSpeciesAnalysis();
        }

        setAllOptions();
    }

    public void setupSpeciesAnalysis() {
        dataPanel.selectAll();
        dataPanel.unlinkModels();
        //TODO tree
        for (PartitionSubstitutionModel model : beautiOptions.getPartitionSubstitutionModels()) {
//        	model.initAllParametersAndOperators();
        }

        int i = tabbedPane.indexOfTab("Trees");
        tabbedPane.removeTabAt(i);
        tabbedPane.insertTab("Trees", null, speciesTreesPanel, "", i);
        speciesTreesPanel.getOptions(beautiOptions);

        beautiOptions.initSpeciesParametersAndOperators();
        beautiOptions.fileNameStem = "LogStem";

        setStatusMessage();
    }

    public void removeSepciesAnalysisSetup() {
        int i = tabbedPane.indexOfTab("Trees");
        tabbedPane.removeTabAt(i);
        if (DataPanel.ALLOW_UNLINKED_TREES) {
            tabbedPane.insertTab("Trees", null, treesPanel, "", i);
        } else {
            tabbedPane.insertTab("Trees", null, oldTreesPanel, "", i);
        }
        
//        beautiOptions.
        
        setStatusMessage();
    }

    public void removeSpecifiedTreePrior(boolean isChecked) { // TipDatesPanel usingTipDates
    	//TODO: wait for new implementation in BEAST
    	   	
	    if (DataPanel.ALLOW_UNLINKED_TREES) {
			treesPanel.setCheckedTipDate(isChecked);
		} else {
			if (isChecked) {
				oldTreesPanel.treePriorCombo.removeItem(TreePrior.YULE);
				oldTreesPanel.treePriorCombo.removeItem(TreePrior.BIRTH_DEATH);
			} else {
				oldTreesPanel.treePriorCombo = new JComboBox(EnumSet.range(TreePrior.CONSTANT, TreePrior.BIRTH_DEATH).toArray());
			}
		}    	
    }
    
    private void importMultiTraits(final File file) throws IOException {
//        if( beautiOptions.taxonList == null ) {
//             JOptionPane.showMessageDialog(this, "No taxa loaded yet - noting done!",
//                            "No taxa loaded", JOptionPane.ERROR_MESSAGE);
//            return;
//        } // move to doImportTraits()

        try {
            Map<String, List<String[]>> traits = Utils.importTraitsFromFile(file, "\t");

            for (Map.Entry<String, List<String[]>> e : traits.entrySet()) {
                final List<String[]> value = e.getValue();
                final Class c = Utils.detectTYpe(e.getValue().get(0)[1]);
                final String label = e.getKey();

                Boolean warningGiven = false;
                for (String[] v : e.getValue()) {
                    final Class c1 = Utils.detectTYpe(v[1]);
                    if (c != c1 && !warningGiven) {
                        JOptionPane.showMessageDialog(this, "Not all values of same type in column" + label,
                                "Incompatible values", JOptionPane.WARNING_MESSAGE);
                        warningGiven = true;
                        // TODO Error - not all values of same type
                    }
                }

                beautiOptions.selecetedTraits.add(label);
                TraitGuesser.TraitType t = (c == Boolean.class || c == String.class) ? TraitGuesser.TraitType.DISCRETE :
                        (c == Integer.class) ? TraitGuesser.TraitType.INTEGER : TraitGuesser.TraitType.CONTINUOUS;

                beautiOptions.traitTypes.put(label, t);

                for (final String[] v : e.getValue()) {
                    final int index = beautiOptions.taxonList.getTaxonIndex(v[0]);
                    if (index >= 0) {
                        // if the taxon isn't in the list then ignore it.
                        // TODO provide a warning of unmatched taxa
                        final Taxon taxon = beautiOptions.taxonList.getTaxon(index);
                        taxon.setAttribute(label, Utils.constructFromString(c, v[1]));
                    }
                }
            }
        } catch (Arguments.ArgumentException e) {
            JOptionPane.showMessageDialog(this, "Error in loading traits file: " + e.getMessage(),
                    "Error Loading file", JOptionPane.ERROR_MESSAGE);
        }
    }

    private void setStatusMessage() {
        String message = "";
        if (beautiOptions.dataPartitions.size() > 0) {
            message += "Data: " + beautiOptions.taxonList.getTaxonCount() + " taxa, " +
                    beautiOptions.dataPartitions.size() +
                    (beautiOptions.dataPartitions.size() > 1 ? " partitions" : " partition");

            if (beautiOptions.isSpeciesAnalysis()) {
                int num = beautiOptions.getSpeciesList().size();
                message += ", " + num + " species"; // species is both singular and plural
            }

            if (beautiOptions.userTrees.size() > 0) {
                message += ", " + beautiOptions.userTrees.size() +
                        (beautiOptions.userTrees.size() > 1 ? " trees" : " tree");
            }
            
            if (beautiOptions.allowDifferentTaxa) {
            	message += " in total";
            }

            if (beautiOptions.isSpeciesAnalysis()) {
                message += ";    Welcome to *BEAST";
            }
        } else if (beautiOptions.userTrees.size() > 0) {
            message += "Trees only : " + beautiOptions.userTrees.size() +
                    (beautiOptions.userTrees.size() > 1 ? " trees, " : " tree, ") +
                    beautiOptions.taxonList.getTaxonCount() + " taxa";
        } else if (beautiOptions.taxonList != null && beautiOptions.taxonList.getTaxonCount() > 0) {
            message += "Taxa only: " + beautiOptions.taxonList.getTaxonCount() + " taxa";
        }

        statusLabel.setText(message);
    }

    public final boolean doGenerate() {

        try {
            generator.checkOptions();
        } catch (IllegalArgumentException iae) {
            JOptionPane.showMessageDialog(this, iae.getMessage(),
                    "Unable to generate file",
                    JOptionPane.ERROR_MESSAGE);
            return false;
        }

        // todo save last dir, sync import and out directory, offer stem as default

        FileNameExtensionFilter filter = new FileNameExtensionFilter("Generate BEAST File...", "xml", "beast");
        chooser.setFileFilter(filter);
      
        final int returnVal = chooser.showSaveDialog(this);
        if( returnVal == JFileChooser.APPROVE_OPTION ) {
            File file = chooser.getSelectedFile();
            try {
                generate(file);

            } catch (IOException ioe) {
                JOptionPane.showMessageDialog(this, "Unable to generate file: " + ioe.getMessage(),
                        "Unable to generate file",
                        JOptionPane.ERROR_MESSAGE);
                return false;
            }
        }

//        FileDialog dialog = new FileDialog(this,
//                "Generate BEAST File...",
//                FileDialog.SAVE);
//
//        dialog.setVisible(true);
//        if (dialog.getFile() != null) {
//            File file = new File(dialog.getDirectory(), dialog.getFile());
//
//            try {
//                generate(file);
//
//            } catch (IOException ioe) {
//                JOptionPane.showMessageDialog(this, "Unable to generate file: " + ioe.getMessage(),
//                        "Unable to generate file",
//                        JOptionPane.ERROR_MESSAGE);
//                return false;
//            }
//        }
        
        clearDirty();
        return true;
    }

    protected void generate(File file) throws IOException {
        getAllOptions();

        FileWriter fw = new FileWriter(file);
        generator.generateXML(fw);
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

    protected AbstractAction importAlignmentAction = new AbstractAction("Import Alignment...") {
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
