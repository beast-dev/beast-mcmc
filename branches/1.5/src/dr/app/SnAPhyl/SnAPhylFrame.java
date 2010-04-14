/*
 * BeautiFrame.java
 *
 * (c) 2002-2005 BEAST Development Core Team
 *
 * This package may be distributed under the
 * Lesser Gnu Public Licence (LGPL)
 */
package dr.app.SnAPhyl;

import dr.app.SnAPhyl.datapanel.DataPanel;
import dr.app.SnAPhyl.initializationpanel.InitializationPanel;
import dr.app.SnAPhyl.mcmcpanel.MCMCPanel;
import dr.app.SnAPhyl.operatorspanel.OperatorsPanel;
import dr.app.SnAPhyl.priorsPanel.DefaultPriorDialog;
import dr.app.SnAPhyl.priorsPanel.PriorsPanel;
import dr.app.SnAPhyl.util.BEAUTiImporter;
import dr.app.beauti.BeautiFrame;
import dr.app.beauti.BeautiPanel;
import dr.app.beauti.components.ComponentFactory;
import dr.app.beauti.components.SequenceErrorModelComponentFactory;
import dr.app.beauti.components.TipDateSamplingComponentFactory;
import dr.app.beauti.generator.BeastGenerator;
import dr.app.beauti.options.BeautiOptions;
import dr.app.java16compat.FileNameExtensionFilter;
import dr.app.util.Utils;
import dr.evolution.io.Importer.ImportException;
import dr.evolution.io.NexusImporter.MissingBlockException;
import org.virion.jam.util.IconUtils;

import javax.swing.*;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import javax.swing.plaf.BorderUIResource;
import java.awt.*;
import java.io.File;
import java.io.IOException;


/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: BeautiFrame.java,v 1.22 2006/09/09 16:07:06 rambaut Exp $
 */
public class SnAPhylFrame extends BeautiFrame {

    private static final long serialVersionUID = 2114148696789612509L;

    private final BeautiOptions beautiOptions;
    private final BeastGenerator generator;

    private final JTabbedPane tabbedPane = new JTabbedPane();
    public final JLabel statusLabel = new JLabel("No data loaded");

    private DataPanel dataPanel;
    private PriorsPanel priorsPanel;
    private InitializationPanel initializationPanel;
    private OperatorsPanel operatorsPanel;
    private MCMCPanel mcmcPanel;

    private BeautiPanel currentPanel;

    private JFileChooser importChooser; // make JFileChooser chooser remember previous path
    private JFileChooser exportChooser; // make JFileChooser chooser remember previous path

    final Icon gearIcon = IconUtils.getIcon(this.getClass(), "images/gear.png");

    public SnAPhylFrame(String title) {
        super(title);

        ComponentFactory[] components = {
                SequenceErrorModelComponentFactory.INSTANCE,
                TipDateSamplingComponentFactory.INSTANCE
        };

        beautiOptions = new BeautiOptions(components);
        generator = new BeastGenerator(beautiOptions, components);
    }

    public void initializeComponents() {

        dataPanel = new DataPanel(this, getImportAction());
        priorsPanel = new PriorsPanel(this, false);
        initializationPanel = new InitializationPanel(this);
        operatorsPanel = new OperatorsPanel(this);
        mcmcPanel = new MCMCPanel(this);

        tabbedPane.addTab("Data Partitions", dataPanel);        
        tabbedPane.addTab("Priors", priorsPanel);
        tabbedPane.addTab("Initialization", initializationPanel);
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
                "NEXUS (*.nex) & BEAST (*.xml) Files", "nex", "nexus", "nx", "xml", "beast"));
        importChooser.setDialogTitle("Import Aligment...");
    }

    /**
     * set all the options for all panels
     */
    public void setAllOptions() {
        dataPanel.setOptions(beautiOptions);
        priorsPanel.setOptions(beautiOptions);
        initializationPanel.setOptions(beautiOptions);
        operatorsPanel.setOptions(beautiOptions);
        mcmcPanel.setOptions(beautiOptions);

        setStatusMessage();
    }

    /**
     * get all the options for all panels
     */
    private void getAllOptions() {
        dataPanel.getOptions(beautiOptions);          
        priorsPanel.getOptions(beautiOptions);
        initializationPanel.getOptions(beautiOptions);
        operatorsPanel.getOptions(beautiOptions);
        mcmcPanel.getOptions(beautiOptions);
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
                return !doGenerateSnAPhyl();
            } else if (option == JOptionPane.CANCEL_OPTION || option == JOptionPane.DEFAULT_OPTION) {
                return false;
            }
            return true;
        }
        return true;
    }

    public String getDefaultFileName() {
        return beautiOptions.fileNameStem + ".beauti";
    }

    public final void doImportSnAPhyl() {
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
                        beautiImporter.importFromFile(file);

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

            setStatusMessage();
            setAllOptions();

//          // @Todo templates are not implemented yet...
////        getOpenAction().setEnabled(true);
////        getSaveAction().setEnabled(true);
            getExportAction().setEnabled(true);
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

    public final boolean doGenerateSnAPhyl() {

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
            		     "Overwrite the exsting file?", JOptionPane.YES_NO_OPTION);
            }

            if (n == JOptionPane.YES_OPTION) {
	            try {
                    getAllOptions();
	                generator.generateXML(file);

	            } catch (IOException ioe) {
	                JOptionPane.showMessageDialog(this, "Unable to generate file: " + ioe.getMessage(),
	                        "Unable to generate file", JOptionPane.ERROR_MESSAGE);
	                return false;
	            } catch (Exception ioe) {
	                JOptionPane.showMessageDialog(this, "Unable to generate file: " + ioe.getMessage(),
	                        "Unable to generate file", JOptionPane.ERROR_MESSAGE);
	                return false;
	            }
            } else {
            	doGenerateSnAPhyl();
            }
        }

        clearDirty();
        return true;
    }

    public Action getImportAction() {
        return importAlignmentAction;
    }

    protected AbstractAction importAlignmentAction = new AbstractAction("Import Alignment...") {
        private static final long serialVersionUID = 3217702096314745005L;

        public void actionPerformed(java.awt.event.ActionEvent ae) {
            doImportSnAPhyl();
        }
    };

    public Action getExportAction() {
        return generateAction;
    }

    protected AbstractAction generateAction = new AbstractAction("Generate BEAST File...", gearIcon) {
        private static final long serialVersionUID = -5329102618630268783L;

        public void actionPerformed(java.awt.event.ActionEvent ae) {
            doGenerateSnAPhyl();
        }
    };

}
