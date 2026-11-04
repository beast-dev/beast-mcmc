/*
 * BeastDialog.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
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
 *
 */

package dr.app.beast;

import beagle.BeagleInfo;
import dr.app.checkpoint.BeastCheckpointer;
import dr.app.gui.FileDrop;
import dr.app.gui.components.WholeNumberField;
import dr.app.util.OSType;
import dr.util.OutputFileHelper;
import jam.html.SimpleLinkListener;
import jam.panels.OptionsPanel;

import javax.swing.*;
import javax.swing.border.Border;
import javax.swing.border.EmptyBorder;
import javax.swing.border.TitledBorder;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyEvent;
import java.io.File;
import java.io.IOException;

/**
 * @author Andrew Rambaut
 * @author Guy Baele
 */

public class BeastDialog {
    private final JFrame frame;

    private final OptionsPanel optionPanel;

    private final WholeNumberField seedText = new WholeNumberField((long) 1, Long.MAX_VALUE);
    private final JCheckBox overwriteCheckBox = new JCheckBox("Allow overwriting of log files");
    private final JCheckBox outputDirectoryCheckBox = new JCheckBox("Write all output files to a folder:");
    private final JTextField outputDirectoryNameText = new JTextField("not selected", 16);
    private final JButton outputDirectoryButton = new JButton("Choose Folder...");
    private final JCheckBox beagleInfoCheckBox = new JCheckBox("Show list of available BEAGLE resources and Quit");
    private final JComboBox beagleResourceCombo = new JComboBox(new Object[]{"CPU", "GPU"});
    private final JCheckBox beagleSSECheckBox = new JCheckBox("Use CPU's SSE extensions when possible");
    private final JComboBox beaglePrecisionCombo = new JComboBox(new Object[]{"Double", "Single"});
    private final JComboBox beagleScalingCombo = new JComboBox(new Object[]{"Default", "Dynamic", "Delayed", "Always", "Never"});

    private final JComboBox threadsCombo = new JComboBox(new Object[]{"Automatic", 0, 1, 2, 3, 4, 5, 6, 7, 8});

    private File inputFile = null;
    private File outputDirectory = null;

    public BeastDialog(final JFrame frame, final String titleString, final Icon icon) {
        this.frame = frame;

        optionPanel = new OptionsPanel(12, 12);

        //this.frame = frame;

        JPanel panel = new JPanel(new BorderLayout());
        panel.setOpaque(false);

        final OptionsPanel optionPanel3 = new OptionsPanel(0, 3);

        final JLabel titleIcon = new JLabel();
        titleIcon.setIcon(icon);

        final JEditorPane titleText = new JEditorPane("text/html", "<html>" + titleString + "</html>");
        titleText.setOpaque(false);
        titleText.setEditable(false);
        titleText.addHyperlinkListener(new SimpleLinkListener());
        optionPanel3.addComponent(titleText);

//        final JButton aboutButton = new JButton("About BEAST...");
//        //aboutButton.setAction();
//        optionPanel3.addComponent(aboutButton);

        optionPanel.addComponents(titleIcon, optionPanel3);

        final JButton inputFileButton = new JButton("Choose File...");
        final JTextField inputFileNameText = new JTextField("not selected", 16);

        inputFileButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent ae) {
                FileDialog dialog = new FileDialog(frame,
                        "Select target file...",
                        FileDialog.LOAD);

                dialog.setVisible(true);
                if (dialog.getFile() == null) {
                    // the dialog was cancelled...
                    return;
                }

                inputFile = new File(dialog.getDirectory(), dialog.getFile());
                inputFileNameText.setText(inputFile.getName());

            }
        });
        inputFileNameText.setEditable(false);

        JPanel panel1 = new JPanel(new BorderLayout(0, 0));
        panel1.add(inputFileNameText, BorderLayout.CENTER);
        panel1.add(inputFileButton, BorderLayout.EAST);
        inputFileNameText.setToolTipText(
                "<html>Drag a BEAST XML file here or use the button to<br>" +
                        "select one from a file dialog box.</html>");
        inputFileButton.setToolTipText(
                "<html>Drag a BEAST XML file here or use the button to<br>" +
                        "select one from a file dialog box.</html>");
        optionPanel.addComponentWithLabel("BEAST XML File: ", panel1);

        Color focusColor = UIManager.getColor("Focus.color");
        Border focusBorder = BorderFactory.createMatteBorder(2, 2, 2, 2, focusColor);
        new FileDrop(null, inputFileNameText, focusBorder, new FileDrop.Listener() {
            public void filesDropped(java.io.File[] files) {
                inputFile = files[0];
                inputFileNameText.setText(inputFile.getName());
            }   // end filesDropped
        }); // end FileDrop.Listener

        overwriteCheckBox.setToolTipText(
                "<html>Specify whether BEAST will overwrite existing log files<br>" +
                      "with the same name.</html>");
        optionPanel.addComponent(overwriteCheckBox);

        outputDirectoryCheckBox.setToolTipText(
                "<html>By default, log, tree and checkpoint files are written to the<br>" +
                        "folder containing the BEAST XML file (or to the folders named<br>" +
                        "in the XML). Select this to write them all to a chosen folder<br>" +
                        "instead (equivalent to the -output command line option).</html>");
        optionPanel.addComponent(outputDirectoryCheckBox);

        outputDirectoryButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent ae) {
                File directory = chooseDirectory(frame, "Select output folder...", outputDirectory);
                if (directory != null) {
                    setOutputDirectory(directory);
                }
            }
        });
        outputDirectoryNameText.setEditable(false);

        JPanel outputDirectoryPanel = new JPanel(new BorderLayout(0, 0));
        outputDirectoryPanel.add(outputDirectoryNameText, BorderLayout.CENTER);
        outputDirectoryPanel.add(outputDirectoryButton, BorderLayout.EAST);
        outputDirectoryNameText.setToolTipText(
                "<html>Drag a folder here or use the button to select one<br>" +
                        "from a file dialog box.</html>");
        outputDirectoryButton.setToolTipText(
                "<html>Drag a folder here or use the button to select one<br>" +
                        "from a file dialog box.</html>");
        final JLabel outputDirectoryLabel = optionPanel.addComponentWithLabel("Output folder: ", outputDirectoryPanel);

        new FileDrop(null, outputDirectoryNameText, focusBorder, new FileDrop.Listener() {
            public void filesDropped(java.io.File[] files) {
                File dropped = files[0];
                if (!dropped.isDirectory()) {
                    // a file was dropped: use the folder it is in
                    dropped = dropped.getAbsoluteFile().getParentFile();
                }
                if (dropped != null) {
                    outputDirectoryCheckBox.setSelected(true);
                    setOutputDirectory(dropped);
                }
            }
        });

        outputDirectoryCheckBox.addChangeListener(new ChangeListener() {
            public void stateChanged(ChangeEvent e) {
                boolean selected = outputDirectoryCheckBox.isSelected();
                outputDirectoryLabel.setEnabled(selected);
                outputDirectoryNameText.setEnabled(selected);
                outputDirectoryButton.setEnabled(selected);
            }
        });
        // the first time the option is turned on, ask for the folder straight away; once
        // a folder has been chosen the option can be turned on and off without changing it
        outputDirectoryCheckBox.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent ae) {
                if (outputDirectoryCheckBox.isSelected() && outputDirectory == null) {
                    File directory = chooseDirectory(frame, "Select output folder...", null);
                    if (directory != null) {
                        setOutputDirectory(directory);
                    } else {
                        // cancelled without choosing a folder so turn the option off again
                        outputDirectoryCheckBox.setSelected(false);
                    }
                }
            }
        });
        outputDirectoryLabel.setEnabled(false);
        outputDirectoryNameText.setEnabled(false);
        outputDirectoryButton.setEnabled(false);

        optionPanel.addSeparator();

        final JButton chkptButton = new JButton("Custom settings");
        chkptButton.setToolTipText("<html>By default, a checkpoint file will be written according to the<br>" +
                "specifications in your XML file. No previous checkpointed file<br>" +
                "will be loaded and no custom file name can be provided.</html>");

        chkptButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                JDialog dialog = new JDialog(frame, "Checkpointing settings");
                dialog.setLocationRelativeTo(chkptButton);
                dialog.setModal(true);
                dialog.setAlwaysOnTop(true);

                OptionsPanel chkptPanel = new OptionsPanel(0, 3);
                //chkptPanel.setLayout(new BorderLayout());
                chkptPanel.setOpaque(false);
                JCheckBox overruleXML = new JCheckBox();
                overruleXML.setToolTipText("This will ignore the checkpointing settings in your XML file.");
                overruleXML.addActionListener(new ActionListener() {
                    @Override
                    public void actionPerformed(ActionEvent actionEvent) {
                        if (overruleXML.isSelected()) {
                            System.setProperty("checkpointOverrule", Boolean.TRUE.toString());
                        } else {
                            System.setProperty("checkpointOverrule", Boolean.FALSE.toString());
                        }
                    }
                });
                chkptPanel.addComponentWithLabel("Overrule XML checkpointing settings", overruleXML);

                chkptPanel.addSeparator();

                JTextField checkpointInput = new JTextField(15);
                checkpointInput.addKeyListener(new java.awt.event.KeyListener() {
                    public void keyTyped(KeyEvent e) {
                    }

                    public void keyPressed(KeyEvent e) {
                    }

                    public void keyReleased(KeyEvent e) {
                        System.setProperty(BeastCheckpointer.LOAD_STATE_FILE, checkpointInput.getText());
                    }
                });

                JTextField checkpointOutput = new JTextField(15);
                checkpointOutput.addKeyListener(new java.awt.event.KeyListener() {
                    public void keyTyped(KeyEvent e) {
                    }

                    public void keyPressed(KeyEvent e) {
                    }

                    public void keyReleased(KeyEvent e) {
                        System.setProperty(BeastCheckpointer.SAVE_STATE_FILE, checkpointOutput.getText());
                    }
                });

                WholeNumberField checkpointEvery = new WholeNumberField(1, Integer.MAX_VALUE);
                checkpointEvery.addKeyListener(new java.awt.event.KeyListener() {
                    public void keyTyped(KeyEvent e) {
                    }

                    public void keyPressed(KeyEvent e) {
                    }

                    public void keyReleased(KeyEvent e) {
                        System.setProperty(BeastCheckpointer.SAVE_STATE_EVERY, checkpointEvery.getValue() + "");
                    }
                });
                checkpointEvery.setValue(1000000);
                checkpointEvery.setColumns(10);

                chkptPanel.addComponentWithLabel("Load previous checkpoint file: ", checkpointInput);
                chkptPanel.addComponentWithLabel("Save new checkpoint file (will overwrite): ", checkpointOutput);
                chkptPanel.addComponentWithLabel("Save checkpoint every: ", checkpointEvery);

                dialog.add(chkptPanel);
                dialog.pack();
                dialog.setVisible(true);
            }
        });

        JPanel chkptPanel = new JPanel(new BorderLayout(0, 0));
        chkptPanel.add(chkptButton, BorderLayout.CENTER);
        optionPanel.addComponentWithLabel("Checkpointing:", chkptPanel);

        optionPanel.addSeparator();

        seedText.setColumns(12);
        seedText.setToolTipText(
                "<html>Specify a particular random number seed to replicate<br>" +
                "precisely the sequence of steps in the MCMC chain. By<br>" +
                "default this uses system information to provide a new<br>" +
                        "seed each run.</html>");

        optionPanel.addComponentWithLabel("Random number seed: ", seedText);

        threadsCombo.setToolTipText(
                "<html>Specify how large a thread pool to use.<br>" +
                      "In most circumstances this should be set to 'automatic'<br>" +
                      "but in some circumstances it may be desirable to restict<br>" +
                      "the number of cores being used. 0 will turn off threading</html>");

        optionPanel.addComponentWithLabel("Thread pool size: ", threadsCombo);

        optionPanel.addSeparator();

        // Only offer the BEAGLE options if the library can be loaded; otherwise explain
        // that it is needed and where to get it.
        final boolean beagleInstalled = BeagleInfo.getVersionNumbers().length > 0;

        if (beagleInstalled) {
            final JLabel beagleLabel = new JLabel("BEAGLE library v" + BeagleInfo.getVersion() + " options:");
            optionPanel.addSpanningComponent(beagleLabel);

            final OptionsPanel optionPanel1 = new OptionsPanel(0, 6);
            optionPanel1.setBorder(BorderFactory.createEmptyBorder());

            OptionsPanel optionPanel2 = new OptionsPanel(0, 3);
            optionPanel2.setBorder(BorderFactory.createEmptyBorder());
            optionPanel2.addComponentWithLabel("Prefer use of: ", beagleResourceCombo);
            optionPanel2.addComponent(beagleSSECheckBox);
            beagleSSECheckBox.setSelected(true);
            optionPanel2.addComponentWithLabel("Prefer precision: ", beaglePrecisionCombo);
            optionPanel2.addComponentWithLabel("Rescaling scheme: ", beagleScalingCombo);
            optionPanel2.addComponent(beagleInfoCheckBox);

            optionPanel1.addComponent(optionPanel2);
            optionPanel.addSpanningComponent(optionPanel1);
        } else {
            final JEditorPane beagleInfo = new JEditorPane("text/html",
                    "<html><div style=\"font-family:'helvetica neue light',helvetica,sans-serif;font-size:12;\">" +
                            "<p><b>The BEAGLE library was not found.</b></p>" +
                            "<p>BEAGLE is a high-performance phylogenetic library that is required to<br>" +
                            "run BEAST and can make use of additional computational resources such<br>" +
                            "as graphics boards. It must be downloaded and installed independently<br>" +
                            "of BEAST:</p>" +
                            "<pre><a href=\"http://github.com/beagle-dev/beagle-lib/\">http://github.com/beagle-dev/beagle-lib/</a></pre></div></html>");
            beagleInfo.setOpaque(false);
            beagleInfo.setEditable(false);
            beagleInfo.addHyperlinkListener(new SimpleLinkListener());
            optionPanel.addSpanningComponent(beagleInfo);
        }

        beagleResourceCombo.setSelectedItem("CPU");
    }

    public boolean showDialog(String title) {

        JOptionPane optionPane = new JOptionPane(optionPanel,
                JOptionPane.PLAIN_MESSAGE,
                JOptionPane.OK_CANCEL_OPTION,
                null,
                new String[]{"Run", "Quit"},
                "Run");
        optionPane.setBorder(new EmptyBorder(12, 12, 12, 12));

        final JDialog dialog = optionPane.createDialog(frame, title);
        //dialog.setResizable(true);
        dialog.pack();

        while (true) {
            optionPane.setValue(JOptionPane.UNINITIALIZED_VALUE);
            dialog.setVisible(true);

            if (optionPane.getValue() == null || !optionPane.getValue().equals("Run")) {
                // Quit was pressed or the window was closed
                dialog.dispose();
                return false;
            }

            String problem = checkSettings();
            if (problem == null) {
                dialog.dispose();
                return true;
            }

            // show the problem and return to the settings dialog
            JOptionPane.showMessageDialog(dialog, problem, "BEAST", JOptionPane.ERROR_MESSAGE);
        }
    }

    /**
     * @return a description of what is wrong with the settings, or null if BEAST can be run with them
     */
    private String checkSettings() {
        if (inputFile == null && !beagleInfoCheckBox.isSelected()) {
            return "No BEAST XML file has been selected.";
        }
        if (!outputDirectoryCheckBox.isSelected()) {
            return null;
        }
        if (outputDirectory == null) {
            return "No output folder has been selected.\nChoose a folder or deselect the option to write output files to a folder.";
        }
        try {
            OutputFileHelper.validateOutputDirectory(outputDirectory);
        } catch (IOException ioe) {
            return ioe.getMessage();
        }
        return null;
    }

    /**
     * Ask the user to pick a folder, using the native folder chooser on macOS.
     *
     * @return the chosen folder or null if the dialog was cancelled
     */
    private static File chooseDirectory(JFrame frame, String title, File current) {
        if (OSType.isMac()) {
            System.setProperty("apple.awt.fileDialogForDirectories", "true");
            try {
                FileDialog dialog = new FileDialog(frame, title, FileDialog.LOAD);
                if (current != null) {
                    dialog.setDirectory(current.getPath());
                }
                dialog.setVisible(true);
                if (dialog.getFile() == null) {
                    return null;
                }
                return new File(dialog.getDirectory(), dialog.getFile());
            } finally {
                System.setProperty("apple.awt.fileDialogForDirectories", "false");
            }
        }

        JFileChooser chooser = new JFileChooser(current);
        chooser.setDialogTitle(title);
        chooser.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
        chooser.setAcceptAllFileFilterUsed(false);
        if (chooser.showOpenDialog(frame) != JFileChooser.APPROVE_OPTION) {
            return null;
        }
        return chooser.getSelectedFile();
    }

    public long getSeed() {
        return seedText.getLongValue();
    }

    public void setSeed(long seed) {
        seedText.setValue(seed);
    }

    public boolean allowOverwrite() {
        return overwriteCheckBox.isSelected();
    }

    public void setAllowOverwrite(boolean allowOverwrite) {
        overwriteCheckBox.setSelected(allowOverwrite);
    }

    public boolean preferBeagleGPU() {
        return beagleResourceCombo.getSelectedItem().equals("GPU");
    }

    public boolean preferBeagleCPU() {
        return (beagleResourceCombo.getSelectedItem().equals("CPU"));
    }

    public void setPreferBeagleGPU() {
        beagleResourceCombo.setSelectedItem("GPU");
    }

    public boolean preferBeagleSSE() {
        return beagleSSECheckBox.isSelected();
    }

    public void setPreferBeagleSSE(boolean preferBeagleSSE) {
        beagleSSECheckBox.setSelected(preferBeagleSSE);
    }

    public boolean preferBeagleSingle() {
        return beaglePrecisionCombo.getSelectedItem().equals("Single");
    }

    public boolean preferBeagleDouble() {
        return beaglePrecisionCombo.getSelectedItem().equals("Double");
    }

    public void setPreferBeagleSingle() {
         beaglePrecisionCombo.setSelectedItem("Single");
    }

    public String scalingScheme() {
        return ((String) beagleScalingCombo.getSelectedItem()).toLowerCase();
    }

    public void setScalingScheme(String scalingScheme) {
        beagleScalingCombo.setSelectedItem(scalingScheme);
    }

    public boolean showBeagleInfo() {
        return beagleInfoCheckBox.isSelected();
    }

    public int getThreadPoolSize() {
        if (threadsCombo.getSelectedIndex() == 0) {
            // Automatic
            return -1;
        }
        return (Integer) threadsCombo.getSelectedItem();
    }

    public File getInputFile() {
        return inputFile;
    }

    /**
     * @return the folder that all output files should be written to, or null if the
     * default locations should be used
     */
    public File getOutputDirectory() {
        if (!outputDirectoryCheckBox.isSelected()) {
            return null;
        }
        return outputDirectory;
    }

    /**
     * Pre-select an output folder (for example one given with the -output command line option)
     * and turn the option on.
     */
    public void setOutputDirectory(File directory) {
        outputDirectory = (directory != null ? directory.getAbsoluteFile() : null);
        if (outputDirectory != null) {
            outputDirectoryNameText.setText(outputDirectory.getPath());
            outputDirectoryNameText.setCaretPosition(outputDirectoryNameText.getText().length());
            outputDirectoryNameText.setToolTipText(outputDirectory.getPath());
            outputDirectoryCheckBox.setSelected(true);
        } else {
            outputDirectoryNameText.setText("not selected");
        }
    }
}