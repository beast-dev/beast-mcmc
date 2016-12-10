/*
 * BeastDialog.java
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

package dr.app.beast;

import dr.app.gui.FileDrop;
import dr.app.gui.components.WholeNumberField;
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
import java.io.File;


public class BeastDialog {
    private final JFrame frame;

    private final OptionsPanel optionPanel;

    private final WholeNumberField seedText = new WholeNumberField((long) 1, Long.MAX_VALUE);
    private final JCheckBox overwriteCheckBox = new JCheckBox("Allow overwriting of log files");
    private final JCheckBox beagleCheckBox = new JCheckBox("Use BEAGLE library if available:");
    private final JCheckBox beagleInfoCheckBox = new JCheckBox("Show list of available BEAGLE resources and Quit");
    private final JComboBox beagleResourceCombo = new JComboBox(new Object[]{"CPU", "GPU"});
    private final JCheckBox beagleSSECheckBox = new JCheckBox("Use CPU's SSE extensions when possible");
    private final JComboBox beaglePrecisionCombo = new JComboBox(new Object[]{"Double", "Single"});
    private final JComboBox beagleScalingCombo = new JComboBox(new Object[]{"Default", "Dynamic", "Delayed", "Always", "Never"});

    private final JComboBox threadsCombo = new JComboBox(new Object[]{"Automatic", 0, 1, 2, 3, 4, 5, 6, 7, 8});

    private File inputFile = null;

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

        optionPanel.addSpanningComponent(beagleCheckBox);
        beagleCheckBox.setSelected(true);

        final OptionsPanel optionPanel1 = new OptionsPanel(0, 6);
//        optionPanel1.setBorder(BorderFactory.createEmptyBorder());
        optionPanel1.setBorder(new TitledBorder(""));

        OptionsPanel optionPanel2 = new OptionsPanel(0, 3);
        optionPanel2.setBorder(BorderFactory.createEmptyBorder());
        final JLabel label1 = optionPanel2.addComponentWithLabel("Prefer use of: ", beagleResourceCombo);
        optionPanel2.addComponent(beagleSSECheckBox);
        beagleSSECheckBox.setSelected(true);
        final JLabel label2 = optionPanel2.addComponentWithLabel("Prefer precision: ", beaglePrecisionCombo);
        final JLabel label3 = optionPanel2.addComponentWithLabel("Rescaling scheme: ", beagleScalingCombo);
        optionPanel2.addComponent(beagleInfoCheckBox);
        optionPanel2.setBorder(BorderFactory.createEmptyBorder());

        optionPanel1.addComponent(optionPanel2);

        final JEditorPane beagleInfo = new JEditorPane("text/html",
                "<html><div style=\"font-family:'helvetica neue light',helvetica,sans-serif;font-size:12;\"><p>BEAGLE is a high-performance phylogenetic library that can make use of<br>" +
                        "additional computational resources such as graphics boards. It must be<br>" +
                        "downloaded and installed independently of BEAST:</p>" +
                        "<pre><a href=\"http://github.com/beagle-dev/beagle-lib/\">http://github.com/beagle-dev/beagle-lib/</a></pre></div></html>");
        beagleInfo.setOpaque(false);
        beagleInfo.setEditable(false);
        beagleInfo.addHyperlinkListener(new SimpleLinkListener());
        optionPanel1.addComponent(beagleInfo);
        optionPanel1.setBorder(BorderFactory.createEmptyBorder());
        optionPanel.addSpanningComponent(optionPanel1);

        beagleInfoCheckBox.setEnabled(false);
        beagleCheckBox.addChangeListener(new ChangeListener() {
            public void stateChanged(ChangeEvent e) {
                beagleInfo.setEnabled(beagleCheckBox.isSelected());
                beagleInfoCheckBox.setEnabled(beagleCheckBox.isSelected());
                label1.setEnabled(beagleCheckBox.isSelected());
                beagleResourceCombo.setEnabled(beagleCheckBox.isSelected());
                beagleSSECheckBox.setEnabled(beagleCheckBox.isSelected());
                label2.setEnabled(beagleCheckBox.isSelected());
                beaglePrecisionCombo.setEnabled(beagleCheckBox.isSelected());
                label3.setEnabled(beagleCheckBox.isSelected());
                beagleScalingCombo.setEnabled(beagleCheckBox.isSelected());
            }
        });

        beagleCheckBox.setSelected(false);
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

        dialog.setVisible(true);

        return (optionPane.getValue() != null ? optionPane.getValue().equals("Run") : false);
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

    public boolean useBeagle() {
        return beagleCheckBox.isSelected();
    }

    public void setUseBeagle(boolean useBeagle) {
         beagleCheckBox.setSelected(useBeagle);
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
}