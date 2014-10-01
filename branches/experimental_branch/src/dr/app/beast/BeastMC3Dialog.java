/*
 * BeastMC3Dialog.java
 *
 * Copyright (C) 2002-2006 Alexei Drummond and Andrew Rambaut
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
import dr.app.gui.components.RealNumberField;
import dr.app.gui.components.WholeNumberField;
import jam.panels.OptionsPanel;

import javax.swing.*;
import javax.swing.border.Border;
import javax.swing.border.EmptyBorder;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyAdapter;
import java.awt.event.KeyEvent;
import java.io.File;
import java.text.DecimalFormat;

/**
 * MCMCMC functionality now rolled into the BEAST Main class.
 * @deprecated
 */
public class BeastMC3Dialog {
    private final JFrame frame;

    private final OptionsPanel optionPanel;

    private final WholeNumberField chainsText = new WholeNumberField((int) 2, 100);
    private final JTextField temperaturesArrayText = new JTextField("-1");
    private final RealNumberField deltaText = new RealNumberField((double) 0.0, Double.MAX_VALUE);
    private final WholeNumberField swapText = new WholeNumberField((int) 1, Integer.MAX_VALUE);
    private final JCheckBox overwriteCheckBox = new JCheckBox("Allow overwriting of log files");
    private final JRadioButton deltaButton = new JRadioButton("Delta: ");
    private final JRadioButton temperaturesButton = new JRadioButton("Temperatures: ");

    private File inputFile = null;

    public BeastMC3Dialog(final JFrame frame, final String titleString, final Icon icon) {
        this.frame = frame;

        optionPanel = new OptionsPanel(12, 12);

        //this.frame = frame;

        JPanel panel = new JPanel(new BorderLayout());
        panel.setOpaque(false);

        final JLabel titleText = new JLabel(titleString);
        titleText.setIcon(icon);
        optionPanel.addSpanningComponent(titleText);
        titleText.setFont(new Font("sans-serif", 0, 12));

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
        optionPanel.addComponentWithLabel("BEAST XML File: ", panel1);

        Color focusColor = UIManager.getColor("Focus.color");
        Border focusBorder = BorderFactory.createMatteBorder(2, 2, 2, 2, focusColor);
        new FileDrop(null, inputFileNameText, focusBorder, new FileDrop.Listener() {
            public void filesDropped(File[] files) {
                inputFile = files[0];
                inputFileNameText.setText(inputFile.getName());
            }   // end filesDropped
        }); // end FileDrop.Listener

        optionPanel.addComponent(overwriteCheckBox);

        optionPanel.addSeparator();

        ActionListener buttonGroupListener = new ActionListener() {
            public void actionPerformed(ActionEvent ae) {
                deltaText.setEnabled(ae.getActionCommand().startsWith("D"));
                temperaturesArrayText.setEnabled(!ae.getActionCommand().startsWith("D"));
            }
        };

        String MESSAGE_DELTA_TEMP = "Either delta or temperatures option should be used, not both";
        deltaButton.setActionCommand("D");
        deltaButton.addActionListener(buttonGroupListener);
        deltaButton.setToolTipText(MESSAGE_DELTA_TEMP);

        temperaturesButton.addActionListener(buttonGroupListener);
        temperaturesButton.setActionCommand("T");
        temperaturesButton.setToolTipText(MESSAGE_DELTA_TEMP);

        ButtonGroup group = new ButtonGroup();
        group.add(deltaButton);
        group.add(temperaturesButton);
        deltaButton.setSelected(true);

        chainsText.setColumns(6);
        optionPanel.addComponentWithLabel("Number of chains: ", chainsText);
        chainsText.setToolTipText("number of chains");

        chainsText.addKeyListener(new KeyAdapter() {
            public void keyReleased(KeyEvent e) {
                temperaturesArrayText.setText(getDefaultTemperaturesArrayString());
            }
        });

        optionPanel.addLabel(MESSAGE_DELTA_TEMP + " :");

        deltaText.setColumns(6);
//        optionPanel.addComponentWithLabel("Delta: ", deltaText);
        deltaText.setToolTipText("temperature increment parameter");
        optionPanel.addComponents(deltaButton, deltaText);

        temperaturesArrayText.setColumns(20);
//        optionPanel.addComponentWithLabel("Temperatures: ", temperaturesArrayText);
        temperaturesArrayText.setToolTipText("a comma-separated list of the hot chain temperatures");
        optionPanel.addComponents(temperaturesButton, temperaturesArrayText);

        swapText.setColumns(10);
        optionPanel.addComponentWithLabel("Swap chain every: ", swapText);
        swapText.setToolTipText("frequency at which chains temperatures will be swapped");

        chainsText.setValue(BeastMC3.HOT_CHAIN_COUNT + 1);
        swapText.setValue(BeastMC3.SWAP_CHAIN_EVERY);

        deltaText.setValue(BeastMC3.DELTA);
        temperaturesArrayText.setText(getDefaultTemperaturesArrayString());
        temperaturesArrayText.setEnabled(false); // choose delta then disable temp...
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

        return optionPane.getValue().equals("Run");
    }

    public boolean allowOverwrite() {
        return overwriteCheckBox.isSelected();
    }

    public int getChains() {
        return chainsText.getValue();
    }

    public double getDelta() {
        return deltaText.getValue();
    }

    public double[] getTemperaturesArray() {
        double[] chainTemperatures = new double[getChains()];
        chainTemperatures[0] = 1.0;

        if (temperaturesButton.isSelected()) {
            String line = temperaturesArrayText.getText();
            if (line != null) {
                String[] tl = line.split(",\\s*");
                double[] hotChainTemperatures = new double[tl.length];
                assert hotChainTemperatures.length == getChains() - 1;

                for (int i = 0; i < tl.length; i++) {
                    hotChainTemperatures[i] = Double.parseDouble(tl[i]);
                }

                System.arraycopy(hotChainTemperatures, 0, chainTemperatures, 1, getChains() - 1);
            } else {
                throw new IllegalArgumentException("Cannot parse hot chain temperatures list (" + line + ")!\n" +
                        "Please use a comma-separated list (e.g. 0.5,0.3)");
            }

        } else {
            for (int i = 1; i < getChains(); i++) {
                chainTemperatures[i] = 1.0 / (1.0 + (getDelta() * i));
            }
        }
        return chainTemperatures;
    }

    public int getSwap() {
        return swapText.getValue();
    }

    public File getInputFile() {
        return inputFile;
    }

    private String getDefaultTemperaturesArrayString() { // exclude chainTemperatures[0] = 1.0
        DecimalFormat df = new DecimalFormat("#0.###");

        String defaultTA = "";
        for (int i = 1; i < getChains(); i++) {
            if (i == 1) {
                defaultTA = df.format(1.0 / (1.0 + (getDelta() * 1)));
            } else {
                defaultTA += "," + df.format(1.0 / (1.0 + (getDelta() * i)));
            }
        }
        return defaultTA;
    }
}