/*
 * RootSequenceEditor.java
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

package dr.app.bss;

import jam.panels.OptionsPanel;

import java.awt.BorderLayout;
import java.awt.Dimension;
import java.awt.Frame;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTextField;
import javax.swing.SwingConstants;
import javax.swing.SwingUtilities;

/**
 * @author Filip Bielejec
 * @version $Id$
 */
public class RootSequenceEditor {

    // Data
    private PartitionDataList dataList = null;
    private int row;

    // Settings
    private OptionsPanel optionPanel;
    //TODO: class that checks for ACTG string values, maybe also with scroller
    private JTextField ancestralSequenceField;

    //Buttons
    private JButton done;
    private JButton cancel;

    // Window
    private JDialog window;
    private Frame owner;

    public RootSequenceEditor(PartitionDataList dataList, int row) {

        this.dataList = dataList;
        this.row = row;

        ancestralSequenceField = new JTextField("", 10);
        window = new JDialog(owner, "Setup root sequence for partition " + (row + 1));
        optionPanel = new OptionsPanel(12, 12, SwingConstants.CENTER);

        optionPanel.setOpaque(false);

        setAncestralSequence();

        // Buttons
        JPanel buttonsHolder = new JPanel();
        buttonsHolder.setOpaque(false);

        cancel = new JButton("Cancel", Utils.createImageIcon(Utils.CLOSE_ICON));
        cancel.addActionListener(new ListenCancel());
        buttonsHolder.add(cancel);

        done = new JButton("Done", Utils.createImageIcon(Utils.CHECK_ICON));
        done.addActionListener(new ListenOk());
        buttonsHolder.add(done);

        // Window
        owner = Utils.getActiveFrame();
        window.setLocationRelativeTo(owner);
        window.getContentPane().setLayout(new BorderLayout());
        window.getContentPane().add(optionPanel, BorderLayout.CENTER);
        window.getContentPane().add(buttonsHolder, BorderLayout.SOUTH);
        window.pack();

    }//END: Constructor

    private void setAncestralSequence() {

        optionPanel.removeAll();

        ancestralSequenceField.setText(dataList.get(row).ancestralSequenceString);
        optionPanel.addComponents(new JLabel("Root sequence:"),
                ancestralSequenceField);

        window.validate();
        window.repaint();
    }// END: setAncestralSequence

    public void collectSettings() {

        dataList.get(row).ancestralSequenceString = ancestralSequenceField.getText();

    }// END: collectSettings

    private class ListenOk implements ActionListener {
        public void actionPerformed(ActionEvent ev) {

            int ancestralSequenceLength = ancestralSequenceField.getText()
                    .length();
            int partitionSiteCount = dataList.get(row)
                    .createPartitionSiteCount();

            if (ancestralSequenceLength == partitionSiteCount) {

                if (ancestralSequenceField.getText().contains("-")) {

                    Utils.showDialog("Ancestral sequence contains one or more gaps.");

                } else {

                    window.setVisible(false);
                    collectSettings();

                }

            } else {

                if (ancestralSequenceLength != 0) {

                    Utils.showDialog("Ancestral sequence length of "
                            + ancestralSequenceLength
                            + " does not match partition site count of "
                            + partitionSiteCount + ".");
                } else {

                    window.setVisible(false);

                }
            }

        }// END: actionPerformed
    }// END: ListenSaveLocationCoordinates

    private class ListenCancel implements ActionListener {
        public void actionPerformed(ActionEvent ev) {

            window.setVisible(false);

        }// END: actionPerformed
    }// END: ListenCancel

    public void showWindow() {
        window.setDefaultCloseOperation(JFrame.HIDE_ON_CLOSE);
        window.setSize(new Dimension(450, 200));
        window.setMinimumSize(new Dimension(100, 100));
        window.setResizable(true);
        window.setModal(true);
        window.setVisible(true);
    }// END: showWindow

    public void launch() {

        if (SwingUtilities.isEventDispatchThread()) {
            showWindow();
        } else {
            SwingUtilities.invokeLater(new Runnable() {
                public void run() {
                    showWindow();
                }
            });
        }// END: edt check

    }// END: launch

}//END: class
