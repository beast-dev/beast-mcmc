/*
 * ViewAlignmentDialog.java
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

package dr.app.beauti.datapanel;

import org.virion.jam.panels.OptionsPanel;

import dr.app.beauti.options.PartitionData;
import dr.evolution.util.Taxon;

import javax.swing.*;
import javax.swing.border.EmptyBorder;

import java.awt.BorderLayout;
import java.awt.event.*;


/**
 * @author Walter Xie
 * @version $Id: ViewAlignmentDialog.java,v 1.5 2009/08/11 13:29:34 rambaut Exp $
 */
public class ViewAlignmentDialog {

    private JFrame frame;
    
     private final JScrollPane scrollPane;
    private final PartitionData partitionData;

    public ViewAlignmentDialog(JFrame frame, PartitionData partitionData) {
        this.frame = frame;
        this.partitionData = partitionData;
                
        
        ViewAligmentPanel panel = new ViewAligmentPanel(partitionData);
                
        scrollPane = new JScrollPane(panel);
        scrollPane.setOpaque(false);
        
        panel.repaint();                
        
    }

    public int showDialog() {

        JOptionPane optionPane = new JOptionPane(scrollPane,
                JOptionPane.PLAIN_MESSAGE,
                JOptionPane.OK_CANCEL_OPTION,
                null,
                null,
                null);
        optionPane.setBorder(new EmptyBorder(12, 12, 12, 12));
        optionPane.setPreferredSize(new java.awt.Dimension(600, 300));

        final JDialog dialog = optionPane.createDialog(frame, "View Alignment Given A Partition Data " + partitionData.getName());
        dialog.pack();

        dialog.setVisible(true);

        int result = JOptionPane.CANCEL_OPTION;
        Integer value = (Integer) optionPane.getValue();
        if (value != null && value != -1) {
            result = value;
        }

        return result;
    }


}