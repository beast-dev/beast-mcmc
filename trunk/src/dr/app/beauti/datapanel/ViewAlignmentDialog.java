/*
 * PriorDialog.java
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
import java.awt.event.*;


/**
 * @author Walter Xie
 * @version $Id: ViewAlignmentDialog.java,v 1.5 2009/08/11 13:29:34 rambaut Exp $
 */
public class ViewAlignmentDialog {

    private JFrame frame;
    
    private final OptionsPanel optionPanel;
    

    public ViewAlignmentDialog(JFrame frame, PartitionData partitionData) {
        this.frame = frame;
                
        optionPanel = new OptionsPanel(12, 12);
        
        JTextField partitionDataNameField = new JTextField(24);
//        if (partitionData.getName().length() > 23) partitionDataNameField.setColumns(partitionData.getName().length() + 5);
        partitionDataNameField.setText(partitionData.getName());
        optionPanel.addComponentWithLabel("Partition data name:", partitionDataNameField);
        
        JTextField sequenceTypeField = new JTextField(16);
        sequenceTypeField.setText(partitionData.getAlignment().getDataType().getDescription());
        optionPanel.addComponentWithLabel("Sequence Type:", sequenceTypeField);    

//        JTextField alignmentNameField = new JTextField(12);
//        alignmentNameField.setText(partitionData.getAlignment().getId());
//        optionPanel.addComponentWithLabel("Alignment name:", alignmentNameField);
        
        JTextField alignmentPartternCountNameField = new JTextField(8);
        alignmentPartternCountNameField.setText(Integer.toString(partitionData.getAlignment().getSiteCount()));
        optionPanel.addComponentWithLabel("Number of sites in this alignment:", alignmentPartternCountNameField);
        
        if (partitionData.getFromSite() > 0) {
	        JTextField partitionDataFromField = new JTextField(8);
	        partitionDataFromField.setText(Integer.toString(partitionData.getFromSite()));
	        optionPanel.addComponentWithLabel("Site starts from", partitionDataFromField);
	        
	        JTextField partitionDataToField = new JTextField(8);
	        partitionDataToField.setText(Integer.toString(partitionData.getToSite()));
	        optionPanel.addComponentWithLabel("to", partitionDataToField);
        }
        
        JTextField alignmentTaxonCountNameField = new JTextField(8);
        alignmentTaxonCountNameField.setText(Integer.toString(partitionData.getAlignment().getTaxonCount()));
        optionPanel.addComponentWithLabel("Number of taxon in this alignment:", alignmentTaxonCountNameField);
        
        JComboBox taxonListCombo = new JComboBox();
        for (int i = 0; i < partitionData.getAlignment().getTaxonCount(); i++) {
            Taxon taxon = partitionData.getAlignment().getTaxon(i);
            taxonListCombo.addItem(taxon.getId());
        }
        optionPanel.addComponentWithLabel("Select a taxon in this alignment:", taxonListCombo);
 
    }

    public int showDialog() {

        JOptionPane optionPane = new JOptionPane(optionPanel,
                JOptionPane.INFORMATION_MESSAGE,
                JOptionPane.OK_CANCEL_OPTION,
                null,
                null,
                null);
        optionPane.setBorder(new EmptyBorder(12, 12, 12, 12));

        final JDialog dialog = optionPane.createDialog(frame, "View Alignment Given A Partition Data");
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