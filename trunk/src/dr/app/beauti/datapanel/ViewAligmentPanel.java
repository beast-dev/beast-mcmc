/*
 * ViewAligmentPanel.java
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
	
import javax.swing.*;

import dr.app.beauti.options.PartitionData;
import dr.evolution.util.Taxon;

import java.awt.*;
import java.awt.event.*;


/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @author Walter Xie
 * @version $Id: ViewAligmentPanel.java,v 1.5 2009/08/11 13:29:34 rambaut Exp$
 */

public class ViewAligmentPanel extends JPanel {
		
	private static final long serialVersionUID = 6129057786946772468L;
	
	private final PartitionData partitionData;
	private final static int FONT_SIZE = 14;
	private final static int GAP = 4;
	private final static int TAXON_NAME_LEN = 300;
	private final static int Y_BOADER = FONT_SIZE + 3 * GAP + 20;
	private int width = TAXON_NAME_LEN;
	private int height = Y_BOADER;
	
	 
	public ViewAligmentPanel (PartitionData partitionData) {
		this.partitionData = partitionData;
		
		getWidthLenth();
	}

	public void paintComponent(Graphics g) {
		super.paintComponent(g);		
		setBackground(Color.white);

//		Graphics2D g2d = (Graphics2D) g;
		
		g.setFont(new Font("SansSerif", Font.PLAIN, FONT_SIZE));
		
		String sequenceString = partitionData.getAlignment().getAlignedSequenceString(0);
        for (int i = 0; i < sequenceString.length(); i++) {
        	if (i == 0 || ((i+1)%10) == 0) {
        		if (i == 0) {
        			g.drawString(Integer.toString(i+1), TAXON_NAME_LEN + i * (FONT_SIZE + GAP), Y_BOADER - (FONT_SIZE + 3 * GAP)); 
        		} else if (i < 100) {
        			g.drawString(Integer.toString(i+1), TAXON_NAME_LEN + i * (FONT_SIZE + GAP) - GAP, Y_BOADER - (FONT_SIZE + 3 * GAP)); 
        		} else {
        			g.drawString(Integer.toString(i+1), TAXON_NAME_LEN + i * (FONT_SIZE + GAP) -  2 * GAP, Y_BOADER - (FONT_SIZE + 3 * GAP)); 
        		} 
        		
        		g.drawLine(TAXON_NAME_LEN + i * (FONT_SIZE + GAP) + GAP, Y_BOADER - (FONT_SIZE + 2 * GAP), 
        				TAXON_NAME_LEN + i * (FONT_SIZE + GAP) + GAP, Y_BOADER - (FONT_SIZE + GAP)); 
        	}
        }
		
		
		for (int i = 0; i < partitionData.getAlignment().getTaxonCount(); i++) {
            
			Taxon taxon = partitionData.getAlignment().getTaxon(i);        
            
            g.drawString(taxon.getId(), 0, Y_BOADER + i * (FONT_SIZE + GAP));
            
            sequenceString = partitionData.getAlignment().getAlignedSequenceString(i);
            
            int charSpace = TAXON_NAME_LEN;            	
            for (char c : sequenceString.toCharArray()) {
            	g.drawString(Character.toString(c), charSpace, Y_BOADER + i * (FONT_SIZE + GAP)); 
            	charSpace = charSpace + FONT_SIZE + GAP;
            }
            
        }
		
		
		
	}
		
	public void setPreferredSize(){			
		super.setPreferredSize(new java.awt.Dimension(width, height));
//		System.out.println(width + "     " + height);
	}
	
	
	private void getWidthLenth() {
		int maxLen = 0;
		
		String sequenceString = null;
		for (int i = 0; i < partitionData.getAlignment().getTaxonCount(); i++) {
            
            sequenceString = partitionData.getAlignment().getAlignedSequenceString(i);
            
            if (sequenceString.length() > maxLen) maxLen = sequenceString.length();
                   
            height = height + FONT_SIZE + GAP;
        }
		
		width = width + maxLen * (FONT_SIZE + GAP);  		
	}
}
