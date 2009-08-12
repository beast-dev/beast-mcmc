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
	
import javax.swing.*;

import dr.app.beauti.options.PartitionData;
import dr.evolution.util.Taxon;

import java.awt.*;
import java.awt.event.*;

import java.awt.font.FontRenderContext;
import java.awt.font.LineBreakMeasurer;
import java.awt.font.TextAttribute;
import java.awt.font.TextLayout;
import java.text.AttributedCharacterIterator;
import java.text.AttributedString;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @author Walter Xie
 * @version $Id: ViewAlignmentDialog.java,v 1.5 2009/08/11 13:29:34 rambaut Exp
 *          $
 */

public class ViewAligmentPanel extends JPanel {
	
	private final PartitionData partitionData;
	 
	public ViewAligmentPanel (PartitionData partitionData) {
		this.partitionData = partitionData;
	}

	public void paintComponent(Graphics g) {
		super.paintComponent(g);
		
		setBackground(Color.white);


//		Graphics2D g2d = (Graphics2D) g;
		
		g.setFont(new Font("SansSerif", Font.PLAIN, 14));
		
		for (int i = 0; i < partitionData.getAlignment().getTaxonCount(); i++) {
            Taxon taxon = partitionData.getAlignment().getTaxon(i);
            
            
            g.drawString(taxon.getId(), 0, i * 18);
            
            g.drawString(partitionData.getAlignment().getAlignedSequenceString(i), 300, i * 18);
        }
		
		
		
	}

}
