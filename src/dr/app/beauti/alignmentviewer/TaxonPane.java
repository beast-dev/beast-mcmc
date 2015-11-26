/*
 * TaxonPane.java
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

package dr.app.beauti.alignmentviewer;

import javax.swing.*;
import java.awt.*;

/**
 * @author Andrew Rambaut
 * @version $Id: TaxonPane.java,v 1.1 2005/11/01 23:52:04 rambaut Exp $
 */
public class TaxonPane extends JPanel {

	private JList taxonList;
	private DefaultListModel taxonListModel;

	public TaxonPane() {

		setLayout(new BorderLayout());

		taxonListModel = new DefaultListModel();
		taxonList = new JList(taxonListModel);
		taxonList.setFont(new Font("sansserif", Font.PLAIN, 10));

		add(taxonList, BorderLayout.CENTER);
	}

	public void setAlignmentBuffer(AlignmentBuffer alignment) {
		taxonListModel.removeAllElements();
		if (alignment != null) {
			for (int i = 0; i < alignment.getSequenceCount(); i++) {
				taxonListModel.addElement(alignment.getTaxonLabel(i));
			}
		}
	}

	public void setRowHeight(int rowHeight) {
		taxonList.setFixedCellHeight(rowHeight);
	}

}
