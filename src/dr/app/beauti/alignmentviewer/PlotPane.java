/*
 * PlotPane.java
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
 * @version $Id: RulerPane.java,v 1.1 2005/11/01 23:52:04 rambaut Exp $
 */
public abstract class PlotPane extends JComponent {
    private AlignmentBuffer alignment = null;
    private ColumnDecorator columnDecorator = null;
    private int columnCount = 0;

    private int colWidth = 0;

    public int getColumnWidth() {

        return colWidth;
    }

    public int getColumnCount() {
        return columnCount;
    }

    public PlotPane() {
        setFont( new Font("sansserif", Font.PLAIN, 10));
    }

    public void setAlignmentBuffer(AlignmentBuffer alignment) {
        this.alignment = alignment;

        if (alignment != null) {
            columnCount = alignment.getSiteCount();
        } else {
            columnCount = 0;
        }

        setup();
    }

    public void setColumnWidth(int columnWidth) {
        this.colWidth = columnWidth;
        setup();
    }

    public void setColumnDecorator(ColumnDecorator columnDecorator) {
        this.columnDecorator = columnDecorator;
    }

    private void setup() {
        setPreferredSize(new Dimension(colWidth * columnCount, 16));
        invalidate();
    }

    public void paint(Graphics graphics) {
        if (alignment == null) return;

        Graphics2D g2 = (Graphics2D)graphics;

        int col1 = 0;
        int col2 = columnCount - 1;

        int y1 = 0;
        int y2 = getHeight();

        if (columnDecorator != null) {
            int x = 0;
            for (int col = col1; col <= col2; col++) {
                Paint columnBackground = columnDecorator.getColumnBackground(col);
                if (columnBackground != null) {
                    g2.setPaint(columnBackground);
                    g2.fillRect(x, y1, colWidth, y2);
                }
                x += colWidth;
            }
        }

		paintPlot(g2);
    }

	public abstract void paintPlot(Graphics2D g2);

}