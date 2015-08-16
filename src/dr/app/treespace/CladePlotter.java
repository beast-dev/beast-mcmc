/*
 * CladePlotter.java
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

package dr.app.treespace;

import javax.swing.*;
import java.awt.*;
import java.awt.geom.Ellipse2D;
import java.awt.geom.Line2D;
import java.util.ArrayList;
import java.util.List;
import dr.app.treespace.CladeSystem.Clade;

/**
 * @author Andrew Rambaut
 * @version $Id$
 */
public class CladePlotter extends JComponent {

    public CladePlotter(final CladeSystem cladeSystem) {
        setCladeSystem(cladeSystem);
    }

    public void setCladeSystem(final CladeSystem cladeSystem) {
        this.cladeSystem = cladeSystem;

        this.clades = cladeSystem.getClades();

        List<Clade> tips = new ArrayList<Clade>();
        for (Clade clade : clades) {
            if (clade.bits.cardinality() == 1) {
                tips.add(clade);
            }
        }

        for (Clade clade : clades) {
            if (clade.bits.cardinality() == tips.size()) {
                rootClade = clade;
            }
        }

        isCalibrated = false;
        repaint();
    }

    public void paint(Graphics g) {
        if (rootClade == null) {
            return;
        }

        if (!isCalibrated) {
            double[] yPosition = { 0.0 };

            calibrateClade((Graphics2D)g, rootClade, yPosition);

            scaleX = ((double)getWidth()) / rootClade.x;
            scaleY = ((double)getHeight()) / yPosition[0];

            isCalibrated = true;
        }

        paintClade((Graphics2D)g, rootClade);
    }

    private void paintClade(Graphics2D g2, Clade clade) {
        if (clade.children != null) {
            for (Clade child : clade.children.keySet()) {
                paintClade(g2, child);
                g2.setColor(Color.black);
                g2.setStroke(new BasicStroke(1.0F));
                g2.draw(new Line2D.Double(transformX(child.x), transformY(child.y), transformX(clade.x), transformY(clade.y)));
            }
        } else {
            paintTip(g2, clade);
        }

        float x = transformX(clade.x);
        float y = transformY(clade.y);

        Shape node = new Ellipse2D.Float(x, y, (float)(nodeWidth * scaleX), (float)(nodeHeight * scaleY));
        g2.setColor(Color.yellow);
        g2.fill(node);
        g2.setColor(Color.black);
        g2.draw(node);
    }

    private void paintTip(Graphics2D g2, Clade clade) {
        float x = transformX(clade.x);
        float y = transformY(clade.y);

        Shape node = new Ellipse2D.Float(x, y, (float)(nodeWidth * scaleX), (float)(nodeHeight * scaleY));
        g2.setColor(Color.yellow);
        g2.fill(node);
        g2.setColor(Color.black);
        g2.draw(node);
        if (clade.label != null) {
            g2.drawString(clade.label, x + xLabelOffset, yLabelOffset);
        }
    }

    private double calibrateClade(Graphics2D g2, Clade clade, double[] yPosition) {
        double y;
        if (clade.children != null) {
            y = 0.0;
            for (Clade child : clade.children.keySet()) {
                y += calibrateClade(g2, child, yPosition);
            }

            // todo: use an average weighted by support
            y /= clade.children.keySet().size();
        } else {
            y = yPosition[0];
            yPosition[0]++;
        }
        clade.x = (float)clade.bits.cardinality();
        clade.y = (float)y;

        return y;
    }

    private float transformX(double x) {
        return (float)getWidth() - (float)(x * scaleX);
    }

    private float transformY(double y) {
        return (float)(y * scaleY);
    }

    CladeSystem cladeSystem;
    Clade rootClade;
    List<Clade> clades;

    private double nodeWidth = 5.0;
    private double nodeHeight = 5.0;

    private float xLabelOffset = 5.0F;
    private float yLabelOffset = 5.0F;
    private boolean isCalibrated = false;

    private double scaleX = 1.0;
    private double scaleY = 1.0;

}
