/*
 * VectorIcons.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
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
 *
 */

package dr.app.beauti.util;

import javax.swing.*;
import java.awt.*;
import java.awt.geom.Area;
import java.awt.geom.Ellipse2D;
import java.awt.geom.Path2D;
import java.awt.geom.Rectangle2D;

/**
 * The small icons of BEAUti's buttons, drawn with Java2D so that they are sharp at any screen resolution
 * (the former 12 and 14 pixel PNGs were blurred on high-resolution displays).
 *
 * @author Andrew Rambaut
 */
public final class VectorIcons {

    private VectorIcons() {
    }

    /**
     * @return the colour to draw the icon in: the component's foreground, or the label colour where there
     * is no component (the macOS screen menu bar paints menu icons with none), greyed when disabled
     */
    private static Color colourFor(Component c, Color enabledColour) {
        if (c != null && !c.isEnabled()) {
            return UIManager.getColor("Label.disabledForeground");
        }
        if (enabledColour != null) {
            return enabledColour;
        }
        return c != null ? c.getForeground() : UIManager.getColor("Label.foreground");
    }

    /**
     * @return a gear (settings) icon in the component's foreground colour
     */
    public static Icon gear(final int size) {
        return new Icon() {
            public void paintIcon(Component c, Graphics g, int x, int y) {
                Graphics2D g2 = (Graphics2D) g.create();
                try {
                    g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
                    g2.translate(x, y);
                    double centre = size / 2.0;
                    double outer = size * 0.36;       // the body, teeth reach the icon's edge
                    double hole = size * 0.13;
                    Area gear = new Area(new Ellipse2D.Double(centre - outer, centre - outer, 2 * outer, 2 * outer));
                    double toothWidth = size * 0.2;
                    for (int i = 0; i < 8; i++) {
                        Path2D tooth = new Path2D.Double(new Rectangle2D.Double(centre - toothWidth / 2, 0.5, toothWidth, centre));
                        tooth.transform(java.awt.geom.AffineTransform.getRotateInstance(i * Math.PI / 4, centre, centre));
                        gear.add(new Area(tooth));
                    }
                    gear.subtract(new Area(new Ellipse2D.Double(centre - hole, centre - hole, 2 * hole, 2 * hole)));
                    g2.setColor(colourFor(c, null));
                    g2.fill(gear);
                } finally {
                    g2.dispose();
                }
            }

            public int getIconWidth() {
                return size;
            }

            public int getIconHeight() {
                return size;
            }
        };
    }

    /**
     * @return a filled circle with a white arrow pointing right (true) or left
     */
    public static Icon circleArrow(final int size, final Color fill, final boolean right) {
        return new Icon() {
            public void paintIcon(Component c, Graphics g, int x, int y) {
                Graphics2D g2 = (Graphics2D) g.create();
                try {
                    g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
                    g2.setRenderingHint(RenderingHints.KEY_STROKE_CONTROL, RenderingHints.VALUE_STROKE_PURE);
                    g2.translate(x, y);
                    g2.setColor(colourFor(c, fill));
                    g2.fill(new Ellipse2D.Double(0.5, 0.5, size - 1, size - 1));

                    double centre = size / 2.0;
                    double half = size * 0.25;                 // half the arrow's length
                    double head = size * 0.18;
                    double tip = right ? centre + half : centre - half;
                    double tail = right ? centre - half : centre + half;
                    double back = right ? tip - head : tip + head;
                    g2.setColor(Color.WHITE);
                    g2.setStroke(new BasicStroke((float) (size * 0.14), BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND));
                    Path2D arrow = new Path2D.Double();
                    arrow.moveTo(tail, centre);
                    arrow.lineTo(tip, centre);
                    arrow.moveTo(back, centre - head);
                    arrow.lineTo(tip, centre);
                    arrow.lineTo(back, centre + head);
                    g2.draw(arrow);
                } finally {
                    g2.dispose();
                }
            }

            public int getIconWidth() {
                return size;
            }

            public int getIconHeight() {
                return size;
            }
        };
    }
}
