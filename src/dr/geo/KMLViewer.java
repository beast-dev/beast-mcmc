/*
 * KMLViewer.java
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

package dr.geo;

import javax.swing.*;
import java.awt.*;
import java.awt.geom.AffineTransform;
import java.awt.geom.Line2D;
import java.awt.geom.Rectangle2D;
import java.awt.image.BufferedImage;

/**
 * @author Alexei Drummond
 */
public class KMLViewer extends JComponent {

    KMLRenderer renderer;
    BufferedImage image;

    public KMLViewer(String kmlFileName) {

        renderer = new KMLRenderer(kmlFileName, Color.green, Color.blue);
    }

    public void paintComponent(Graphics g) {

        if (image == null || image.getWidth() != getWidth() || image.getHeight() != getHeight()) {
            image = new BufferedImage(getWidth(), getHeight(), BufferedImage.TYPE_INT_ARGB);
            renderer.render(image);
        }

        g.drawImage(image, 0, 0, null);

        drawGrid(5, 5, (Graphics2D) g, renderer.viewTransform.getTransform());
    }

    private void drawGrid(int dLat, int dLong, Graphics2D g2d, AffineTransform transform) {

        for (double longitude = -180; longitude < 180; longitude += dLong) {
            Line2D line = new Line2D.Double(longitude, -90, longitude, 90);
            g2d.draw(transform.createTransformedShape(line));
        }
        for (double lat = -90; lat < 90; lat += dLat) {
            Line2D line = new Line2D.Double(-180, lat, 180, lat);
            g2d.draw(transform.createTransformedShape(line));
        }
    }

    public static void main(String[] args) {

        String filename = args[0];

        JFrame frame = new JFrame("KMLViewer - " + filename);

        KMLViewer viewer = new KMLViewer(filename);
        Rectangle2D viewport = viewer.renderer.getBounds();

        frame.getContentPane().add(BorderLayout.CENTER, viewer);

        int width;
        int height;
        if (viewport.getHeight() > viewport.getWidth()) {
            height = 900;
            width = (int) (height * viewport.getWidth() / viewport.getHeight());
        } else {
            width = 900;
            height = (int) (width * viewport.getHeight() / viewport.getWidth());
        }
        System.out.println("Height = " + height + ", Width = " + width);

        frame.setSize(width, height);
        frame.setVisible(true);
    }
}