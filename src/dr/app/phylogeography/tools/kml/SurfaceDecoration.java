/*
 * SurfaceDecoration.java
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

package dr.app.phylogeography.tools.kml;

import java.awt.*;

/**
 * @author Andrew Rambaut
 * @author Philippe Lemey
 * @version $Id$
 */
public class SurfaceDecoration {
    public SurfaceDecoration() {
    }

    public SurfaceDecoration(String colorProperty, boolean visible, Color startColor, Color endColor, double opacity) {
        this.colorProperty = colorProperty;
        isVisible = visible;
        this.startColor = startColor;
        this.endColor = endColor;
        this.opacity = opacity;
    }

    public String getColorProperty() {
        return colorProperty;
    }

    public void setColorProperty(String colorProperty) {
        this.colorProperty = colorProperty;
    }

    public boolean isVisible() {
        return isVisible;
    }

    public void setVisible(boolean visible) {
        isVisible = visible;
    }

    public Color getStartColor() {
        return startColor;
    }

    public void setStartColor(Color startColor) {
        this.startColor = startColor;
    }

    public Color getEndColor() {
        return endColor;
    }

    public void setEndColor(Color endColor) {
        this.endColor = endColor;
    }

    public double getOpacity() {
        return opacity;
    }

    public void setOpacity(double opacity) {
        this.opacity = opacity;
    }

    private String colorProperty = "";
    private boolean isVisible = true;
    private Color startColor = Color.green;
    private Color endColor = Color.yellow;
    private double opacity = 0.5;

}