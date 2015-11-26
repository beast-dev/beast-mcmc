/*
 * BranchStyle.java
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
public class BranchStyle {
    public BranchStyle() {
    }

    public BranchStyle(String widthProperty, String colorProperty) {
        this.widthProperty = widthProperty;
        this.colorProperty = colorProperty;
    }

    public String getWidthProperty() {
        return widthProperty;
    }

    public void setWidthProperty(String widthProperty) {
        this.widthProperty = widthProperty;
    }

    public double getWidthPropertyMinimum() {
        return widthPropertyMinimum;
    }

    public void setWidthPropertyMinimum(double widthPropertyMinimum) {
        this.widthPropertyMinimum = widthPropertyMinimum;
    }

    public double getWidthPropertyMaximum() {
        return widthPropertyMaximum;
    }

    public void setWidthPropertyMaximum(double widthPropertyMaximum) {
        this.widthPropertyMaximum = widthPropertyMaximum;
    }

    public double getWidthBase() {
        return widthBase;
    }

    public void setWidthBase(double widthBase) {
        this.widthBase = widthBase;
    }

    public double getWidthScale() {
        return widthScale;
    }

    public void setWidthScale(double widthScale) {
        this.widthScale = widthScale;
    }

    public String getColorProperty() {
        return colorProperty;
    }

    public void setColorProperty(String colorProperty) {
        this.colorProperty = colorProperty;
    }

    public double getColorPropertyMinimum() {
        return colorPropertyMinimum;
    }

    public void setColorPropertyMinimum(double colorPropertyMinimum) {
        this.colorPropertyMinimum = colorPropertyMinimum;
    }

    public double getColorPropertyMaximum() {
        return colorPropertyMaximum;
    }

    public void setColorPropertyMaximum(double colorPropertyMaximum) {
        this.colorPropertyMaximum = colorPropertyMaximum;
    }

    public Color getColorStart() {
        return colorStart;
    }

    public void setColorStart(Color colorStart) {
        this.colorStart = colorStart;
    }

    public Color getColorFinish() {
        return colorFinish;
    }

    public void setColorFinish(Color colorFinish) {
        this.colorFinish = colorFinish;
    }

    private String widthProperty = null;
    private double widthPropertyMinimum = 0.0;
    private double widthPropertyMaximum = 1.0;

    private double widthBase = 4.0;
    private double widthScale = 3.0;

    private String colorProperty = null;
    private double colorPropertyMinimum = 0.0;
    private double colorPropertyMaximum = 1.0;

    private Color colorStart = Color.blue;
    private Color colorFinish = Color.white;
}
