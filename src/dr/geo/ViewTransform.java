/*
 * ViewTransform.java
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

import java.awt.geom.AffineTransform;
import java.awt.geom.Point2D;
import java.awt.geom.Rectangle2D;

/**
 * @author Alexei Drummond
 */
public class ViewTransform {

    double scaleX;
    double scaleY;

    int width;
    int height;

    Point2D topLeft;
    Point2D bottomRight;

    AffineTransform transform;

    public ViewTransform(Rectangle2D bounds, int width, int height) {

        set(bounds, width, height);
    }

    public AffineTransform getTransform() {
        return transform;
    }

    public Rectangle2D getViewport() {
        return new Rectangle2D.Double(topLeft.getX(), topLeft.getY(),
                bottomRight.getX() - topLeft.getX(), bottomRight.getY() - topLeft.getY());
    }


    public void set(Rectangle2D bounds, int width, int height) {
        topLeft = new Point2D.Double(bounds.getMinX(), bounds.getMinY());
        bottomRight = new Point2D.Double(bounds.getMaxX(), bounds.getMaxY());
        this.width = width;
        this.height = height;
        computeScalesAndTransform();
    }

    public void setBounds(Rectangle2D bounds) {
        topLeft = new Point2D.Double(bounds.getMinX(), bounds.getMinY());
        bottomRight = new Point2D.Double(bounds.getMaxX(), bounds.getMaxY());
        computeScalesAndTransform();
    }

    public void setSize(int width, int height) {
        this.width = width;
        this.height = height;
        computeScalesAndTransform();
    }

    private void computeScalesAndTransform() {
        scaleX = width / (bottomRight.getX() - topLeft.getX());
        scaleY = height / (bottomRight.getY() - topLeft.getY());

        transform = computeFullTransform();
    }

    private AffineTransform getTranslate() {
        return AffineTransform.getTranslateInstance(-topLeft.getX(), -bottomRight.getY());
    }

    private AffineTransform getScale() {
        return AffineTransform.getScaleInstance(scaleX, -scaleY);
    }

    private AffineTransform computeFullTransform() {
        AffineTransform transform = getScale();
        transform.concatenate(getTranslate());
        return transform;
    }
}
