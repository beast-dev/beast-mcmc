/*
 * ChartLayout.java
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

package dr.app.gui.chart;

import java.awt.*;

/**
 * A layout manager similar to BorderLayout but with specific regions for
 * parts of a chart (plot, axis, titles, etc.).
 *
 * @author Andrew Rambaut
 */

public class ChartLayout implements LayoutManager2 {
    int hgap;
    int vgap;

    Component title = null;
    Component xLabel = null;
    Component yLabel = null;
    Component chart = null;

    /**
     * Constructs a new ChartLayout.
     */
    public ChartLayout() {
        this(0, 0);
    }

    /**
     * Constructs a ChartLayout with the specified gaps.
     *
     * @param hgap the horizontal gap
     * @param vgap the vertical gap
     */
    public ChartLayout(int hgap, int vgap) {
        this.hgap = hgap;
        this.vgap = vgap;
    }

    /**
     * @return the horizontal gap between components.
     */
    public int getHgap() {
        return hgap;
    }

    /**
     * Sets the horizontal gap between components.
     *
     * @param hgap the horizontal gap between components
     */
    public void setHgap(int hgap) {
        this.hgap = hgap;
    }

    /**
     * @return the vertical gap between components.
     */
    public int getVgap() {
        return vgap;
    }

    /**
     * Sets the vertical gap between components.
     *
     * @param vgap the vertical gap between components
     */
    public void setVgap(int vgap) {
        this.vgap = vgap;
    }

    /**
     * Adds the specified named component to the layout.
     *
     * @param comp the component to be added
     */
    public void addLayoutComponent(Component comp, Object constraints) {
        synchronized (comp.getTreeLock()) {
            if ((constraints == null) || (constraints instanceof String)) {
                addLayoutComponent((String) constraints, comp);
            } else {
                throw new IllegalArgumentException("cannot add to layout: constraint must be a string (or null)");
            }
        }
    }

    /**
     * @deprecated replaced by <code>addLayoutComponent(Component, Object)</code>.
     */
    public void addLayoutComponent(String name, Component comp) {
        synchronized (comp.getTreeLock()) {
            /* Special case:  treat null the same as "Chart". */
            if (name == null) {
                name = "Chart";
            }

            if ("Chart".equals(name) || "Table".equals(name)) {
                chart = comp;
            } else if ("Title".equals(name)) {
                title = comp;
            } else if ("XLabel".equals(name)) {
                xLabel = comp;
            } else if ("YLabel".equals(name)) {
                yLabel = comp;
            } else {
                throw new IllegalArgumentException("cannot add to layout: unknown constraint: " + name);
            }
        }
    }

    /**
     * Removes the specified component from the layout.
     *
     * @param comp the component to be removed
     */
    public void removeLayoutComponent(Component comp) {
        synchronized (comp.getTreeLock()) {
            if (comp == chart) {
                chart = null;
            } else if (comp == title) {
                title = null;
            } else if (comp == xLabel) {
                xLabel = null;
            } else if (comp == yLabel) {
                yLabel = null;
            }
        }
    }

    /**
     * Returns the minimum dimensions needed to layout the components
     * contained in the specified target container.
     *
     * @param target the Container on which to do the layout
     * @see Container
     * @see #preferredLayoutSize
     */
    public Dimension minimumLayoutSize(Container target) {
        synchronized (target.getTreeLock()) {
            Dimension dim = new Dimension(0, 0);

            if ((chart != null) && chart.isVisible()) {
                Dimension d = chart.getMinimumSize();
                dim.width = d.width;
                dim.height = d.height;
            }

            if ((xLabel != null) && xLabel.isVisible()) {
                Dimension d = xLabel.getMinimumSize();
                dim.width = Math.max(d.width, dim.width);
                dim.height += d.height + vgap;
            }

            if ((yLabel != null) && yLabel.isVisible()) {
                Dimension d = yLabel.getMinimumSize();
                dim.width += d.width + hgap;
                dim.height = Math.max(d.height, dim.height);
            }

            if ((title != null) && title.isVisible()) {
                Dimension d = title.getMinimumSize();
                dim.width = Math.max(d.width, dim.width);
                dim.height += d.height + vgap;
            }

            Insets insets = target.getInsets();
            dim.width += insets.left + insets.right;
            dim.height += insets.top + insets.bottom;

            return dim;
        }
    }

    /**
     * Returns the preferred dimensions for this layout given the components
     * in the specified target container.
     *
     * @param target the component which needs to be laid out
     * @see Container
     * @see #minimumLayoutSize
     */
    public Dimension preferredLayoutSize(Container target) {
        synchronized (target.getTreeLock()) {
            Dimension dim = new Dimension(0, 0);

            if ((chart != null) && chart.isVisible()) {
                Dimension d = chart.getPreferredSize();
                dim.width = d.width;
                dim.height = d.height;
            }

            if ((xLabel != null) && xLabel.isVisible()) {
                Dimension d = xLabel.getPreferredSize();
                dim.width = Math.max(d.width, dim.width);
                dim.height += d.height + vgap;
            }

            if ((yLabel != null) && yLabel.isVisible()) {
                Dimension d = yLabel.getPreferredSize();
                dim.width += d.width + hgap;
                dim.height = Math.max(d.height, dim.height);
            }

            if ((title != null) && title.isVisible()) {
                Dimension d = title.getPreferredSize();
                dim.width = Math.max(d.width, dim.width);
                dim.height += d.height + vgap;
            }

            Insets insets = target.getInsets();
            dim.width += insets.left + insets.right;
            dim.height += insets.top + insets.bottom;

            return dim;
        }
    }

    /**
     * Returns the maximum dimensions for this layout given the components
     * in the specified target container.
     *
     * @param target the component which needs to be laid out
     */
    public Dimension maximumLayoutSize(Container target) {
        return new Dimension(Integer.MAX_VALUE, Integer.MAX_VALUE);
    }

    /**
     * Returns the alignment along the x axis.  This specifies how
     * the component would like to be aligned relative to other
     * components.  The value should be a number between 0 and 1
     * where 0 represents alignment along the origin, 1 is aligned
     * the furthest away from the origin, 0.5 is centered, etc.
     */
    public float getLayoutAlignmentX(Container parent) {
        return 0.5f;
    }

    /**
     * Returns the alignment along the y axis.  This specifies how
     * the component would like to be aligned relative to other
     * components.  The value should be a number between 0 and 1
     * where 0 represents alignment along the origin, 1 is aligned
     * the furthest away from the origin, 0.5 is centered, etc.
     */
    public float getLayoutAlignmentY(Container parent) {
        return 0.5f;
    }

    /**
     * Invalidates the layout, indicating that if the layout manager
     * has cached information it should be discarded.
     */
    public void invalidateLayout(Container target) {
    }

    /**
     * Lays out the specified container. This method will actually reshape the
     * components in the specified target container in order to satisfy the
     * constraints of the CompassLayout object.
     *
     * @param target the component being laid out
     * @see Container
     */
    public void layoutContainer(Container target) {
        synchronized (target.getTreeLock()) {
            Insets insets = target.getInsets();
            Dimension d;

            d = target.getSize();

            int top = insets.top;
            int bottom = d.height - insets.bottom;
            int left = insets.left;
            int right = d.width - insets.right;

            // defines the centre box.
            int top1 = top;
            int bottom1 = bottom;
            int left1 = left;

            if ((xLabel != null) && xLabel.isVisible()) {
                d = xLabel.getPreferredSize();
                bottom1 -= d.height + vgap;
            }
            if ((yLabel != null) && yLabel.isVisible()) {
                d = yLabel.getPreferredSize();
                left1 += d.width + hgap;
            }
            if ((title != null) && title.isVisible()) {
                d = title.getPreferredSize();
                top1 += d.height + vgap;
            }

            if ((xLabel != null) && xLabel.isVisible())
                xLabel.setBounds(left1, bottom1, right - left1, bottom - bottom1);

            if ((yLabel != null) && yLabel.isVisible())
                yLabel.setBounds(left, top1, left1 - left, bottom1 - top1);

            if ((title != null) && title.isVisible())
                title.setBounds(left, top, right - left, top1 - top);

            if ((chart != null) && chart.isVisible())
                chart.setBounds(left1, top1, right - left1, bottom1 - top1);
        }
    }

    /**
     * Returns the String representation of this CompassLayout's values.
     */
    public String toString() {
        return getClass().getName() + "[hgap=" + hgap + ",vgap=" + vgap + "]";
    }
}
