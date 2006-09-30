/*
 * CompassLayout.java
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

package org.virion.jam.layouts;

import java.awt.*;

/**
 * A layout manager similar to BorderLayout but with 8 compass directions. It will
 * layout a container using members named "North", "NorthEast", "East", etc. and
 * "Center".
 * <p/>
 * The "North", "South", "East" and "West" components get layed out
 * according to their preferred sizes and the constraints of the
 * container's size. The "Center" component will get any space left
 * over. The corner components are layed out accordingly
 *
 * @author Andrew Rambaut
 */

public class CompassLayout implements LayoutManager2 {
    int hgap;
    int vgap;

    Component north, northWest;
    Component west, southWest;
    Component east, northEast;
    Component south, southEast;
    Component center;

    /**
     * Constructs a new CompassLayout.
     */
    public CompassLayout() {
        this(0, 0);
    }

    /**
     * Constructs a CompassLayout with the specified gaps.
     *
     * @param hgap the horizontal gap
     * @param vgap the vertical gap
     */
    public CompassLayout(int hgap, int vgap) {
        this.hgap = hgap;
        this.vgap = vgap;
    }

    /**
     * Returns the horizontal gap between components.
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
     * Returns the vertical gap between components.
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
     * @param name the String name
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
            /* Special case:  treat null the same as "Center". */
            if (name == null) {
                name = "Center";
            }

            if ("Center".equals(name)) {
                center = comp;
            } else if ("North".equals(name)) {
                north = comp;
            } else if ("South".equals(name)) {
                south = comp;
            } else if ("East".equals(name)) {
                east = comp;
            } else if ("West".equals(name)) {
                west = comp;
            } else if ("NorthEast".equals(name)) {
                northEast = comp;
            } else if ("SouthEast".equals(name)) {
                southEast = comp;
            } else if ("NorthWest".equals(name)) {
                northWest = comp;
            } else if ("SouthWest".equals(name)) {
                southWest = comp;
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
            if (comp == center) {
                center = null;
            } else if (comp == north) {
                north = null;
            } else if (comp == south) {
                south = null;
            } else if (comp == east) {
                east = null;
            } else if (comp == west) {
                west = null;
            } else if (comp == northEast) {
                northEast = null;
            } else if (comp == southEast) {
                southEast = null;
            } else if (comp == northWest) {
                northWest = null;
            } else if (comp == southWest) {
                southWest = null;
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

            if ((east != null) && east.isVisible()) {
                Dimension d = east.getMinimumSize();
                dim.width += d.width + hgap;
                dim.height = Math.max(d.height, dim.height);
            }
            if ((west != null) && west.isVisible()) {
                Dimension d = west.getMinimumSize();
                dim.width += d.width + hgap;
                dim.height = Math.max(d.height, dim.height);
            }
            if ((center != null) && center.isVisible()) {
                Dimension d = center.getMinimumSize();
                dim.width += d.width;
                dim.height = Math.max(d.height, dim.height);
            }
            if ((north != null) && north.isVisible()) {
                Dimension d = north.getMinimumSize();
                dim.width = Math.max(d.width, dim.width);
                dim.height += d.height + vgap;
            }
            if ((south != null) && south.isVisible()) {
                Dimension d = south.getMinimumSize();
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

            if ((east != null) && east.isVisible()) {
                Dimension d = east.getPreferredSize();
                dim.width += d.width + hgap;
                dim.height = Math.max(d.height, dim.height);
            }
            if ((west != null) && west.isVisible()) {
                Dimension d = west.getPreferredSize();
                dim.width += d.width + hgap;
                dim.height = Math.max(d.height, dim.height);
            }
            if ((center != null) && center.isVisible()) {
                Dimension d = center.getPreferredSize();
                dim.width += d.width;
                dim.height = Math.max(d.height, dim.height);
            }
            if ((north != null) && north.isVisible()) {
                Dimension d = north.getPreferredSize();
                dim.width = Math.max(d.width, dim.width);
                dim.height += d.height + vgap;
            }
            if ((south != null) && south.isVisible()) {
                Dimension d = south.getPreferredSize();
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
            Dimension d, dn, ds, de, dw;

            d = target.getSize();

            int top = insets.top;
            int bottom = d.height - insets.bottom;
            int left = insets.left;
            int right = d.width - insets.right;

            // defines the centre box.
            int top1 = top;
            int bottom1 = bottom;
            int left1 = left;
            int right1 = right;

            if ((north != null) && north.isVisible()) {
                d = north.getPreferredSize();
                top1 += d.height + vgap;
            }
            if ((south != null) && south.isVisible()) {
                d = south.getPreferredSize();
                bottom1 -= d.height + vgap;
            }
            if ((east != null) && east.isVisible()) {
                d = east.getPreferredSize();
                right1 -= d.width + hgap;
            }
            if ((west != null) && west.isVisible()) {
                d = west.getPreferredSize();
                left1 += d.width + hgap;
            }

            if ((north != null) && north.isVisible())
                north.setBounds(left1, top, right1 - left1, top1 - top);

            if ((south != null) && south.isVisible())
                south.setBounds(left1, bottom1, right1 - left1, bottom - bottom1);

            if ((east != null) && east.isVisible())
                east.setBounds(right1, top1, right - right1, bottom1 - top1);

            if ((west != null) && west.isVisible())
                west.setBounds(left, top, left1 - left, bottom - top);

            if ((center != null) && center.isVisible())
                center.setBounds(left1, top1, right1 - left1, bottom1 - top1);

            if ((northWest != null) && northWest.isVisible())
                northWest.setBounds(left, top, left1 - left, top1 - top);

            if ((southWest != null) && southWest.isVisible())
                southWest.setBounds(left, bottom1, left1 - left, bottom - bottom1);

            if ((northEast != null) && northEast.isVisible())
                northEast.setBounds(right1, top, right - right1, top1 - top);

            if ((southEast != null) && southEast.isVisible())
                southEast.setBounds(right1, bottom1, right - right1, bottom - bottom1);
        }
    }

    /**
     * Returns the String representation of this CompassLayout's values.
     */
    public String toString() {
        return getClass().getName() + "[hgap=" + hgap + ",vgap=" + vgap + "]";
    }
}
