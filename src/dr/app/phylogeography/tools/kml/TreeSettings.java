/*
 * TreeSettings.java
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

/**
 * @author Andrew Rambaut
 * @author Philippe Lemey
 * @version $Id$
 */
public class TreeSettings {
    public TreeSettings(String name, String description, TreeType treeType) {
        this.name = name;
        this.description = description;
        this.treeType = treeType;
    }

    public String getName() {
        return name;
    }

    public void setName(String name) {
        this.name = name;
    }

    public String getDescription() {
        return description;
    }

    public void setDescription(String description) {
        this.description = description;
    }

    public TreeType getTreeType() {
        return treeType;
    }

    public void setTreeType(TreeType treeType) {
        this.treeType = treeType;
    }

    public double getRootFraction() {
        return rootFraction;
    }

    public void setRootFraction(double rootFraction) {
        this.rootFraction = rootFraction;
    }

    public double getArcAltitude() {
        return arcAltitude;
    }

    public void setArcAltitude(double arcAltitude) {
        this.arcAltitude = arcAltitude;
    }

    public boolean isTaxaVisible() {
        return taxaVisible;
    }

    public void setTaxaVisible(boolean taxaVisible) {
        this.taxaVisible = taxaVisible;
    }

    public BranchStyle getBranchStyle() {
        return branchStyle;
    }

    public SurfaceDecoration getLocationContourStyle() {
        return locationContourStyle;
    }

    public SurfaceDecoration getGroundContourStyle() {
        return groundContourStyle;
    }

    public void setGroundContourStyle(SurfaceDecoration groundContourStyle) {
        this.groundContourStyle = groundContourStyle;
    }

    public BranchStyle getProjectionStyle() {
        return projectionStyle;
    }

    private String name;
    private String description;

    private TreeType treeType;

    private double rootFraction = 0.05; // additional fraction of the total plotHeight for the root branch

    private double arcAltitude = 10000; // this is the factor with which to multiply the time of the branch to get the altitude for that branch in the surface Tree

    private boolean taxaVisible = false;

    private BranchStyle branchStyle = new BranchStyle();

    private SurfaceDecoration locationContourStyle = new SurfaceDecoration();
    private SurfaceDecoration groundContourStyle = new SurfaceDecoration();

    private BranchStyle projectionStyle = new BranchStyle();
}
