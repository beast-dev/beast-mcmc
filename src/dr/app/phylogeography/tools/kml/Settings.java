/*
 * Settings.java
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

import java.util.List;
import java.util.ArrayList;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id$
 */
public class Settings {
    public Settings(AnalysisType analysisType) {
        this.analysisType = analysisType;
    }

    public AnalysisType getAnalysisType() {
        return analysisType;
    }

    public void setAnalysisType(AnalysisType analysisType) {
        this.analysisType = analysisType;
    }

    public String getTraitName() {
        return traitName;
    }

    public void setTraitName(String traitName) {
        this.traitName = traitName;
    }

    public String getLatitudeName() {
        return latitudeName;
    }

    public void setLatitudeName(String latitudeName) {
        this.latitudeName = latitudeName;
    }

    public String getLongitudeName() {
        return longitudeName;
    }

    public void setLongitudeName(String longitudeName) {
        this.longitudeName = longitudeName;
    }

    public double getMostRecentDate() {
        return mostRecentDate;
    }

    public void setMostRecentDate(double mostRecentDate) {
        this.mostRecentDate = mostRecentDate;
    }

    public double getAgeCutOff() {
        return ageCutOff;
    }

    public void setAgeCutOff(double ageCutOff) {
        this.ageCutOff = ageCutOff;
    }

    public double getPlotAltitude() {
        return plotAltitude;
    }

    public void setPlotAltitude(double plotAltitude) {
        this.plotAltitude = plotAltitude;
    }

    public double getColumnRadius() {
        return columnRadius;
    }

    public void setColumnRadius(double columnRadius) {
        this.columnRadius = columnRadius;
    }

    public int getTimeDivisionCount() {
        return timeDivisionCount;
    }

    public void setTimeDivisionCount(int timeDivisionCount) {
        this.timeDivisionCount = timeDivisionCount;
    }

    public TreeSettings getAltitudeTreeSettings() {
        return altitudeTreeSettings;
    }

    public TreeSettings getGroundTreeSettings() {
        return groundTreeSettings;
    }

    public SurfaceDecoration getGroundContours() {
        return groundContours;
    }

    public SurfaceDecoration getProjections() {
        return projections;
    }

    public SurfaceDecoration getTaxonLabels() {
        return taxonLabels;
    }

    public SurfaceDecoration getLocationLabels() {
        return locationLabels;
    }

    private AnalysisType analysisType;

    private String traitName = "location";
    private String latitudeName = "location1";
    private String longitudeName = "location2";

    private double ageCutOff = 0.0; // upper bound for when time starts - 0 to include entire tree
    private double mostRecentDate;

    private double plotAltitude = 0;
    private double columnRadius = 1;

    private final TreeSettings altitudeTreeSettings = new TreeSettings("altitudeTree", "", TreeType.RECTANGLE_TREE);
    private final TreeSettings groundTreeSettings = new TreeSettings("groundTree", "", TreeType.SURFACE_TREE);
    private final SurfaceDecoration groundContours = new SurfaceDecoration();
    private final SurfaceDecoration projections = new SurfaceDecoration();
    private final SurfaceDecoration taxonLabels = new SurfaceDecoration();
    private final SurfaceDecoration locationLabels = new SurfaceDecoration();

    private int timeDivisionCount = 0;
}
