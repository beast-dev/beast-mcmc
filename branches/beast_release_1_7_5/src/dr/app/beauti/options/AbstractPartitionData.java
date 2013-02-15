/*
 * AbstractPartitionData.java
 *
 * Copyright (c) 2002-2011 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.app.beauti.options;

import dr.evolution.alignment.Patterns;
import dr.evolution.datatype.DataType;
import dr.evolution.distance.DistanceMatrix;
import dr.evolution.distance.JukesCantorDistanceMatrix;
import dr.evolution.util.TaxonList;

import java.io.Serializable;
import java.util.List;

/**
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @author Walter Xie
 */
public abstract class AbstractPartitionData implements Serializable {


    protected String fileName;
    protected String name;
    protected List<TraitData> traits;

    protected BeautiOptions options;
    protected PartitionSubstitutionModel model;
    protected PartitionClockModel clockModel;
    protected PartitionTreeModel treeModel;

    protected double meanDistance;

    protected DistanceMatrix distances;

    public AbstractPartitionData(BeautiOptions options, String name, String fileName) {
        this.options = options;
        this.name = name;
        this.fileName = fileName;
    }

    protected void calculateMeanDistance(Patterns patterns) {
        if (patterns != null) {
            distances = new JukesCantorDistanceMatrix(patterns);
            meanDistance = distances.getMeanDistance();
        } else {
            distances = null;
            meanDistance = 0.0;
        }
    }

    public String getFileName() {
        return fileName;
    }

    public String getName() {
        return name;
    }

    public void setName(String name) {
        this.name = name;
    }

    public String toString() {
        return getName();
    }

    public List<TraitData> getTraits() {
        return traits;
    }

    public double getMeanDistance() {
        return meanDistance;
    }

    public DistanceMatrix getDistances() {
        return distances;
    }

    public void setPartitionSubstitutionModel(PartitionSubstitutionModel model) {
        this.model = model;
    }

    public PartitionSubstitutionModel getPartitionSubstitutionModel() {
        return this.model;
    }

    public void setPartitionClockModel(PartitionClockModel clockModel) {
        this.clockModel = clockModel;
    }

    public PartitionClockModel getPartitionClockModel() {
        return clockModel;
    }

    public PartitionTreeModel getPartitionTreeModel() {
        return treeModel;
    }

    public void setPartitionTreeModel(PartitionTreeModel treeModel) {
        this.treeModel = treeModel;
    }

    public int getTaxonCount() {
        if (getTaxonList() != null) {
            return getTaxonList().getTaxonCount();
        } else {
            // is a trait
            return -1;
        }
    }

    public boolean isCreatedFromTrait() {
       return traits != null;
    }

    public abstract String getPrefix(); // be careful of microsatellite PartitionPattern

    public abstract TaxonList getTaxonList();

    public abstract int getSiteCount();

    public abstract DataType getDataType();

    public abstract String getDataDescription();

}
