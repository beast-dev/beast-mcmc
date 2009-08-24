/*
 * Copyright (C) 2002-2009 Alexei Drummond and Andrew Rambaut
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
 * BEAST is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package dr.app.beauti.options;

import java.util.ArrayList;
import java.util.List;

import dr.evolution.alignment.Alignment;
import dr.evolution.alignment.Patterns;
import dr.evolution.alignment.SiteList;
import dr.evolution.distance.DistanceMatrix;
import dr.evolution.distance.JukesCantorDistanceMatrix;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @author Walter Xie
 */
public class PartitionData extends PartitionOptions {

    private final String fileName;
    private final Alignment alignment;
    private final double meanDistance;

    private String name;
    private boolean coding;

    private int fromSite;
    private int toSite;
    private int every = 1;

    private PartitionSubstitutionModel model;
    private PartitionClockModel clockModel;
    private PartitionTreeModel treeModel;

    public PartitionData(String name, String fileName, Alignment alignment) {
        this(name, fileName, alignment, -1, -1, 1);
    }

    public PartitionData(String name, String fileName, Alignment alignment, int fromSite, int toSite, int every) {
        this.name = name;
        this.fileName = fileName;
        this.alignment = alignment;
        this.coding = false;

        this.fromSite = fromSite;
        this.toSite = toSite;
        this.every = every;

        List<PartitionData> p = new ArrayList<PartitionData>();
        p.add(this);
        meanDistance = calculateMeanDistance(p);
        	
//        Patterns patterns = new Patterns(alignment);
//        DistanceMatrix distances = new JukesCantorDistanceMatrix(patterns);
//        meanDistance = distances.getMeanDistance();
//        meanDistance = 0.0;
    }
    
    public double getMeanDistance() {
        return meanDistance;
    }

    public String getFileName() {
        return fileName;
    }

    public Alignment getAlignment() {
        return alignment;
    }

    public String getName() {
        return name;
    }

    public void setName(String name) {
        this.name = name;
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

    public boolean isCoding() {
        return coding;
    }

    public void setCoding(boolean coding) {
        this.coding = coding;
    }

    public int getFromSite() {
        return fromSite;
    }

    public int getToSite() {
        return toSite;
    }

    public int getEvery() {
        return every;
    }

    public int getSiteCount() {
        int from = getFromSite();
        if (from < 1) {
            from = 1;
        }
        int to = getToSite();
        if (to < 1) {
            to = alignment.getSiteCount();
        }
        return (to - from + 1) / every;
    }

    public int getTaxaCount() {
        int n = alignment.getSequenceCount();

        if (n > 0) {
            return n;
        } else {
            return 0;
        }
    }

    public String toString() {
        return getName();
    }

	@Override
	public void selectOperators(List<Operator> ops) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void selectParameters(List<Parameter> params) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public String getPrefix() {
		// TODO Auto-generated method stub
		return null;
	}

}
