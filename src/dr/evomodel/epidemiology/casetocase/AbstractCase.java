/*
 * AbstractCase.java
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

package dr.evomodel.epidemiology.casetocase;

import dr.evolution.util.Date;
import dr.evolution.util.Taxa;
import dr.inference.model.AbstractModel;
import dr.inference.model.Parameter;

import java.util.ArrayList;

/**
 * Abstract class for outbreak; best implemented as an inner class in implementations of AbstractOutbreak
 *
 * User: Matthew Hall
 * Date: 15/04/13
 */


public abstract class AbstractCase extends AbstractModel {

    public AbstractCase(String name){
        super(name);
    }

    protected String caseID;
    protected Taxa associatedTaxa;
    protected double endOfInfectiousTime;
    protected boolean wasEverInfected;
    protected Parameter infectionBranchPosition;

    public String getName(){
        return caseID;
    }

    public Taxa getAssociatedTaxa() {
        return associatedTaxa;
    }

    public double getEndTime(){
        return endOfInfectiousTime;
    }

    public abstract boolean noninfectiousYet(double time);

    public String toString(){
        return caseID;
    }

    public abstract double[] getCoords();

    public boolean wasEverInfected(){
        return wasEverInfected;
    }

    public Parameter getInfectionBranchPosition(){
        return infectionBranchPosition;
    }

    public abstract ArrayList<Date> getExaminationTimes();

    public void setInfectionBranchPosition(double value){
        infectionBranchPosition.setParameterValue(0, value);
    }

    public void setEverInfected(boolean value){
        wasEverInfected = value;
    }


}

