/*
 * CountableRealizationsParameter.java
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

package dr.evomodel.branchmodel.lineagespecific;

import java.util.LinkedList;

import dr.inference.model.Bounds;
import dr.inference.model.CompoundParameter;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.inference.model.VariableListener;

/**
 * @author Marc Suchard
 * @author Filip Bielejec
 * 
 */
@SuppressWarnings("serial")
public class CountableRealizationsParameter extends Parameter.Abstract implements VariableListener {

	private final Parameter categoriesParameter;
	private final CompoundParameter uniquelyRealizedParameters;
	
	
	private final int dim;
	private final int realizationCount;

	private final int uniqueRealizationCount;
	
    private final LinkedList<Parameter> paramList;
    private Bounds<Double> bounds = null;
	
    public CountableRealizationsParameter(Parameter categoriesParameter, //
    		 CompoundParameter uniquelyRealizedParameters //
    		) {
    	
    	this.categoriesParameter = categoriesParameter;
    	this.uniquelyRealizedParameters = uniquelyRealizedParameters;
    	
    	dim = uniquelyRealizedParameters.getParameter(0).getDimension();
    	// TODO Make sure all parameters have same dimension
    	
    	realizationCount = categoriesParameter.getDimension();
    	
    	paramList = new LinkedList<Parameter>();
    	paramList.add(categoriesParameter);
    	paramList.add(uniquelyRealizedParameters);
    	
//        this.paramList = parameter;
//        for (Parameter p : paramList) {
//            p.addVariableListener(this);
//        }
    	
    	uniqueRealizationCount = uniquelyRealizedParameters.getDimension();
    
    }//END: Constructor

    public int getDimension() {
    	return dim * realizationCount; //  paramList.get(0).getDimension(); // Unwritten contract
    }

    protected void storeValues() {
        for (Parameter p : paramList) {
            p.storeParameterValues();
        }
    }

    protected void restoreValues() {
        for (Parameter p : paramList) {
            p.restoreParameterValues();
        }
    }

    protected void acceptValues() {
        for (Parameter p : paramList) {
            p.acceptParameterValues();
        }
    }

    protected void adoptValues(Parameter source) {
        throw new RuntimeException("Not implemented");
    }

    public double getParameterValue(int index) {
    	
//    	int whichCategoryIndex = index % dim; // TODO Maybe?
//    	int whichDimIndex = index - whichCategoryIndex * dim; // TODO Maybe?
//    	Parameter param = uniquelyRealizedParameters.getParameter((int) categoriesParameter.getParameterValue(whichCategoryIndex));
//    	return param.getParameterValue(whichDimIndex);     
    	
    	int whichCategoryIndex = categoriesParameter.getValue(index).intValue();//
    	int whichDimIndex = 0;
    	
    	Parameter param = uniquelyRealizedParameters.getParameter(whichCategoryIndex);
    	
    	return param.getParameterValue(whichDimIndex);     
    }//END: getParameterValue

	public void setParameterValue(int dim, double value) {

		int whichCategoryIndex = (int) categoriesParameter.getParameterValue(dim);
		uniquelyRealizedParameters.setParameterValue(whichCategoryIndex, value);

	}//END: setParameterValue

    public void setParameterValueQuietly(int dim, double value) {
        throw new RuntimeException("Not implemented");
    }

    public void setParameterValueNotifyChangedAll(int dim, double value){
        throw new RuntimeException("Not implemented");
    }

    public String getParameterName() {
        if (getId() == null) {
            StringBuilder sb = new StringBuilder("product");
            for (Parameter p : paramList) {
                sb.append(".").append(p.getId());
            }
            setId(sb.toString());
        }
        return getId();
    }

    public void addBounds(Bounds<Double> bounds) {
        this.bounds = bounds;
    }

    public Bounds<Double> getBounds() {
        if (bounds == null) {
            return paramList.get(0).getBounds(); // TODO
        } else {
            return bounds;
        }
    }

    public void addDimension(int index, double value) {
        throw new RuntimeException("Not yet implemented.");
    }

    public double removeDimension(int index) {
        throw new RuntimeException("Not yet implemented.");
    }

    public void variableChangedEvent(Variable variable, int index, ChangeType type) {
        fireParameterChangedEvent(index,type);
    }

}//END: class
