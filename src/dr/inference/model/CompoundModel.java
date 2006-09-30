/*
 * CompoundModel.java
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

package dr.inference.model;

import java.util.ArrayList;

/**
 * An interface that describes a model of some data.
 *
 * @version $Id: CompoundModel.java,v 1.8 2005/05/24 20:25:59 rambaut Exp $
 *
 * @author Alexei Drummond
 * @author Andrew Rambaut
 */

public class CompoundModel implements Model {

	public static final String COMPOUND_MODEL = "compoundModel";

	public CompoundModel(String name) { this.name = name; }
	
	
	public void addModel(Model model) {
		
		if ( !models.contains(model) ) {
		
			models.add(model);
		}
	}
	
	public int getModelCount() {
		return models.size();
	}

	public final Model getModel(int i) { 
		return (Model)models.get(i); 
	}

	public void addModelListener(ModelListener listener) {
		throw new IllegalArgumentException("Compound models don't have listeners");
	}
	
	public void removeModelListener(ModelListener listener) {
		throw new IllegalArgumentException("Compound models don't have listeners");
	}

	public void storeModelState() {
		for (int i = 0; i < models.size(); i++) {
			((Model)models.get(i)).storeModelState();
		}
	}

	public void restoreModelState() {
		for (int i = 0; i < models.size(); i++) {
			((Model)models.get(i)).restoreModelState();
		}
	}

	public void acceptModelState() {
		for (int i = 0; i < models.size(); i++) {
			((Model)models.get(i)).acceptModelState();
		}
	}

	public boolean isValidState() { 
	
		for (int i = 0; i < models.size(); i++) {
			if (!getModel(i).isValidState()) {
				return false;
			}
		}
		
		return true;
	}
	
	public int getParameterCount() { return 0; }
		
	public Parameter getParameter(int index) {
		throw new IllegalArgumentException("Compound models don't have parameters");
	}
	
	public Parameter getParameter(String name) {
		throw new IllegalArgumentException("Compound models don't have parameters");
	}
	
	// **************************************************************
    // Identifiable IMPLEMENTATION
    // **************************************************************
	
	private String id = null;
		
	public void setId(String id) { this.id = id; }
	
	public String getId() { return id; }
	
	/**
	 * @return the name of this model
	 */
	public String getModelName() { return name; }
	
	/* AER - do we need a parser? 
	public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {
		
		public String getParserName() { return COMPOUND_MODEL; }
				
		public Object parseXMLObject(XMLObject xo) throws XMLParseException {
				
			CompoundModel compoundModel = new CompoundModel("model");
				
			int childCount = xo.getChildCount();
			for (int i = 0; i < childCount; i++) {
				Object xoc = xo.getChild(i);
				if (xoc instanceof Model) {
					compoundModel.addModel((Model)xoc);	
				} 
			}
			return compoundModel;
		}
		
		//************************************************************************
		// AbstractXMLObjectParser implementation
		//************************************************************************
		
		public String getParserDescription() {
			return "This element represents a combination of models.";
		}
		
		public Class getReturnType() { return CompoundModel.class; }
		
		public XMLSyntaxRule[] getSyntaxRules() { return rules; }
		
		private XMLSyntaxRule[] rules = new XMLSyntaxRule[] {
			new OneOrMoreRule(Model.class)
		};
	};*/
	
	String name = null;
	ArrayList models = new ArrayList();
}
		
