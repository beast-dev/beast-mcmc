/*
 * CompoundModel.java
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

package dr.inference.model;

import java.util.ArrayList;
import java.util.List;

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
            // add all listeners to this model
            for( ModelListener ml : listeners ) {
                 model.addModelListener(ml);
            }
        }
	}

	public int getModelCount() {
		return models.size();
	}

	public final Model getModel(int i) {
		return models.get(i);
	}

    public boolean isUsed() {
        return listeners.size() > 0;
    }

    public void addModelListener(ModelListener listener) {
        // add listener to all models comprizing this compund model - a change in any one of them
        // means the compund model changed

        listeners.add(listener);
        for( Model m : models ) {
            m.addModelListener(listener);
        }

        //throw new IllegalArgumentException("Compound models don't have listeners");
	}

	public void removeModelListener(ModelListener listener) {
        for( Model m : models ) {
            m.removeModelListener(listener);
        }
        listeners.remove(listener);
       // throw new IllegalArgumentException("Compound models don't have listeners");
	}

	public void storeModelState() {
        for (Model model : models) {
            model.storeModelState();
        }
    }

	public void restoreModelState() {
        for (Model model : models) {
            model.restoreModelState();
        }
    }

	public void acceptModelState() {
        for (Model model : models) {
            model.acceptModelState();
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

	public int getVariableCount() { return 0; }

	public Variable getVariable(int index) {
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

	private String name = null;
	private final ArrayList<Model> models = new ArrayList<Model>();
    private final List<ModelListener> listeners = new ArrayList<ModelListener>();
}

