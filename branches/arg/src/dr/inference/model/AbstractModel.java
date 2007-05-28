/*
 * AbstractModel.java
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

import dr.inference.parallel.MPISerializable;
import org.w3c.dom.Document;
import org.w3c.dom.Element;

import java.util.ArrayList;

/**
 * A model that brings together a number of model components
 *
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @version $Id: AbstractModel.java,v 1.13 2006/08/17 15:30:08 rambaut Exp $
 */
public abstract class AbstractModel implements Model, ModelListener,
		ParameterListener, StatisticList, MPISerializable {

	public AbstractModel(String name) {
		this.name = name;
	}

	/**
	 * Adds a sub-model to this model. If the model is already in the
	 * list then it does nothing.
	 */
	public void addModel(Model model) {

		if (!models.contains(model)) {

			models.add(model);
			model.addModelListener(this);
		}
	}

	public void removeModel(Model model) {
		models.remove(model);
		model.removeModelListener(this);
	}

	public int getModelCount() {
		return models.size();
	}

	public final Model getModel(int i) {
		return models.get(i);
	}

	public final void addParameter(Parameter parameter) {
		parameters.add(parameter);
		parameter.addParameterListener(this);

		// parameters are also statistics
		addStatistic(parameter);
	}

	public final void removeParameter(Parameter parameter) {
		parameters.remove(parameter);
		parameter.removeParameterListener(this);

		// parameters are also statistics
		removeStatistic(parameter);
	}

	/**
	 * @param parameter
	 * @return true of the given parameter is contained in this model
	 */
	public final boolean hasParameter(Parameter parameter) {
		return parameters.contains(parameter);
	}

	/**
	 * Adds a model listener.
	 */
	public void addModelListener(ModelListener listener) {
		listenerHelper.addModelListener(listener);
	}

	/**
	 * remove a model listener.
	 */
	public void removeModelListener(ModelListener listener) {
		listenerHelper.removeModelListener(listener);
	}

	/**
	 * Fires a model changed event.
	 */
	public void fireModelChanged() {
		listenerHelper.fireModelChanged(this, this, -1);
	}

	public void fireModelChanged(Object object) {
		listenerHelper.fireModelChanged(this, object, -1);
	}

	public void fireModelChanged(Object object, int index) {
		listenerHelper.fireModelChanged(this, object, index);
	}

	public final int getParameterCount() {
		return parameters.size();
	}

	public final Parameter getParameter(int i) {
		return parameters.get(i);
	}

	/**
	 * @return the parameter of the component that is called name
	 */
	public final Parameter getParameter(String name) {

		int i, n = getParameterCount();
		Parameter parameter;
		for (i = 0; i < n; i++) {
			parameter = getParameter(i);
			String paramName = parameter.getParameterName();
			if (paramName.equals(name)) {
				return parameter;
			}
		}

		return null;
		//throw new IllegalArgumentException("Parameter named " + name + " not found.");
	}

	// **************************************************************
	// ModelListener IMPLEMENTATION
	// **************************************************************

	public final void modelChangedEvent(Model model, Object object, int index) {

//		String message = "  model: " + getModelName() + "/" + getId() + "  component: " + model.getModelName();
//		if (object != null) {
//			message += " object: " + object;
//		}
//		if (index != -1) {
//			message += " index: " + index;
//		}
//		System.out.println(message);

		handleModelChangedEvent(model, object, index);
	}

	abstract protected void handleModelChangedEvent(Model model, Object object, int index);

	// **************************************************************
	// ParameterListener IMPLEMENTATION
	// **************************************************************

	public final void parameterChangedEvent(Parameter parameter, int index) {
		handleParameterChangedEvent(parameter, index);
		listenerHelper.fireModelChanged(this, parameter);
	}

	/**
	 * This method is called whenever a parameter is changed.
	 * Typically the model component sets a flag to recalculate intermediates
	 * Recalculation is typically done when the modelComponent is asked for some information
	 * that requires them. This mechanism is 'lazy' so that this method
	 * can be safely called a lot of times, without excessive calculation occurring.
	 */
	protected abstract void handleParameterChangedEvent(Parameter parameter, int index);

	// **************************************************************
	// Model IMPLEMENTATION
	// **************************************************************

	public final void storeModelState() {
		if (isValidState) {
			//System.out.println("STORE MODEL: " + getModelName() + "/" + getId());

			for (int i = 0; i < models.size(); i++) {
				getModel(i).storeModelState();
			}
			for (Parameter parameter : parameters) {
				parameter.storeParameterValues();
			}

			storeState();
			isValidState = false;
		}
	}

	public final void restoreModelState() {
		if (!isValidState) {
			//System.out.println("RESTORE MODEL: " + getModelName() + "/" + getId());

			for (int i = 0; i < parameters.size(); i++) {
				getParameter(i).restoreParameterValues();
			}
			for (int i = 0; i < models.size(); i++) {
				getModel(i).restoreModelState();
			}

			restoreState();
			isValidState = true;
		}
	}

	public final void acceptModelState() {
		if (!isValidState) {
			//System.out.println("ACCEPT MODEL: " + getModelName() + "/" + getId());

			for (int i = 0; i < parameters.size(); i++) {
				getParameter(i).acceptParameterValues();
			}
			for (int i = 0; i < models.size(); i++) {
				getModel(i).acceptModelState();
			}

			acceptState();
			isValidState = true;
		}
	}

	public boolean isValidState() {
		return isValidState;
	}

	public final String getModelName() {
		return name;
	}


	/**
	 * Additional state information, outside of the sub-model is stored by this call.
	 */
	protected abstract void storeState();

	/**
	 * After this call the model is guaranteed to have returned its extra state information to
	 * the values coinciding with the last storeState call.
	 * Sub-models are handled automatically and do not need to be considered in this method.
	 */
	protected abstract void restoreState();

	/**
	 * This call specifies that the current state is accept. Most models will not need to do anything.
	 * Sub-models are handled automatically and do not need to be considered in this method.
	 */
	protected abstract void acceptState();

	// **************************************************************
	// StatisticList IMPLEMENTATION
	// **************************************************************

	public final void addStatistic(Statistic statistic) {

		statistics.add(statistic);
	}

	public final void removeStatistic(Statistic statistic) {

		statistics.remove(statistic);
	}

	/**
	 * @return the number of statistics of this component.
	 */
	public int getStatisticCount() {

		return statistics.size();
	}

	/**
	 * @return the ith statistic of the component
	 */
	public Statistic getStatistic(int i) {

		return statistics.get(i);
	}

	public final Statistic getStatistic(String name) {

		int i, n = getStatisticCount();
		Statistic statistic;
		for (i = 0; i < n; i++) {
			statistic = getStatistic(i);
			if (name.equals(statistic.getStatisticName())) {
				return statistic;
			}
		}

		return null;
	}

	// **************************************************************
	// Identifiable IMPLEMENTATION
	// **************************************************************

	private String id = null;

	public void setId(String id) {
		this.id = id;
	}

	public String getId() {
		return id;
	}


	public String toString() {
		if (id != null) {
			return id;
		} else if (name != null) {
			return name;
		}
		return super.toString();
	}

	// **************************************************************
	// XMLElement IMPLEMENTATION
	// **************************************************************

	public Element createElement(Document d) {
		throw new RuntimeException("Not implemented!");
	}

	// **************************************************************
	// MPI IMPLEMENTATION
	// **************************************************************


	public void sendState(int toRank) {

		// Iterate through child models
		for (Model model : models) {
			((AbstractModel) model).sendState(toRank);
		}
		// Send current model parameters
		for (Parameter parameter : parameters) {
			((Parameter.Abstract) parameter).sendState(toRank);
		}

	}


	public void sendStateNoParameters(int toRank) {

		// Iterate through child models
		for (Model model : models) {
			((AbstractModel) model).sendState(toRank);
		}

	}

	public void receiveStateNoParameters(int fromRank) {
		for (Model model : models) {
			((AbstractModel) model).receiveState(fromRank);
		}
	}

	public void receiveState(int fromRank) {
		for (Model model : models) {
			((AbstractModel) model).receiveState(fromRank);
		}
		// Send current model parameters
		for (Parameter parameter : parameters) {
			((Parameter.Abstract) parameter).receiveState(fromRank);
		}


	}

	boolean isValidState = true;

	protected Model.ListenerHelper listenerHelper = new Model.ListenerHelper();

	private ArrayList<Model> models = new ArrayList<Model>();
	private ArrayList<Parameter> parameters = new ArrayList<Parameter>();
	private ArrayList<Statistic> statistics = new ArrayList<Statistic>();

	private String name;
}