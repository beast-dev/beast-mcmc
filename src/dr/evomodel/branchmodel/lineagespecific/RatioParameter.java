/*
 * RatioParameter.java
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

import java.util.List;

import dr.inference.model.Bounds;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.inference.model.VariableListener;

public class RatioParameter extends Parameter.Abstract implements
		VariableListener {

	private Parameter parameter1;
	private Parameter parameter2;
	private Bounds bounds = null;

	public RatioParameter(Parameter parameter1, Parameter parameter2) {

		this.parameter1 = parameter1;
		this.parameter2 = parameter2;

		this.parameter1.addVariableListener(this);
		this.parameter2.addVariableListener(this);

	}// END: Constructor

	@Override
	public double getParameterValue(int dim) {

		double value = parameter1.getParameterValue(dim)
				/ parameter2.getParameterValue(dim);

		return value;
	}// END: getParameterValue

	@Override
	public void setParameterValue(int dim, double value) {
		throw new RuntimeException("Not implemented");
	}

	@Override
	public void setParameterValueQuietly(int dim, double value) {
		throw new RuntimeException("Not implemented");
	}

	@Override
	public void setParameterValueNotifyChangedAll(int dim, double value) {
		throw new RuntimeException("Not implemented");
	}

	@Override
	public String getParameterName() {
		if (getId() == null) {

			StringBuilder sb = new StringBuilder("ratio");
			sb.append(parameter1.getId()).append(".")
					.append(parameter2.getId());
			setId(sb.toString());
		}

		return getId();
	}// END: getParameterName

	@Override
	public void addBounds(Bounds<Double> bounds) {
		this.bounds = bounds;
	}

	@Override
	public Bounds<Double> getBounds() {
		return bounds;
	}

	@Override
	public void addDimension(int index, double value) {
		throw new RuntimeException("Not yet implemented.");
	}

	@Override
	public double removeDimension(int index) {
		throw new RuntimeException("Not yet implemented.");
	}

	@Override
	public void variableChangedEvent(Variable variable, int index,
			dr.inference.model.Variable.ChangeType type) {
		fireParameterChangedEvent(index, type);
	}

	@Override
	protected void storeValues() {
		parameter1.storeParameterValues();
		parameter2.storeParameterValues();
	}

	@Override
	protected void restoreValues() {
		parameter1.restoreParameterValues();
		parameter2.restoreParameterValues();
	}

	@Override
	protected void acceptValues() {
		parameter1.acceptParameterValues();
		parameter2.acceptParameterValues();
	}

	@Override
	protected void adoptValues(Parameter source) {
		throw new RuntimeException("Not implemented");
	}

}// END: class
