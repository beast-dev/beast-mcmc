/*
 * ExposeRateCategoriesWrapper.java
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

package dr.evomodel.substmodel;

import dr.inference.model.Bounds;
import dr.inference.model.Model;
import dr.inference.model.ModelListener;
import dr.inference.model.Parameter;

/**
 * @author Marc A. Suchard
 * @author Rory E. Wasiolek
 */
public class ExposeRateCategoriesWrapper extends Parameter.Abstract implements ModelListener {

    public ExposeRateCategoriesWrapper(MarkovModulatedSubstitutionModel mmSubstModel) {
        super();
        this.mmSubstModel = mmSubstModel;
    }

    public int getDimension() {
        return mmSubstModel.getNumBaseModel();
    }

    @Override
    protected void storeValues() {
        // Do nothing
    }

    @Override
    protected void restoreValues() {
        // Do nothing
    }

    @Override
    protected void acceptValues() {
        // Do nothing
    }

    @Override
    protected void adoptValues(Parameter source) {
        // Do nothing
    }

    /**
     * @param dim the index of the parameter dimension of interest
     * @return the parameter's scalar value in the given dimension
     */
    @Override
    public double getParameterValue(int dim) {
        return mmSubstModel.getModelRateScalar(dim);
    }

    /**
     * sets the scalar value in the given dimension of this parameter
     *
     * @param dim   the index of the dimension to set
     * @param value the value to set
     */
    @Override
    public void setParameterValue(int dim, double value) {
        throw new RuntimeException("Not implemented for wrapper.");
    }

    /**
     * sets the scalar value in the given dimensin of this parameter to val, without firing any events
     *
     * @param dim   the index of the dimension to set
     * @param value the value to set
     */
    @Override
    public void setParameterValueQuietly(int dim, double value) {
        throw new RuntimeException("Not implemented for wrapper.");
    }

    /**
     * sets the scalar value in the given dimensin of this parameter to val,
     * and notifies that values in all dimension have been changed
     *
     * @param dim   the index of the dimension to set
     * @param value the value to set
     */
    @Override
    public void setParameterValueNotifyChangedAll(int dim, double value) {
        throw new RuntimeException("Not implemented for wrapper.");
    }

    /**
     * @return the name of this parameter
     */
    @Override
    public String getParameterName() {
        return "useful name";
    }

    /**
     * Adds new bounds to this parameter
     *
     * @param bounds to add
     */
    @Override
    public void addBounds(Bounds<Double> bounds) {
        throw new RuntimeException("Not implemented for wrapper.");
    }

    /**
     * @return the intersection of all bounds added to this parameter
     */
    @Override
    public Bounds<Double> getBounds() {
        throw new RuntimeException("Not implemented for wrapper.");
    }

    /**
     * Adds an extra dimension at the given index
     *
     * @param index Index of the dimension to add
     * @param value value to save at end of new array
     */
    @Override
    public void addDimension(int index, double value) {
        throw new RuntimeException("Not implemented for wrapper.");
    }

    /**
     * Removes the specified dimension from parameter
     *
     * @param index Index of dimension to lose
     * @return the value of the dimension removed
     */
    @Override
    public double removeDimension(int index) {
        throw new RuntimeException("Not implemented for wrapper.");
    }

    private final MarkovModulatedSubstitutionModel mmSubstModel;


    /**
     * The model has changed. The model firing the event can optionally
     * supply a reference to an object and an index if appropriate. Use
     * of this extra information will be contingent on recognising what
     * model it was that fired the event.
     */
    @Override
    public void modelChangedEvent(Model model, Object object, int index) {
        if (model == mmSubstModel) { // TODO limit to passing along only exposed value changes
            fireParameterChangedEvent();
        }
    }

    /**
     * The model has been restored.
     * Required only for notification of non-models (say pure likelihoods) which depend on
     * models.
     *
     * @param model
     */
    @Override
    public void modelRestored(Model model) {
        // Do nothing
    }
}
