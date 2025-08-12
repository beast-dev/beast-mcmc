/*
 * PolymorphismAwareSubstitutionModel.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
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
 *
 */

package dr.evomodel.substmodel;

import dr.evolution.datatype.DataType;
import dr.inference.model.*;
import dr.util.Citable;
import dr.util.Citation;

import java.util.List;

/**
 * PoMo model of sequence evolution
 *
 * @author Xiang Ji
 * @author Nicola De Maio
 * @author Ben Redelings
 * @author Marc A. Suchard
 */
public class PolymorphismAwareSubstitutionModel extends AbstractModel implements ActionEnabledSubstitution, Citable {

    private final SubstitutionModel baseSubstitutionModel;
    private final int virtualPopSize;

    public PolymorphismAwareSubstitutionModel(SubstitutionModel baseSubstitutionModel, int virtualPopSize) {
        super("PoMoSubstitutionModel");
        this.baseSubstitutionModel = baseSubstitutionModel;
        this.virtualPopSize = virtualPopSize;
        addModel(baseSubstitutionModel);
    }


    @Override
    public int getNonZeroEntryCount() {
        return 0;
    }

    @Override
    public void getNonZeroEntries(int[] rowIndices, int[] colIndices, double[] values) {

    }

    @Override
    public void getTransitionProbabilities(double distance, double[] matrix) {
        throw new RuntimeException("Not yet implemented!");
    }

    @Override
    public EigenDecomposition getEigenDecomposition() {
        throw new RuntimeException("Not yet implemented!");
    }

    @Override
    public FrequencyModel getFrequencyModel() {
        return null;
    }

    @Override
    public void getInfinitesimalMatrix(double[] matrix) {

    }

    @Override
    public DataType getDataType() {
        return null;
    }

    @Override
    public boolean canReturnComplexDiagonalization() {
        return false;
    }

    @Override
    public void addModelListener(ModelListener listener) {

    }

    @Override
    public void removeModelListener(ModelListener listener) {

    }

    @Override
    public boolean isValidState() {
        return false;
    }

    @Override
    public int getModelCount() {
        return 0;
    }


    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {

    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {

    }


    @Override
    protected void storeState() {

    }

    @Override
    protected void restoreState() {

    }

    @Override
    protected void acceptState() {

    }

    @Override
    public boolean isUsed() {
        return false;
    }

    @Override
    public Citation.Category getCategory() {
        return null;
    }

    @Override
    public String getDescription() {
        return "";
    }

    @Override
    public List<Citation> getCitations() {
        return null;
    }

    @Override
    public String getId() {
        return "";
    }

    @Override
    public void setId(String id) {

    }
}
