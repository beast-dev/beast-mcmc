/*
 * MutationDeathModel.java
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

/**
 * Package: MutationDeathModel
 * Description:
 *
 *
 * Created by
 * Alexander V. Alekseyenko (alexander.alekseyenko@gmail.com)
 * Date: Feb 19, 2008
 * Time: 12:41:01 PM
 */

package dr.evomodel.substmodel;

import dr.evolution.datatype.DataType;
import dr.oldevomodelxml.substmodel.MutationDeathModelParser;
import dr.inference.model.Model;
import dr.inference.model.Parameter;

import java.util.Arrays;

public class MutationDeathModel extends ComplexSubstitutionModel { //BaseSubstitutionModel {
    private SubstitutionModel CTMCModel;
    private Parameter delParameter = null;
    protected double[] trMatrix;
    private Parameter baseSubModelFreq;
    private Parameter thisSubModelFreq;
    private Parameter mutationRate;

    private double[] baseModelMatrix;


    public MutationDeathModel(Parameter delParameter, DataType dT, SubstitutionModel evoModel,
                              FrequencyModel freqModel, Parameter mutationRate) {
        super(MutationDeathModelParser.MD_MODEL, dT, freqModel, null);

        CTMCModel = evoModel;
        stateCount = freqModel.getFrequencyCount();
        this.delParameter = delParameter;
        this.dataType = dT;
        this.mutationRate = mutationRate;

        addVariable(delParameter);
        delParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));

        addVariable(mutationRate);
//        addModel(freqModel);
        if (evoModel != null) {
            addModel(evoModel.getFrequencyModel());
            addModel(evoModel);
        }
        trMatrix = new double[(stateCount - 1) * (stateCount - 1)];

        if (CTMCModel != null) {
            baseSubModelFreq = CTMCModel.getFrequencyModel().getFrequencyParameter();
        } else {
            baseSubModelFreq = new Parameter.Default(new double[]{1.0});
        }
        thisSubModelFreq = getFrequencyModel().getFrequencyParameter();

        double total = 0;
        for (int i = 0; i < baseSubModelFreq.getDimension(); i++) {
            double value = thisSubModelFreq.getParameterValue(i);
            total += value;
            baseSubModelFreq.setParameterValue(i, value);
        }

        for (int i = 0; i < baseSubModelFreq.getDimension(); i++) {
            baseSubModelFreq.setParameterValue(i, baseSubModelFreq.getParameterValue(i) / total);
        }

        thisSubModelFreq.setParameterValue(thisSubModelFreq.getDimension() - 1, 0.0);

        copyFrequencies();

        frequenciesChanged();
        ratesChanged();

        if (CTMCModel != null) {
            final int baseStateCount = CTMCModel.getDataType().getStateCount();
            baseModelMatrix = new double[baseStateCount * baseStateCount];
        } else {
            baseModelMatrix = new double[1];
        }
        setDoNormalization(false);

//        setupRelativeRates();
    }

    public Parameter getDeathParameter() {
        return this.delParameter;
    }

    protected void handleModelChangedEvent(Model model, Object object, int index) {

        if (object == baseSubModelFreq) {
            copyFrequencies();
            fireModelChanged(object, index);
        } else if (model == CTMCModel) {
            fireModelChanged(object, index);
        }
    }

    private void copyFrequencies() {
        for (int i = 0; i < baseSubModelFreq.getDimension(); i++)
            thisSubModelFreq.setParameterValueQuietly(i, baseSubModelFreq.getParameterValue(i));

    }

    @Override
    protected void frequenciesChanged() {
    }

    @Override
    protected void ratesChanged() {
    }

    @Override
    protected void setupQMatrix(double[] rates, double[] pi, double[][] matrix) {

        // Zero matrix
        for (int i = 0; i < matrix.length; ++i) {
            Arrays.fill(matrix[i], 0.0);
        }

        int baseStateCount = 1;
        if (CTMCModel != null) {
            CTMCModel.getInfinitesimalMatrix(baseModelMatrix);
            baseStateCount = CTMCModel.getDataType().getStateCount();
        }

        double deathR = delParameter.getParameterValue(0);
        double mutationR = 2 * mutationRate.getParameterValue(0); // TODO Double-check

        int k = 0;
        for (int i = 0; i < baseStateCount; ++i) {
            for (int j = 0; j < baseStateCount; ++j) {
                matrix[i][j] = baseModelMatrix[k] - deathR;
                k++;
            }
            matrix[i][baseStateCount] = deathR;
        }
    }

    @Override
    protected void setupRelativeRates(double[] rates) {

    }

}