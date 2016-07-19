/*
 * ProductChainFrequencyModel.java
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

import dr.inference.model.Model;
import dr.inference.model.Parameter;

import java.util.List;

/**
 * @author Marc A. Suchard
 * @author Vladimir Minin
 *         <p/>
 *         A class for implementing a kronecker sum of CTMC models in BEAST using BEAGLE
 *         This work is supported by NSF grant 0856099
 *         <p/>
 *         O'Brien JD, Minin VN and Suchard MA (2009) Learning to count: robust estimates for labeled distances between
 *         molecular sequences. Molecular Biology and Evolution, 26, 801-814
 */

public class ProductChainFrequencyModel extends FrequencyModel {

    public ProductChainFrequencyModel(String name, List<FrequencyModel> freqModels) {
        super(name);
        this.freqModels = freqModels;
        int freqCount = 1;
        numBaseModel = freqModels.size();
        stateSizes = new int[numBaseModel];
        for (int i = 0; i < numBaseModel; i++) {
            int size = freqModels.get(i).getFrequencyCount();
            stateSizes[i] = size;
            freqCount *= size;
            addModel(freqModels.get(i));
        }
        tmp = new int[numBaseModel];
        totalFreqCount = freqCount;
    }

    protected void handleModelChangedEvent(Model model, Object object, int index) {
        fireModelChanged(model);
    }

    public void setFrequency(int i, double value) {
        throw new RuntimeException("Not implemented");
    }

    public double getFrequency(int index) {
        double freq = 1.0;
        decomposeEntry(index, tmp);       
        for (int i = 0; i < numBaseModel; i++) {
            freq *= freqModels.get(i).getFrequency(tmp[i]);
        }
        return freq;
    }

    public int[] decomposeEntry(int index) {
        int[] tmp = new int[numBaseModel];
        decomposeEntry(index, tmp);
        return tmp;
    }
    
    private void decomposeEntry(int index, int[] decomposition) {
        int current = index;
        for (int i = numBaseModel - 1; i >= 0; --i) {           
            decomposition[i] = current % stateSizes[i];
            current /= stateSizes[i];
        }
    }

    public int getFrequencyCount() {
        return totalFreqCount;
    }

    public Parameter getFrequencyParameter() {
        throw new RuntimeException("Not implemented");
    }
     
    private List<FrequencyModel> freqModels;

    private final int numBaseModel;
    private final int totalFreqCount;
    private final int[] stateSizes;
    private final int[] tmp;
}
