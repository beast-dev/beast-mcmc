/*
 * Analysis.java
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

package dr.app.mapper.application.mapper;

import dr.app.mapper.application.MapperDocument;

import java.util.ArrayList;
import java.util.List;

/**
 * Date: 09/03/2013
 * Time: 09:16
 *
 * @author rambaut
 */
public class Analysis {
    public void addDataSet(MapperDocument dataSet) {
        dataSets.add(dataSet);
        fireDataSetChanged();
    }

    public int getDataSetCount() {
        return dataSets.size();
    }

    public MapperDocument getDataSet(int index) {
        return dataSets.get(index);
    }

    /**
     * Return a copy of the list of data sets.
     * @return the data sets
     */
    public List<MapperDocument> getDataSets() {
        return new ArrayList<MapperDocument>(dataSets);
    }

    // Listeners and broadcasting
    public void addListener(AnalysisListener listener) {
        listeners.add(listener);
    }

    public void removeListener(AnalysisListener listener) {
        listeners.remove(listener);
    }

    private void fireDataSetChanged() {
      for (AnalysisListener listener : listeners) {
          listener.analysisChanged();
      }
    }

    private final List<AnalysisListener> listeners = new ArrayList<AnalysisListener>();
    private final List<MapperDocument> dataSets = new ArrayList<MapperDocument>();
}
