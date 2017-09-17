/*
 * CorrelationData.java
 *
 * Copyright (c) 2002-2017 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.app.gui.util;

import java.util.HashMap;
import java.util.List;
import java.util.Set;

/**
 * @author Guy Baele
 */
public class CorrelationData {

    private HashMap<String, List<Double>> data;

    public CorrelationData() {
        data = new HashMap<String, List<Double>>();
    }

    public void add(String traceName, List values) {
        if (data.containsKey(traceName)) {
            List temp = data.get(traceName);
            temp.addAll(values);
            data.put(traceName, temp);
        } else {
            data.put(traceName, values);
        }
    }

    public Set<String> getTraceNames() {
        return data.keySet();
    }

    public List getDataForKey(String traceName) {
        return data.get(traceName);
    }

    public int numberOfEntries() {
        return data.size();
    }

    public void clear() {
        this.data.clear();
    }

}
