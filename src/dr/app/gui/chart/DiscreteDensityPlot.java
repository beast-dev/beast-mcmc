/*
 * CategoryDensityPlot.java
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

package dr.app.gui.chart;

import dr.inference.trace.TraceDistribution;
import dr.stats.FrequencyCounter;
import dr.stats.Variate;
import dr.util.FrequencyDistribution;

import java.awt.*;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

public class DiscreteDensityPlot extends ColumnPlot {
    public DiscreteDensityPlot(FrequencyCounter<Integer> frequencies) {
        super();

        setData(frequencies);
    }

    /**
     * Set data, all integers
     */
    private void setData(FrequencyCounter<Integer> frequencies) {
        Variate.D xData = new Variate.D();
        Variate.D yData = new Variate.D();

        for (int i : frequencies.getUniqueValues()) {
            int index = i;

//            if (orderMap != null && orderMap.size() > 0) {
//                if (i >= orderMap.size()) {
//                    System.out.println("oops");
//                }
//                index = orderMap.get(i);
//            }
            
            xData.add((double)index);
            yData.add(frequencies.getProbability(index));

        }

        setData(xData, yData);
    }
}


