/*
 * PercentColumn.java
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

package dr.inference.loggers;

/**
 *
 * Percent column - multiplies by 100 and prepends a '%'.
 *
 * Values outside [0,1] are shown as is.
 *
 * @author Joseph Heled
 *         Date: 4/06/2008
 */
public class PercentColumn extends NumberColumn {
    private final NumberColumn column;

    public PercentColumn(NumberColumn col) {
        super(col.getLabel());
        this.column = col;
    }

    public void setSignificantFigures(int sf) {
        column.setSignificantFigures(sf);
    }
    
    public int getSignificantFigures() {
        return column.getSignificantFigures();
    }

    public void setMinimumWidth(int minimumWidth) {
      column.setMinimumWidth(minimumWidth);
    }

    public int getMinimumWidth() {
        return column.getMinimumWidth();
    }

    public String getFormattedValue() {
        double val = column.getDoubleValue();
        if( val >= 0 && val <= 1 ) {
            return column.formatValue(val * 100) + "%";
        }
        return column.getFormattedValue();
    }

    public double getDoubleValue() {
        return column.getDoubleValue();
    }
}
