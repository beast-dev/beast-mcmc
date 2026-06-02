/*
 * XY.java
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

package dr.app.gui.chart;

/**
 * make x and y together as the key of map
 *
 * @author Walter Xie
 */
public class XY implements Comparable<XY> {
    final double x, y;

    XY(double x, double y) {
        this.x = x;
        this.y = y;
    }

    public boolean equals(Object o){
        if(this==o)return true;
        if(o==null)return false;
        if(o instanceof XY){
            XY xy = (XY)o;
            return this.x==xy.x && this.y==xy.y;
        }
        return false;
    }

    public int hashCode(){
        return (int) (x * 31 + y);
    }

    //http://stackoverflow.com/questions/9307751/override-compareto-and-sort-using-two-strings
    @Override
    public int compareTo(final XY o) {
        int cmp = Double.compare(this.x, o.x);
        if (cmp == 0) cmp = Double.compare(this.y, o.y);
        return cmp;
    }

//    public String toString() {
//        return Double.toString(x) + "|" + Double.toString(y);
//    }
}
