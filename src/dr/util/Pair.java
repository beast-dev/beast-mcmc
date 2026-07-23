/*
 * Pair.java
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

package dr.util;

/**
 * Date: 15/06/2016
 * Time: 17:57
 *
 * @author rambaut
 */
public class Pair<A, B> {
    public final A first;
    public final B second;

    public Pair(A first, B second) {
        this.first = first;
        this.second = second;
    }

    public A getFirst() {
        return first;
    }

    public B getSecond() {
        return second;
    }

    public String toString() {
        return "Pair[" + this.first + "," + this.second + "]";
    }

    private boolean equals(Object var0, Object var1) {
        return var0 == null && var1 == null || var0 != null && var0.equals(var1);
    }

    public boolean equals(Object var1) {
        return var1 instanceof Pair && equals(this.first, ((Pair) var1).first) && equals(this.second, ((Pair) var1).second);
    }

    public int hashCode() {
        return this.first == null ? (this.second == null ? 0 : this.second.hashCode() + 1) : (this.second == null ? this.first.hashCode() + 2 : this.first.hashCode() * 17 + this.second.hashCode());
    }
}
