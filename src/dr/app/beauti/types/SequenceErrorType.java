/*
 * SequenceErrorType.java
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

package dr.app.beauti.types;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 */
public enum SequenceErrorType {

    NO_ERROR("Off"),
    AGE_TRANSITIONS("Age-dependent: Transitions only"),
    AGE_ALL("Age-dependent: All substitutions"),
    BASE_TRANSITIONS("Age-independent: Transitions only"),
    BASE_ALL("Age-independent: All substitutions"),
    HYPERMUTATION_HA3G("RT Hypermutation: hA3G"),
    HYPERMUTATION_HA3F("RT Hypermutation: hA3F"),
    HYPERMUTATION_BOTH("RT Hypermutation: hA3G + hA3F"),
    HYPERMUTATION_ALL("RT Hypermutation: G->A");

    SequenceErrorType(String displayName) {
        this.displayName = displayName;
    }

    public String toString() {
        return displayName;
    }

    private final String displayName;
}