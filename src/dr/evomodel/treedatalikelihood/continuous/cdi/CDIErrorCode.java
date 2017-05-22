/*
 * CDIErrorCode.java
 *
 * Copyright (c) 2002-2016 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodel.treedatalikelihood.continuous.cdi;

/**
 * Created by msuchard on 9/14/16.
 */
public enum CDIErrorCode {

    NO_ERROR(0, "no error"),
    GENERAL_ERROR(-1, "general error"),
    OUT_OF_MEMORY_ERROR(-2, "out of memory error"),
    UNIDENTIFIED_EXCEPTION_ERROR(-3, "unidentified exception error"),
    UNINITIALIZED_INSTANCE_ERROR(-4, "uninitialized instance error"),
    OUT_OF_RANGE_ERROR(-5, "One of the indices specified exceeded the range of the array"),
    NO_RESOURCE_ERROR(-6, "No resource matches requirements"),
    NO_IMPLEMENTATION_ERROR(-7, "No implementation matches requirements"),
    FLOATING_POINT_ERROR(-8, "Floating-point range exceeded");

    CDIErrorCode(int errCode, String meaning) {
        this.errCode = errCode;
        this.meaning = meaning;
    }

    public int getErrCode() {
        return errCode;
    }

    public String getMeaning() {
        return meaning;
    }

    private final int errCode;
    private final String meaning;
}
