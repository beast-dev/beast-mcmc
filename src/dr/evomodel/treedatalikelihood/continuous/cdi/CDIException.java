/*
 * CDIException.java
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
public class CDIException extends RuntimeException {
    public CDIException(String functionName, int errCode) {
        this.functionName = functionName;
        this.errCode = errCode;
    }

    @Override
    public String getMessage() {
        CDIErrorCode code = null;
        for (CDIErrorCode c : CDIErrorCode.values()) {
            if (c.getErrCode() == errCode) {
                code = c;
            }
        }
        if (code == null) {
            return "CDI function, " + functionName + ", returned error code " + errCode + " (unrecognized error code)";
        }
        return "CDI function, " + functionName + ", returned error code " + errCode + " (" + code.getMeaning() + ")";
    }

    public String getFunctionName() {
        return functionName;
    }

    public int getErrCode() {
        return errCode;
    }

    private final String functionName;
    private final int errCode;
}
