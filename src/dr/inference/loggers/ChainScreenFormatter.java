/*
 * ChainScreenFormatter.java
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

package dr.inference.loggers;

/**
 * The destination of the screen log of one chain when it is shown with convergence diagnostics
 * or when several chains are run together and their screen logs are shown as one.
 *
 * @author Andrew Rambaut
 */
public interface ChainScreenFormatter extends LogFormatter {

    /**
     * @return true if convergence diagnostics are being computed for the run, so
     * {@link ChainDiagnosticColumn}s can be shown
     */
    boolean hasDiagnostics();

    /**
     * @return the number of independent chains being run
     */
    int getChainCount();

    /**
     * Gives the columns of the chain's screen log, in order (without the state).
     */
    void setColumns(LogColumn[] columns);
}
