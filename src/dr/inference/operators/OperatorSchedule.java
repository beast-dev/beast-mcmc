/*
 * OperatorSchedule.java
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

package dr.inference.operators;

import java.io.Serializable;
import java.util.List;

/**
 * An interface the defines an operator schedule for use in
 * choosing the next operator during an MCMC run.
 *
 * @author Alexei Drummond
 * @version $Id: OperatorSchedule.java,v 1.3 2005/05/24 20:26:00 rambaut Exp $
 */
public interface OperatorSchedule extends Serializable {
    /**
     * @return Choose the next operator.
     */
    int getNextOperatorIndex();

    /**
     * @return Total number of operators
     */
    int getOperatorCount();

    /**
     * @param index
     * @return the index'th operator
     */
    MCMCOperator getOperator(int index);

    void addOperator(MCMCOperator op);

    void addOperators(List<MCMCOperator> v);

    /**
     * Should be called after operators weight is externally changed.
     */
    void operatorsHasBeenUpdated();

    /**
     * @return the optimization schedule
     */
    double getOptimizationTransform(double d);

    /**
     * @return the minimum number of times an operator has been called
     */
    long getMinimumAcceptAndRejectCount();

    public enum OptimizationTransform {
        DEFAULT("default"),
        LOG("log"),
        SQRT("sqrt"),
        LINEAR("linear");

        OptimizationTransform(String name) {
            this.name = name;
        }

        @Override
        public String toString() {
            return name;
        }

        private final String name;
    };
}
