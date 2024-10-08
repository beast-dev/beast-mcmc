/*
 * OperatorSchedule.java
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

package dr.inference.operators;

import java.io.Serializable;
import java.util.List;

/**
 * An interface the defines an operator schedule for use in
 * choosing the next operator during an MCMC run.
 *
 * @author Alexei Drummond
 */
public interface OperatorSchedule extends Serializable {

    String SHOW_OPERATORS = "show_operators";
    OptimizationTransform DEFAULT_TRANSFORM = OptimizationTransform.LOG;
    //int POWERB = 1000000; //TODO discuss with Luiz why this value doesn't seem to work
    int POWERB = 1;
    double POWERC = 0.5556; // c = 5/9

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
    OptimizationTransform getOptimizationTransform();

    /**
     * @return the minimum number of times an operator has been called
     */
    long getMinimumAcceptAndRejectCount();

    enum OptimizationTransform {

        LOG("log") {
            @Override public double transform(double d) { return Math.log(d); }
        },
        SQRT("sqrt") {
            @Override public double transform(double d) { return Math.sqrt(d); }
        },
        LINEAR("linear") {
            @Override public double transform(double d) { return d; }
        },
        POWER("power") {
            @Override public double transform(double d) {
                return Math.pow(d/((double)POWERB), POWERC);
            }
        };

        OptimizationTransform(String name) {
            this.name = name;
        }

        public abstract double transform(double d);

        @Override
        public String toString() {
            return name;
        }

        private final String name;
    }
}
