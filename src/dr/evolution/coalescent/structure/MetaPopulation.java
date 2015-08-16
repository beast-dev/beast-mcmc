/*
 * MetaPopulation.java
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

package dr.evolution.coalescent.structure;

/**
 * Created by IntelliJ IDEA.
 * User: alexei
 * Date: Jul 13, 2009
 * Time: 9:52:38 PM
 * To change this template use File | Settings | File Templates.
 */
public interface MetaPopulation {
    int getPopulationCount();

    double[] getPopulationSizes(double time);/* returns value of demographic function at time t  (population size; one entry of double[] getPopulationSizes)
    * (This function mirrors an equivalent function in DemographicFunction)
    */

    double getDemographic(double time, int population);/* calculates the integral 1/N(x) dx from start to finish, for one of the populations
    * (This function mirrors an equivalent function in DemographicFunction)
    */

    double getIntegral(double start, double finish, int population);
}
