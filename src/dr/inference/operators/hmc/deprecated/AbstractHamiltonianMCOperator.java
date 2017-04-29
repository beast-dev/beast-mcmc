/*
 * AbstractHamiltonianMCOperator.java
 *
 * Copyright (c) 2002-2017 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.inference.operators.hmc.deprecated;

import dr.inference.operators.AbstractCoercableOperator;
import dr.inference.operators.CoercionMode;
import dr.math.distributions.NormalDistribution;

/**
 * Created by max on 12/3/15.
 */
@Deprecated
public abstract class AbstractHamiltonianMCOperator extends AbstractCoercableOperator {
    public AbstractHamiltonianMCOperator(CoercionMode mode, double momentumSd) {
        super(mode);
        this.momentumSd=momentumSd;
    }

    protected double getMomentumSd()
    {return momentumSd;}

    protected void setMomentumSd(double momentum){
        momentumSd=momentum;
    }

    private double momentumSd;
    protected double[] momentum;

    protected void drawMomentum(int size){
        momentum=new double[size];
        for (int i = 0; i <size ; i++) {
            momentum[i]= (Double) (new NormalDistribution(0.0, momentumSd)).nextRandom();
        }
    }
}
