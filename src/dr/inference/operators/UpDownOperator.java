/*
 * UpDownOperator.java
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

import dr.math.MathUtils;

public class UpDownOperator extends AbstractAdaptableOperator {

    private Scalable[] upParameter = null;
    private Scalable[] downParameter = null;
    private double scaleFactor;

    public UpDownOperator(Scalable[] upParameter, Scalable[] downParameter,
                          double scale, double weight, AdaptationMode mode) {

        super(mode);
        setWeight(weight);

        this.upParameter = upParameter;
        this.downParameter = downParameter;
        this.scaleFactor = scale;
    }

    public final double getScaleFactor() {
        return scaleFactor;
    }

    public final void setScaleFactor(double sf) {
        if( (sf > 0.0) && (sf < 1.0) ) {
            scaleFactor = sf;
        } else {
            throw new IllegalArgumentException("scale must be between 0 and 1");
        }
    }

    /**
     * change the parameter and return the hastings ratio.
     */
    public final double doOperation() {

        final double scale = (scaleFactor + (MathUtils.nextDouble() * ((1.0 / scaleFactor) - scaleFactor)));
        int goingUp = 0, goingDown = 0;

        if( upParameter != null ) {
            for( Scalable up : upParameter ) {
                goingUp += up.scale(scale, -1, false);
//                try {
//                    goingUp += up.scale(scale, -1, true);
//                } catch (RuntimeException re) {
//                    return Double.NEGATIVE_INFINITY;
//                }
            }
            for( Scalable up : upParameter ) {
                if (!up.testBounds()) {
                    return Double.NEGATIVE_INFINITY;
                }
            }
        }

        if( downParameter != null ) {
            for( Scalable dn : downParameter ) {
                goingDown += dn.scale(1.0 / scale, -1, false);
//                try {
//                    goingDown += dn.scale(1.0 / scale, -1, true);
//                } catch (RuntimeException re) {
//                    return Double.NEGATIVE_INFINITY;
//                }
            }
            for( Scalable dn : downParameter ) {
                if (!dn.testBounds()) {
                    return Double.NEGATIVE_INFINITY;
                }
            }
        }


        return (goingUp - goingDown - 2) * Math.log(scale);
    }

    public String getAdaptableParameterName() {
        return "scaleFactor";
    }

    public final String getOperatorName() {
        String name = "";
        if( upParameter != null ) {
            name = "up:";
            for( Scalable up : upParameter ) {
                name = name + up.getName() + " ";
            }
        }

        if( downParameter != null ) {
            name += "down:";
            for( Scalable dn : downParameter ) {
                name = name + dn.getName() + " ";
            }
        }
        return name;
    }

    @Override
    protected double getAdaptableParameterValue() {
        return Math.log(1.0 / scaleFactor - 1.0) / Math.log(10);
    }

    @Override
    public void setAdaptableParameterValue(double value) {
        scaleFactor = 1.0 / (Math.pow(10.0, value) + 1.0);
    }

    @Override
    public double getRawParameter() {
        return scaleFactor;
    }

}


