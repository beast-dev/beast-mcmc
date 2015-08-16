/*
 * LogisticGrowth.java
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

package dr.evolution.coalescent;

/**
 * This class models logistic growth.
 *
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @version $Id: LogisticGrowth.java,v 1.15 2005/05/24 20:25:56 rambaut Exp $
 */
public class LogisticGrowth extends ExponentialGrowth {

    /**
     * Construct demographic model with default settings
     */
    public LogisticGrowth(Type units) {

        super(units);
    }

    public void setShape(double value) {
        c = value;
    }

    public double getShape() {
        return c;
    }

    double lowLimit = 0; // 1e-6;

    /**
     * An alternative parameterization of this model. This
     * function sets the time at which there is a 0.5 proportion
     * of N0.
     *
     * The general form for any k where tk is the time at which Nt = N0/k:
     *		c = (k - 1.0) / (exp(r * tk) - k);
     */
    public void setTime50(double time50) {
        c = 1.0 / (Math.exp(getGrowthRate() * time50) - 2.0);
    }

    public void setShapeFromTimeAtAlpha(double time, double alpha) {

        // New parameterization of logistic shape to be the time at which the
        // population reached some proportion alpha:
        double ert = Math.exp(- getGrowthRate() * time);
        c = ((1.0 - alpha) * ert) / (ert - alpha);
    }
    // Implementation of abstract methods

    /**
     * Gets the value of the demographic function N(t) at time t.
     *
     * @param t the time
     * @return the value of the demographic function N(t) at time t.
     */
    public double getDemographic(double t) {

        double nZero = getN0();
        double r = getGrowthRate();
        double c = getShape();

//		return nZero * (1 + c) / (1 + (c * Math.exp(r*t)));
//		AER rearranging this to use exp(-rt) may help
// 		with some overflow situations...

        double expOfMRT = Math.exp(-r * t);
        return lowLimit + (nZero * (1 + c) * expOfMRT) / (c + expOfMRT);
    }

    public double getLogDemographic(double t) {
        final double d = getDemographic(t);
        if( d == 0.0 && lowLimit == 0.0 ) {
            double nZero = getN0();
            double r = getGrowthRate();
            double c = getShape();
            int sign = c > 0 ? 1 : -1;

            final double v1 = Math.log(c * sign) + r * t;
            double ld = Math.log(nZero);
            if( v1 < 600 ) {
                double v = sign * Math.exp(v1);

                if( c > -1 ) {
                    ld += Math.log1p(c) - Math.log1p(v);
                } else {
                    ld += Math.log((1+c)/(1+v));
                }
            } else {
                 ld += Math.log1p(c) - sign * v1;
            }
            return ld;
        }
//        if(  ! (Math.abs(Math.log(d) - ld) < 1e-12) ) {
//           return Math.log(d);
//        }
        return Math.log(d);
    }
    /**
     * Returns value of demographic intensity function at time t
     * (= integral 1/N(x) dx from 0 to t).
     */
    public double getIntensity(double t) {

        double nZero = getN0();
        double r = getGrowthRate();
        double c = getShape();

        double ert = Math.exp(r * t);
        if( lowLimit == 0 ) {
       // double emrt = Math.exp(-r * t);
          return (c * (ert - 1)/r + t)  / ((1+c) * nZero);
        }
        double z = lowLimit;
        return (r*t*z + (1 + c)*nZero*Math.log(nZero + c*nZero + z + c*ert*z))/(r*z*(nZero + c*nZero + z));
    }

    /**
     * Returns value of demographic intensity function at time t
     * (= integral 1/N(x) dx from 0 to t).
     */
    public double getInverseIntensity(double x) {

        throw new RuntimeException("Not implemented!");
    }

    public double getIntegral(double start, double finish) {
        if( lowLimit > 0 ) {
            double v1 = getNumericalIntegral(start, finish);
            final double v2 = getIntensity(finish) - getIntensity(start);
            return v2;
        }
        double intervalLength = finish - start;

        double nZero = getN0();
        double r = getGrowthRate();
        double c = getShape();
        double expOfMinusRT = Math.exp(-r * start);
        double expOfMinusRG = Math.exp(-r * intervalLength);

        double term1 = nZero * (1.0 + c);
        if (term1 == 0.0) {
            return Double.POSITIVE_INFINITY;
        }

        double term2 = c * (1.0 - expOfMinusRG);

        double term3 = (term1 * expOfMinusRT) * r * expOfMinusRG;
        double term2over3;
        if (term3 == 0.0) {
            double l1 = expOfMinusRG < 1e-8 ?  -r * intervalLength : Math.log1p(expOfMinusRG);
            final int sign = c > 0 ? 1 : -1;
            term2over3 = (sign/term1) * Math.exp(l1 + r * start + Math.log(c*sign) - Math.log(r));

           // throw new RuntimeException("Infinite integral!");
        } else {
//            if (term3 != 0.0 && term2 == 0.0) {
//                term2over3 = 0.0;
//            } else if (term3 == 0.0 && term2 == 0.0) {
//                throw new RuntimeException("term3 and term2 are both zeros. N0=" + getN0() + " growthRate=" + getGrowthRate() + "c=" + c);
//            } else {
                term2over3 = term2 / term3;
//            }
        }

        final double term5 = intervalLength / term1;

//        double v0 = 1/term1 * (finish + c * (Math.exp(r*finish) - 1) /r);
//        double v1 = 1/term1 * (start + c * (Math.exp(r*start) - 1) /r);
//        double v =  1/term1 * ((finish + c * (Math.exp(r*finish) - 1) /r)  - (start + c * (Math.exp(r*start) - 1) /r));
//        double v2 = 1/term1 * ((finish-start) + (c * (Math.exp(r*finish) - 1) /r)  - (c * (Math.exp(r*start) - 1) /r));
//        double v3 = 1/term1 * ((finish-start) + (c/r)* ( Math.exp(r*finish) - 1)   - (c/r) * (Math.exp(r*start) - 1) );
//        double v4 =  1/term1 * ((finish-start) + (c/r) * (Math.exp(r*finish) - Math.exp(r*start) ) );
       // double v = ( (c * (Math.exp(r*finish) - Math.exp(r*start)) / r) + (start - finish)) / term1;

        return term5 + term2over3;
    }

    public int getNumArguments() {
        return 3;
    }

    public String getArgumentName(int n) {
        switch (n) {
            case 0:
                return "N0";
            case 1:
                return "r";
            case 2:
                return "c";
        }
        throw new IllegalArgumentException("Argument " + n + " does not exist");
    }

    public double getArgument(int n) {
        switch (n) {
            case 0:
                return getN0();
            case 1:
                return getGrowthRate();
            case 2:
                return getShape();
        }
        throw new IllegalArgumentException("Argument " + n + " does not exist");
    }

    public void setArgument(int n, double value) {
        switch (n) {
            case 0:
                setN0(value);
                break;
            case 1:
                setGrowthRate(value);
                break;
            case 2:
                setShape(value);
                break;
            default:
                throw new IllegalArgumentException("Argument " + n + " does not exist");

        }
    }

    public double getLowerBound(int n) {
        return 0.0;
    }

    public double getUpperBound(int n) {
        return Double.POSITIVE_INFINITY;
    }

    public DemographicFunction getCopy() {
        LogisticGrowth df = new LogisticGrowth(getUnits());
        df.setN0(getN0());
        df.setGrowthRate(getGrowthRate());
        df.c = c;

        return df;
    }

    //
    // private stuff
    //

    private double c;
}
