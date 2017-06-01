/*
 * MsatPatternStatistic.java
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

package dr.evolution.alignment;


/**
 * @author Chieh-Hsi
 */
public class MicrosatellitePatternStatistic {
    Patterns msatPattern;
    double msatLengthVar;
    String mode;


    public MicrosatellitePatternStatistic(Patterns msatPattern){
        this(msatPattern, "variance");
    }

    public MicrosatellitePatternStatistic(Patterns msatPattern, String mode){
        this.msatPattern = msatPattern;
        this.msatLengthVar = computeMsatLengthVariance();
        this.mode = mode;
    }

    public double computeMsatLengthVariance(){
        double var = 0.0;
        int[] msatPat = msatPattern.getPattern(0);
        double mean = 0.0;
        for(int i = 0; i < msatPat.length; i++){
            mean += msatPat[i];
        }
        mean = mean/msatPat.length;
        for(int i = 0; i < msatPat.length; i++){
            var+=(msatPat[i] - mean)*(msatPat[i] - mean);
        }
        var = var/msatPat.length;
        System.out.println(2*var);
        return var;

    }

    public String toString(){
        if(mode.equals("thetaV")){
            return ""+2*msatLengthVar;
        }else{
            return ""+msatLengthVar;
        }

    }


}
