/*
 * TruncatedDistributionLikelihood.java
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

package dr.inference.model;


import dr.inference.distribution.DistributionLikelihood;
import dr.util.Attribute;

import java.util.ArrayList;

/**
 * @author Max Tolkoff
 */

public class TruncatedDistributionLikelihood extends DistributionLikelihood {
    public TruncatedDistributionLikelihood(DistributionLikelihood distribution, Parameter low, Parameter high){
        super(distribution.getDistribution());
        this.dataList = (ArrayList<Attribute<double[]>>) distribution.getDataList();
//        for(Attribute<double[]> data : list){
//            System.out.println("here");
//            this.addData(data);
//        }
        this.low = low;
        this.high = high;
    }

    @Override
    protected double getLogPDF(double value, int i){
        if(high.getParameterValue(i % high.getDimension()) == low.getParameterValue(i % low.getDimension()))
            return 0;
        if(value > low.getParameterValue(i % low.getDimension()) && value < high.getParameterValue(i % high.getDimension())){
            double p1 = getDistribution().logPdf(value);
//            System.out.println(p1);
            double p2 = 0;
            if(!Double.isInfinite(low.getParameterValue(i % low.getDimension())))
                p2 = getDistribution().cdf(low.getParameterValue(i % low.getDimension()));
//            System.out.println(p2);
            double p3 = 0;
            if(!Double.isInfinite(high.getParameterValue(i % high.getDimension())))
                p3 = 1 - getDistribution().cdf(high.getParameterValue(i % high.getDimension()));
//            System.out.println(p3);
            return p1 - Math.log(1 - (p3 + p2));
        }
        else{
//            System.out.println(high.getParameterValue(i % high.getDimension()));
//            System.out.println(low.getParameterValue(i % low.getDimension()));
            return Double.NEGATIVE_INFINITY;
        }
    }

    Parameter low;
    Parameter high;
}
