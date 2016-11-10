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
    protected double getLogPDF(double value){
        if(value > low.getParameterValue(0) && value < high.getParameterValue(0)){
            double p1 = getDistribution().logPdf(value);
//            System.out.println(p1);
            double p2 = 0;
            if(!Double.isInfinite(low.getParameterValue(0)))
                p2 = getDistribution().cdf(low.getParameterValue(0));
//            System.out.println(p2);
            double p3 = 1;
            if(!Double.isInfinite(high.getParameterValue(0)))
                p3 = 1 - getDistribution().cdf(high.getParameterValue(0));
//            System.out.println(p3);

            return p1 - Math.log(p3 + p2);
        }
        else
            return Double.NEGATIVE_INFINITY;
    }

    Parameter low;
    Parameter high;
}
