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
