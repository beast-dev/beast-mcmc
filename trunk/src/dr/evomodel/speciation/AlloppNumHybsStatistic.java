package dr.evomodel.speciation;

import dr.inference.model.Statistic;

/**
 * Created with IntelliJ IDEA.
 * User: Graham
 * Date: 08/10/12
 * Time: 18:24
 * To change this template use File | Settings | File Templates.
 */
public class AlloppNumHybsStatistic  extends Statistic.Abstract {
    AlloppSpeciesNetworkModel aspnet;

    public AlloppNumHybsStatistic(AlloppSpeciesNetworkModel aspnet) {
        super("NumHybs");
        this.aspnet = aspnet;
    }


    @Override
    public int getDimension() {
        return 1;
    }

    @Override
    public double getStatisticValue(int dim) {
        return aspnet.getNumberOfTetraTrees();
    }
}
