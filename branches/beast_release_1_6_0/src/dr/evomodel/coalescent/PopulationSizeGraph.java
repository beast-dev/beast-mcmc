package dr.evomodel.coalescent;

import dr.inference.model.Statistic;

/**
 * @author Joseph Heled
 * @version $Id$
 *          Created by IntelliJ IDEA.
 *          User: joseph Date: 5/02/2007 Time: 16:17:11
 */
public class PopulationSizeGraph extends Statistic.Abstract {

    private double tm = 0;
    private VariableDemographicModel vdm = null;

    public PopulationSizeGraph(VariableDemographicModel vdm, double  tm) {
        super("popGraph");
        this.vdm = vdm;
        this.tm = tm;
    }

    public int getDimension() { return 1; }

    public double getStatisticValue(int dim) {
        return vdm.getDemographicFunction().getDemographic(tm);
    }
}
