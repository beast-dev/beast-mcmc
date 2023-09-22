package dr.evomodel.siteratemodel;

import dr.inference.model.Model;

/**
 * @author Andrew Rambaut
 * @version $
 */
public interface SiteRateDelegate extends Model {
     int getCategoryCount();
     
     void getCategories(double[] categoryRates, double[] categoryProportions);
}
