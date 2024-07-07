package dr.app.tools.newtreeannotator;

import dr.evolution.util.Taxon;

import java.util.List;

/**
 * @author Andrew Rambaut
 * @version $
 */
interface Clade {
    int getCount();

    void setCount(int count);

    double getCredibility();

    void setCredibility(double credibility);

    int getSize();

    Taxon getTaxon();

    Clade getBestLeft();

    Clade getBestRight();

    Object getKey();

    void addAttributeValues(Object[] values);

    List<Object[]> getAttributeValues();
}
