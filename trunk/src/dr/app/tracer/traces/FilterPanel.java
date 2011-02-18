package dr.app.tracer.traces;

import javax.swing.*;

/**
 * @author Walter Xie
 */
public abstract class FilterPanel extends JPanel {

    abstract Object[] getSelectedValues();

    boolean containsNullValue() {
        if (getSelectedValues() == null) return true;
        for (Object ob : getSelectedValues()) {
            if (ob == null || ob.toString().equals("")) return true;
        }
        return false;
    }
//
//    FilterPanel getInstance(String[] a, String[] b, TraceFactory.TraceType traceType) {
//        if (traceType == TraceFactory.TraceType.DOUBLE) {
//            return new FilterContinuousPanel(a, b);
//        } else {
//            return new FilterDiscretePanel(a, b);
//        }
//    }
}
