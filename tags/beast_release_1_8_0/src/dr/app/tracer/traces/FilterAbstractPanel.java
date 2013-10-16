package dr.app.tracer.traces;

import javax.swing.*;

/**
 * @author Walter Xie
 */
public abstract class FilterAbstractPanel extends JPanel {

    abstract String[] getSelectedValues();

    boolean containsNullValue() {
        if (getSelectedValues() == null || getSelectedValues().length < 1) return true;
        for (String ob : getSelectedValues()) {
            if (ob == null || ob.equals("")) return true;
        }
        return false;
    }
//
//    FilterAbstractPanel getInstance(String[] a, String[] b, TraceFactory.TraceType traceType) {
//        if (traceType == TraceFactory.TraceType.DOUBLE) {
//            return new FilterContinuousPanel(a, b);
//        } else {
//            return new FilterDiscretePanel(a, b);
//        }
//    }
}
