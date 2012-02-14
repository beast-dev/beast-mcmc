package dr.inference.operators;

import dr.math.MathUtils;

import java.util.List;
import java.util.Vector;

/**
 * Package: CombinedOperatorSchedule
 * Description:
 * <p/>
 * <p/>
 * Created by
 *
 * @author Alexander V. Alekseyenko (alexander.alekseyenko@gmail.com)
 *         Date: Dec 23, 2009
 *         Time: 10:17:28 AM
 */
public class CombinedOperatorSchedule implements OperatorSchedule {

    List<OperatorSchedule> operatorSchedules = null;

    public CombinedOperatorSchedule() {
        operatorSchedules = new Vector<OperatorSchedule>();
    }

    public void addOperatorSchedule(OperatorSchedule os) {
        operatorSchedules.add(os);
    }

    public int getScheduleCount() {
        return operatorSchedules.size();
    }

    public int getNextOperatorIndex() {

        final int v = MathUtils.nextInt(operatorSchedules.size());

        int offset = 0;
        for (int i = 0; i < v; ++i) {
            offset += operatorSchedules.get(i).getOperatorCount();
        }

        return offset + operatorSchedules.get(v).getNextOperatorIndex();
    }

    public void reset() {
        for (OperatorSchedule os : operatorSchedules) {
            for (int i = 0; i < os.getOperatorCount(); ++i) {
                os.getOperator(i).reset();
            }
        }
    }

    public int getOperatorCount() {
        int operatorCount = 0;
        for (OperatorSchedule os : operatorSchedules) {
            operatorCount += os.getOperatorCount();
        }
        return operatorCount;
    }

    public MCMCOperator getOperator(int index) {
        for (OperatorSchedule os : operatorSchedules) {
            int opCount = os.getOperatorCount();
            if (index < opCount) {
                return os.getOperator(index);
            } else {
                index -= opCount;
            }
        }
        // if we reach here the index must be out of bounds return null
        return null;
    }

    public void addOperator(MCMCOperator op) {
        System.err.println("ERROR: addOperator() should not be called on CombinedOperatorSchedule, " +
                "use individual operator schedules to add operators. This call will have no affect!");
    }

    public void addOperators(List<MCMCOperator> v) {
        System.err.println("ERROR: addOperators() should not be called on CombinedOperatorSchedule, " +
                "use individual operator schedules to add operators. This call will have no affect!");
    }

    public void operatorsHasBeenUpdated() {
        for (OperatorSchedule os : operatorSchedules) {
            os.operatorsHasBeenUpdated();
        }
    }

    public double getOptimizationTransform(double d) {
        if (operatorSchedules.size() > 0)
            return operatorSchedules.get(0).getOptimizationTransform(d);
        else
            return 0;
    }

    public int getMinimumAcceptAndRejectCount() {
        int minCount = Integer.MAX_VALUE;
        for (OperatorSchedule os : operatorSchedules) {
            if (os.getMinimumAcceptAndRejectCount() < minCount) {
                minCount = os.getMinimumAcceptAndRejectCount();
            }
        }
        return minCount;
    }


}
