package dr.inference.operators;

import dr.inference.model.Parameter;
import dr.inference.model.Bounds;

/**
 * A generic interface for objects capabale of scaling.
 *
 * A default impelementation for any parameter.
 *
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @version $Id: UpDownOperator.java,v 1.25 2005/06/14 10:40:34 rambaut Exp $
 */

public interface Scalable {
    /**
     *
     * @param factor  scaling factor
     * @return  Number of dimentions.
     *
     * @throws OperatorFailedException
     */
    int scale(double factor) throws OperatorFailedException;

    /**
     *
     * @return  Name for display purposes.
     */
    String getName();

    public class Default implements Scalable {
        private Parameter parameter;

        Default(Parameter p) {
            this.parameter = p;
        }

        public int scale(double factor) throws OperatorFailedException {
            final int dimension = parameter.getDimension();

            for(int i = 0; i < dimension; i++) {
                parameter.setParameterValue(i, parameter.getParameterValue(i) * factor);
            }

            final Bounds bounds = parameter.getBounds();
            
            for(int i = 0; i < dimension; i++) {
                final double value = parameter.getParameterValue(i);
                if( value < bounds.getLowerLimit(i) || value > bounds.getUpperLimit(i) ) {
                    throw new OperatorFailedException("proposed value outside boundaries");
                }
            }
            return dimension;
        }

        public String getName() {
            return parameter.getParameterName();
        }
    }
}
