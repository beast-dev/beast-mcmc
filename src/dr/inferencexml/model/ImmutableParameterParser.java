package dr.inferencexml.model;

import dr.inference.model.Bounds;
import dr.inference.model.Parameter;
import dr.inference.model.Statistic;
import dr.xml.*;

/**
 * Created by Guy Baele on 18/12/15.
 */
public class ImmutableParameterParser extends AbstractXMLObjectParser {

    public static final String IMMUTABLE_PARAMETER = "immutableParameter";

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        final Statistic statistic = (Statistic) xo.getChild(Statistic.class);

        Parameter.Abstract immutableParameter = new Parameter.Abstract() {
            public void setParameterValueNotifyChangedAll(int dim, double value) {
                throw new RuntimeException("Forbidden call to ImmutableParameter.");
            }
            public void setParameterValueQuietly(int dim, double value) {
                throw new RuntimeException("Forbidden call to ImmutableParameter.");
            }
            public void storeValues() {
                //do nothing
            }
            public void restoreValues() {
                //do nothing
            }
            public void acceptValues() {
                //do nothing
            }

            public int getDimension() {
                return statistic.getDimension();
            }

            public void setParameterValue(int dim, double value) {
                throw new RuntimeException("Forbidden call to ImmutableParameter.");
            }
            public double getParameterValue(int dim) {
                return statistic.getStatisticValue(dim);
            }
            public String getParameterName() {
                if (getId() == null)
                    return "immutable." + statistic.getStatisticName();
                return getId();
            }
            public void adoptValues(Parameter source) {
                throw new RuntimeException("Forbidden call to ImmutableParameter.");
            }
            public void addDimension(int index, double value) {
                throw new RuntimeException("Forbidden call to ImmutableParameter.");
            }
            public double removeDimension(int index) {
                throw new RuntimeException("Forbidden call to ImmutableParameter.");
            }
            public void addBounds(Bounds<Double> bounds) {
                throw new RuntimeException("Forbidden call to ImmutableParameter.");
            }
            public Bounds<Double> getBounds() {
                throw new RuntimeException("Forbidden call to ImmutableParameter.");
            }
        };

        return immutableParameter;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(Statistic.class),
    };

    public String getParserDescription() {
        return "An immutable parameter generated from a statistic.";
    }

    public Class getReturnType() {
        return Parameter.class;
    }

    public String getParserName() {
        return IMMUTABLE_PARAMETER;
    }

}
