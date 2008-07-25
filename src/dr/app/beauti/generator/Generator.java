package dr.app.beauti.generator;

import dr.app.beauti.PriorType;
import dr.app.beauti.XMLWriter;
import dr.app.beauti.options.BeautiOptions;
import dr.app.beauti.options.ModelOptions;
import dr.inference.loggers.Columns;
import dr.inference.model.ParameterParser;
import dr.util.Attribute;

import java.util.ArrayList;

/**
 * @author Alexei Drummond
 */
public abstract class Generator extends ModelOptions {

    BeautiOptions options;

    public Generator(BeautiOptions options) {
        this.options = options;
    }

    /**
     * fix a parameter
     *
     * @param id    the id
     * @param value the value
     */
    public void fixParameter(String id, double value) {
        dr.app.beauti.options.Parameter parameter = parameters.get(id);
        if (parameter == null) {
            throw new IllegalArgumentException("parameter with name, " + id + ", is unknown");
        }
        parameter.isFixed = true;
        parameter.initial = value;
    }

    /**
     * write a parameter
     *
     * @param id     the id
     * @param writer the writer
     */
    public void writeParameter(String id, XMLWriter writer) {
        dr.app.beauti.options.Parameter parameter = parameters.get(id);
        if (parameter == null) {
            throw new IllegalArgumentException("parameter with name, " + id + ", is unknown");
        }
        if (parameter.isFixed) {
            writeParameter(id, 1, parameter.initial, Double.NaN, Double.NaN, writer);
        } else {
            if (parameter.priorType == PriorType.UNIFORM_PRIOR || parameter.priorType == PriorType.TRUNC_NORMAL_PRIOR) {
                writeParameter(id, 1, parameter.initial, parameter.uniformLower, parameter.uniformUpper, writer);
            } else {
                writeParameter(id, 1, parameter.initial, parameter.lower, parameter.upper, writer);
            }
        }
    }

    /**
     * write a parameter
     *
     * @param id        the id
     * @param dimension the dimension
     * @param writer    the writer
     */
    public void writeParameter(String id, int dimension, XMLWriter writer) {
        dr.app.beauti.options.Parameter parameter = parameters.get(id);
        if (parameter == null) {
            throw new IllegalArgumentException("parameter with name, " + id + ", is unknown");
        }
        if (parameter.isFixed) {
            writeParameter(id, dimension, parameter.initial, Double.NaN, Double.NaN, writer);
        } else
        if (parameter.priorType == PriorType.UNIFORM_PRIOR || parameter.priorType == PriorType.TRUNC_NORMAL_PRIOR) {
            writeParameter(id, dimension, parameter.initial, parameter.uniformLower, parameter.uniformUpper, writer);
        } else {
            writeParameter(id, dimension, parameter.initial, parameter.lower, parameter.upper, writer);
        }
    }


    /**
     * write a parameter
     *
     * @param id        the id
     * @param dimension the dimension
     * @param value     the value
     * @param lower     the lower bound
     * @param upper     the upper bound
     * @param writer    the writer
     */
    public void writeParameter(String id, int dimension, double value, double lower, double upper, XMLWriter writer) {
        ArrayList<Attribute.Default> attributes = new ArrayList<Attribute.Default>();
        attributes.add(new Attribute.Default<String>("id", id));
        if (dimension > 1) {
            attributes.add(new Attribute.Default<String>("dimension", dimension + ""));
        }
        if (!Double.isNaN(value)) {
            attributes.add(new Attribute.Default<String>("value", multiDimensionValue(dimension, value)));
        }
        if (!Double.isNaN(lower)) {
            attributes.add(new Attribute.Default<String>("lower", multiDimensionValue(dimension, lower)));
        }
        if (!Double.isNaN(upper)) {
            attributes.add(new Attribute.Default<String>("upper", multiDimensionValue(dimension, upper)));
        }

        Attribute[] attrArray = new Attribute[attributes.size()];
        for (int i = 0; i < attrArray.length; i++) {
            attrArray[i] = attributes.get(i);
        }

        writer.writeTag(ParameterParser.PARAMETER, attrArray, true);
    }

    void writeSumStatisticColumn(XMLWriter writer, String name, String label) {
        writer.writeOpenTag(Columns.COLUMN,
                new Attribute[]{
                        new Attribute.Default<String>(Columns.LABEL, label),
                        new Attribute.Default<String>(Columns.DECIMAL_PLACES, "0"),
                        new Attribute.Default<String>(Columns.WIDTH, "12")
                }
        );
        writer.writeTag("sumStatistic", new Attribute.Default<String>("idref", name), true);
        writer.writeCloseTag(Columns.COLUMN);
    }

    private String multiDimensionValue(int dimension, double value) {
        String multi = "";

        multi += value + "";
        for (int i = 2; i <= dimension; i++)
            multi += " " + value;

        return multi;
    }

}
