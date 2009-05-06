package dr.app.beauti.generator;

import dr.app.beauti.XMLWriter;
import dr.app.beauti.components.ComponentFactory;
import dr.app.beauti.options.BeautiOptions;
import dr.app.beauti.options.ModelOptions;
import dr.app.beauti.priorsPanel.PriorType;
import dr.inference.loggers.Columns;
import dr.inference.model.ParameterParser;
import dr.util.Attribute;
import dr.xml.XMLParser;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * @author Alexei Drummond
 */
public abstract class Generator {

    protected final BeautiOptions options;


    protected Generator(BeautiOptions options) {
        this.options = options;
    }

    public Generator(BeautiOptions options, ComponentFactory[] components) {
        this.options = options;
        if (components != null) {
            for (ComponentFactory component : components) {
                this.components.add(component.getGenerator(options));
            }
        }
    }

    /**
     * fix a parameter
     *
     * @param id    the id
     * @param value the value
     */
    public void fixParameter(String id, double value) {
        dr.app.beauti.options.Parameter parameter = options.getParameter(id);
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
    public void writeParameterRef(String wrapperName, String id, XMLWriter writer) {
        writer.writeOpenTag(wrapperName);
        writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default<String>(XMLParser.IDREF, id), true);
        writer.writeCloseTag(wrapperName);
    }

    /**
     * write a parameter
     *
     * @param id     the id
     * @param writer the writer
     */
    public void writeParameter(String id, ModelOptions options, XMLWriter writer) {
        dr.app.beauti.options.Parameter parameter = options.getParameter(id);
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
     * @param id     the id
     * @param writer the writer
     */
    public void writeParameter(String wrapperName, String id, int dimension, XMLWriter writer) {
        writer.writeOpenTag(wrapperName);
        writeParameter(id, dimension, writer);
        writer.writeCloseTag(wrapperName);
    }

    /**
     * write a parameter
     *
     * @param id     the id
     * @param writer the writer
     */
    public void writeParameter(String wrapperName, String id, ModelOptions options, XMLWriter writer) {
        writer.writeOpenTag(wrapperName);
        writeParameter(id, options, writer);
        writer.writeCloseTag(wrapperName);
    }

    /**
     * write a parameter
     *
     * @param id        the id
     * @param dimension the dimension
     * @param writer    the writer
     */
    public void writeParameter(String id, int dimension, XMLWriter writer) {
        dr.app.beauti.options.Parameter parameter = options.getParameter(id);
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
        attributes.add(new Attribute.Default<String>(XMLParser.ID, id));
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
    public void writeParameter(String wrapperName, String id, int dimension, double value, double lower, double upper, XMLWriter writer) {
        writer.writeOpenTag(wrapperName);
        writeParameter(id, dimension, value, lower, upper, writer);
        writer.writeCloseTag(wrapperName);
    }


    void writeSumStatisticColumn(XMLWriter writer, String name, String label) {
        writer.writeOpenTag(Columns.COLUMN,
                new Attribute[]{
                        new Attribute.Default<String>(Columns.LABEL, label),
                        new Attribute.Default<String>(Columns.DECIMAL_PLACES, "0"),
                        new Attribute.Default<String>(Columns.WIDTH, "12")
                }
        );
        writer.writeTag("sumStatistic", new Attribute.Default<String>(XMLParser.IDREF, name), true);
        writer.writeCloseTag(Columns.COLUMN);
    }

    private String multiDimensionValue(int dimension, double value) {
        String multi = "";

        multi += value + "";
        for (int i = 2; i <= dimension; i++)
            multi += " " + value;

        return multi;
    }

    public void writeReferenceComment(String[] lines, XMLWriter writer) {
        for (String line : lines)
            writer.writeComment(line);
    }

    protected void generateInsertionPoint(final ComponentGenerator.InsertionPoint ip, final XMLWriter writer) {
        generateInsertionPoint(ip, null, writer);
    }

    protected void generateInsertionPoint(final ComponentGenerator.InsertionPoint ip, final Object item, final XMLWriter writer) {
        for (ComponentGenerator component : components) {
            if (component.usesInsertionPoint(ip)) {
                component.generateAtInsertionPoint(ip, item, writer);
            }
        }
    }

    private final List<ComponentGenerator> components = new ArrayList<ComponentGenerator>();
}
