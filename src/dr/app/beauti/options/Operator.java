package dr.app.beauti.options;

/**
 * @author Alexei Drummond
 * @author Andrew Rambaut
 */
public class Operator {

    public Operator(String name,
                    String description,
                    Parameter parameter,
                    OperatorType operatorType,
                    double tuning,
                    double weight,
                    ModelOptions options) {

        this.baseName = name;
        this.description = description;
        this.parameter1 = parameter;
        this.parameter2 = null;

        this.type = operatorType;
        this.tuningEdited = false;
        this.tuning = tuning;
        this.weight = weight;

        this.options = options;

        this.inUse = true;
    }

    public Operator(String name, String description,
                    Parameter parameter1, Parameter parameter2,
                    OperatorType operatorType, double tuning, double weight, ModelOptions options) {
        this.baseName = name;
        this.description = description;
        this.parameter1 = parameter1;
        this.parameter2 = parameter2;

        this.type = operatorType;
        this.tuningEdited = false;
        this.tuning = tuning;
        this.weight = weight;

        this.inUse = true;
    }

    public String getDescription() {
        if (description == null || description.length() == 0) {
            String prefix = "";
            if (type == OperatorType.SCALE) {
                prefix = "Scales the ";
            } else if (type == OperatorType.RANDOM_WALK) {
                prefix = "A random-walk on the ";
            }
            return prefix + parameter1.getDescription();
        } else {
            return description;
        }
    }

    public boolean isTunable() {
        return tuning > 0;
    }

    public void setPrefix(String prefix) {
        this.prefix = prefix;
    }

    public String getName() {
        String name = baseName;
        if (prefix != null) name = prefix + "." + baseName;
        return name;
    }

    public ModelOptions getModelOptions() {
        return options;
    }

    //public void addAttribute(String name, String value) {
    //    attributes.put(name, value);
    //}

    private final String baseName;
    public String prefix = null;

    public final String description;

    public final OperatorType type;
    public boolean tuningEdited;
    public double tuning;
    public double weight;
    public boolean inUse;

    public final Parameter parameter1;
    public final Parameter parameter2;

    private ModelOptions options;

    //private Map<String, String> attributes = new HashMap<String, String>();

}
