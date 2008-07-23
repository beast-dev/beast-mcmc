package dr.app.beauti.options;

/**
 * @author Alexei Drummond
 * @author Andrew Rambaut
 */
public class Operator {
    public Operator(String name, String description, Parameter parameter, String operatorType, double tuning, double weight) {
        this.name = name;
        this.description = description;
        this.parameter1 = parameter;
        this.parameter2 = null;

        this.type = operatorType;
        this.tuningEdited = false;
        this.tuning = tuning;
        this.weight = weight;

        this.inUse = true;
    }

    public Operator(String name, String description,
                    Parameter parameter1, Parameter parameter2,
                    String operatorType, double tuning, double weight) {
        this.name = name;
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
            if (type.equals(BeautiOptions.SCALE)) {
                prefix = "Scales the ";
            } else if (type.equals(BeautiOptions.RANDOM_WALK)) {
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

    public final String name;
    public final String description;

    public final String type;
    public boolean tuningEdited;
    public double tuning;
    public double weight;
    public boolean inUse;

    public final Parameter parameter1;
    public final Parameter parameter2;

}
