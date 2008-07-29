package dr.app.beauti.options;

import dr.app.beauti.PriorType;
import dr.evolution.util.TaxonList;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 */
public class Parameter {

    /**
     * A constructor for "special" parameters which are not user-configurable
     *
     * @param name        the name
     * @param description the description
     */
    public Parameter(String name, String description) {
        this.baseName = name;
        this.description = description;
        this.scale = BeautiOptions.NONE;
        this.isNodeHeight = false;
        this.isStatistic = false;
        this.taxa = null;
        this.priorType = PriorType.NONE;
        this.initial = Double.NaN;
        this.lower = Double.NaN;
        this.upper = Double.NaN;
    }

    public Parameter(String name, String description, int scale,
                     double initial, double lower, double upper) {
        this.baseName = name;
        this.description = description;
        this.initial = initial;
        this.isNodeHeight = false;
        this.isStatistic = false;

        this.taxa = null;

        this.priorType = PriorType.UNIFORM_PRIOR;
        this.scale = scale;
        this.priorEdited = false;
        this.lower = lower;
        this.upper = upper;

        uniformLower = lower;
        uniformUpper = upper;
    }

    public Parameter(TaxonList taxa, String description) {
        this.taxa = taxa;
        this.baseName = null;
        this.description = description;

        this.isNodeHeight = true;
        this.isStatistic = true;
        this.priorType = PriorType.NONE;
        this.scale = BeautiOptions.TIME_SCALE;
        this.priorEdited = false;
        this.lower = 0.0;
        this.upper = Double.MAX_VALUE;

        uniformLower = lower;
        uniformUpper = upper;
    }

    public Parameter(String name, String description, boolean isDiscrete) {
        this.taxa = null;

        this.baseName = name;
        this.description = description;

        this.isNodeHeight = false;
        this.isStatistic = true;
        this.isDiscrete = isDiscrete;
        this.priorType = PriorType.UNIFORM_PRIOR;
        this.scale = BeautiOptions.NONE;
        this.priorEdited = false;
        this.initial = Double.NaN;
        this.lower = Double.NaN;
        this.upper = Double.NaN;
    }

    public Parameter(String name, String description, double lower, double upper) {
        this.taxa = null;

        this.baseName = name;
        this.description = description;

        this.isNodeHeight = false;
        this.isStatistic = true;
        this.isDiscrete = false;
        this.priorType = PriorType.UNIFORM_PRIOR;
        this.scale = BeautiOptions.NONE;
        this.priorEdited = false;
        this.initial = Double.NaN;
        this.lower = lower;
        this.upper = upper;

        uniformLower = lower;
        uniformUpper = upper;
    }

    public Parameter(String name, String description, boolean isNodeHeight,
                     double initial, double lower, double upper) {
        this.baseName = name;
        this.description = description;
        this.initial = initial;

        this.taxa = null;

        this.isNodeHeight = isNodeHeight;
        this.isStatistic = false;
        this.priorType = PriorType.NONE;
        this.scale = BeautiOptions.TIME_SCALE;
        this.priorEdited = false;
        this.lower = lower;
        this.upper = upper;

        uniformLower = lower;
        uniformUpper = upper;
    }

    public void setPrefix(String prefix) {
        this.prefix = prefix;
    }

    public String getFullName() {
        if (prefix != null) return prefix + baseName;
        return baseName;
    }

    public String getName() {
        if (taxa != null) {
            return "tmrca(" + taxa.getId() + ")";
        } else {
            return getFullName();
        }
    }

    public String getXMLName() {
        if (taxa != null) {
            return "tmrca_" + taxa.getId();
        } else {
            return getFullName();
        }
    }

    public String getDescription() {
        if (taxa != null) {
            return description + taxa.getId();
        } else if (prefix != null) {
            return description + " of partition " + prefix;
        }
        return description;
    }

    private final String baseName;

    private String prefix = null;

    private final String description;
    public double initial;

    public final TaxonList taxa;

    public boolean isDiscrete = false;

    public boolean isFixed = false;
    public final boolean isNodeHeight;
    public final boolean isStatistic;

    public PriorType priorType;
    public boolean priorEdited;
    public final int scale;
    public double lower;
    public double upper;

    public double uniformUpper = 0.0;
    public double uniformLower = 0.0;
    public double exponentialMean = 1.0;
    public double exponentialOffset = 0.0;
    public double normalMean = 1.0;
    public double normalStdev = 1.0;
    public double logNormalMean = 0.0;
    public double logNormalStdev = 1.0;
    public double logNormalOffset = 0.0;
    public double gammaAlpha = 1.0;
    public double gammaBeta = 1.0;
    public double gammaOffset = 0.0;
    public double poissonMean = 1.0;
    public double poissonOffset = 0.0;
}
