package dr.app.beauti.types;


/**
 * @author Alexei Drummond
 * @author Walter Xie
 */
public enum PriorScaleType {

	NONE,
	TIME_SCALE,
    GROWTH_RATE_SCALE,
    BIRTH_RATE_SCALE,
    SUBSTITUTION_RATE_SCALE,
    LOG_STDEV_SCALE,
    SUBSTITUTION_PARAMETER_SCALE,
    T50_SCALE,
    ROOT_RATE_SCALE,
    LOG_VAR_SCALE,
    ORIGIN_SCALE,
    GAMMA_STAR_BEAST;


//    PriorScaleType(String name) {
//        this.name = name;
//    }
//
//    public String toString() {
//        return name;
//    }
//
//    private final String name;

}
