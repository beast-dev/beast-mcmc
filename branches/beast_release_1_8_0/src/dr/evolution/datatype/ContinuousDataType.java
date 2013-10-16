package dr.evolution.datatype;

/**
 * @author Andrew Rambaut
 *
 * Continuous data type. This is a place holder to allow mixing of continuous with
 * discrete traits. None of the methods will return anything useful.
 */
public class ContinuousDataType extends DataType {

    public static final String DESCRIPTION = "continuous";
    public static final ContinuousDataType INSTANCE = new ContinuousDataType();

    /**
     * Constructor
     */
    public ContinuousDataType(){
        stateCount = 0;
        ambiguousStateCount = 0;
    }

    @Override
    public char[] getValidChars() {
        return null;
    }

    /**
     * @return the description of the data type
     */
    public String getDescription() {
		return DESCRIPTION;
	}

    public int getType(){
        return CONTINUOUS;
    }

}
