package dr.evolution.datatype;

/**
 * @author JT 
 * Borrowed from ContinuousDataType
 *
 * Continuous data type. This is a place holder to allow mixing of continuous with
 * discrete traits. None of the methods will return anything useful.
 */
public class IntegerDataType extends DataType {

    public static final String DESCRIPTION = "integer";
    public static final IntegerDataType INSTANCE = new IntegerDataType();

    /**
     * Constructor
     */
    public IntegerDataType(){
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
