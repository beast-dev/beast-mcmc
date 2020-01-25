package dr.evolution.datatype;

/**
 * @author Gabriel Hassler
 * <p>
 * Dummy data type. This is a place holder to allow a data type with no data.
 * None of the methods will return anything useful.
 */

public class DummyDataType extends DataType {


    public DummyDataType() {

    }

    @Override
    public char[] getValidChars() {
        return new char[0];
    }

    @Override
    public String getDescription() {
        return "NA";
    }

    @Override
    public int getType() {
        return 9;
    }
}
