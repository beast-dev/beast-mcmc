package dr.evolution.datatype;

/**
 * @author Andrew Rambaut
 * <p>
 * Tree data type. This is a placeholder data type for when a tree is the data itself.
 */

public class TreeDataType extends DataType {

    public static final TreeDataType INSTANCE = new TreeDataType();

    public TreeDataType() {

    }

    @Override
    public char[] getValidChars() {
        return new char[0];
    }

    @Override
    public String getDescription() {
        return "Tree";
    }

    @Override
    public int getType() {
        return TREE;
    }
}
