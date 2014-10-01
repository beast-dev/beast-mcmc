package dr.evoxml.util;

import dr.evolution.datatype.DataType;
import dr.evolution.datatype.GeneticCode;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;

/**
 * @author Alexei Drummond
 *
 * @version $Id$
 */
public class DataTypeUtils {

    public static DataType getDataType(XMLObject xo) throws XMLParseException {

        DataType dataType = null;

        if (xo.hasAttribute(DataType.DATA_TYPE)) {
            String dataTypeStr = xo.getStringAttribute(DataType.DATA_TYPE);

            if (xo.hasAttribute(GeneticCode.GENETIC_CODE)) {
                dataTypeStr += "-" + xo.getStringAttribute(GeneticCode.GENETIC_CODE);
            }

            dataType = DataType.getRegisteredDataTypeByName(dataTypeStr);
        }

        for (int i = 0; i < xo.getChildCount(); i++) {

            Object child = xo.getChild(i);
            if (child instanceof DataType) {
                if (dataType != null) {
                    throw new XMLParseException("Multiple dataTypes defined for alignment element");
                }

                dataType = (DataType) child;
            }
        }

        return dataType;
    }

}
