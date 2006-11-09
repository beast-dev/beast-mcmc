package dr.evoxml;

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

            if (dataType == null) dataType = (DataType)xo.getChild(DataType.class);

        return dataType;
    }

}
