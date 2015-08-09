/*
 * DataTypeUtils.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
 *
 * This file is part of BEAST.
 * See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * BEAST is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 *  BEAST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

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
