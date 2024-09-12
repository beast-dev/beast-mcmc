/*
 * BeautiModelIDProvider.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
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
 *
 */

package dr.app.beauti.components;

import dr.app.beauti.util.XMLWriter;

public interface BeautiModelIDProvider {

    /**
     * Returns the parser xml tag for writing.
     *
     * @return the ID
     */
    String getParserTag();

    /**
     * Returns the default ID for an xml component.
     *
     * @param modelName the model name
     * @return the ID
     */
    String getId(String modelName);

    /**
     * Write the ID ref of the ID provider from the model name.
     *
     * @param writer the XML writer
     * @param modelName the model name
     */
    default void writeIDrefFromName(XMLWriter writer, String modelName) {
        writer.writeIDref(this.getParserTag(), this.getId(modelName));
    }

    /**
     * Write the ID ref of the ID provider directly from the ID.
     *
     * @param writer the XML writer
     * @param modelID the model name
     */
    default void writeIDrefFromID(XMLWriter writer, String modelID) {
        writer.writeIDref(this.getParserTag(), modelID);
    }

    default BeautiParameterIDProvider getBeautiParameterIDProvider(String parameterKey) {
        throw new IllegalArgumentException("No parameter to access.");
    }

}
