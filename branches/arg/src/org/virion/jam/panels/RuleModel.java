/*
 * RuleModel.java
 *
 * Copyright (C) 2002-2006 Alexei Drummond and Andrew Rambaut
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

package org.virion.jam.panels;


/**
 * RuleModel.
 *
 * @author Andrew Rambaut
 * @version $Id: RuleModel.java,v 1.2 2006/09/09 18:16:16 rambaut Exp $
 */


public interface RuleModel {

    /**
     * Returns an array of strings to be presented as a combo box which
     * are available fields to define rules on.
     *
     * @return the field names
     */
    Object[] getFields();

    /**
     * Returns an array of strings to be presented as a combo box which
     * are possible rule conditions for the specified field.
     *
     * @return the condition names
     */
    Object[] getConditions(Object field);

    /**
     * Returns an array of strings to be presented as a combo box which
     * are possible values for the field. Should return null if a text
     * box is required.
     *
     * @return the values
     */
    Object[] getValues(Object field, Object condition);

}