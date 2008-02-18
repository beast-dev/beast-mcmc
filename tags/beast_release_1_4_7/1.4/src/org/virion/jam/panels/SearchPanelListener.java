/*
 * SearchPanelListener.java
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
 * An interface for listeners to the SearchPanel class.
 * @author Andrew Rambaut
 * Date: Jul 26, 2004
 * Time: 5:37:15 PM
 */
public interface SearchPanelListener {

    /**
     * Called when the user requests a search by pressing return having
     * typed a search string into the text field. If the continuousUpdate
     * flag is true then this method is called when the user types into
     * the text field.
     * @param searchString the user's search string
     */
    void searchStarted(String searchString);

    /**
     * Called when the user presses the cancel search button or presses
     * escape while the search is in focus.
     */
    void searchStopped();

}
