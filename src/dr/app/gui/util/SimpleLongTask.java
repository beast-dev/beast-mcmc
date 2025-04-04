/*
 * SimpleLongTask.java
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

package dr.app.gui.util;

public abstract class SimpleLongTask extends LongTask {

    boolean background = false;
    private SwingWorker worker = null;
    public int current = 0;
    public int length = 1;
    public boolean pleaseStop = false;
    public String message = "";
    public String description = "";

    /**
     * Called to find out how much work needs
     * to be done.
     */
    public int getLengthOfTask() {
        return length;
    }

    /**
     * Called to find out how much has been done.
     */
    public int getCurrent() {
        return current;
    }

    /**
     * Called to stop task.
     */
    public void stop() {
        pleaseStop = true;
    }

    public String getMessage() {
        return message;
    }

    public String getDescription() {
        return description;
    }
}
