/*
 * SimpleLongTask.java
 *
 * Copyright (c) 2009 JAM Development Team
 *
 * This package is distributed under the Lesser Gnu Public Licence (LGPL)
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
