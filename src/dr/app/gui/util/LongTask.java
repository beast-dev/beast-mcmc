/*
 * LongTask.java
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

package dr.app.gui.util;

public abstract class LongTask {

    private SwingWorker worker = null;
    Object answer;
    boolean finished = false;

    /**
     * Does the actual work and returns some kind of result.
     */
    public abstract Object doWork() throws java.lang.Exception;

    /**
     * Called to start the task.
     */
    public final void go() {
        worker = new SwingWorker() {
            public Object construct() {

                try {
                    answer = doWork();
                } catch (java.lang.Exception e) {
                    throw new RuntimeException(e.toString());
                }
                finished = true;
                return answer;
            }
        };
        worker.start();
    }

    public final Object getAnswer() {
        return answer;
    }

    /**
     * Called to find out how much work needs
     * to be done.
     */
    public abstract int getLengthOfTask();

    /**
     * Called to find out how much has been done.
     */
    public abstract int getCurrent();

    /**
     * Called to stop task.
     */
    public void stop() {
        finished = true;
    }

    /**
     * Called to find out if the task has completed.
     */
    public boolean done() {
        return finished;
    }

    /**
     * Called to get the current message of the task.
     */
    public abstract String getMessage();

    /**
     * Called to get the description of this task.
     */
    public String getDescription() {
        return "Running a long task...";
    }
}
