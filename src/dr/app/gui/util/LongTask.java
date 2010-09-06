/*
 * LongTask.java
 *
 * Copyright (c) 2009 JAM Development Team
 *
 * This package is distributed under the Lesser Gnu Public Licence (LGPL)
 *
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
