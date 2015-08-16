/*
 * LongTaskMonitor.java
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

import javax.swing.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

public class LongTaskMonitor {

    public final static int ONE_SECOND = 1000;

    private ProgressMonitor progressMonitor;
    private Timer timer;
    private LongTask task;
    private JComponent parent;
    private TaskListener listener;

    public LongTaskMonitor(LongTask task, JComponent parent, TaskListener listener) {

        this.task = task;
        this.parent = parent;
        this.listener = listener;

        //Create a timer.
        timer = new Timer(ONE_SECOND, new TimerListener());

        progressMonitor = new ProgressMonitor(parent,
                task.getDescription(),
                "", 0, task.getLengthOfTask());
        progressMonitor.setProgress(0);
        progressMonitor.setMillisToPopup(ONE_SECOND);

        task.go();
        timer.start();
    }

    /**
     * The actionPerformed method in this class
     * is called each time the Timer "goes off".
     */
    class TimerListener implements ActionListener {
        public void actionPerformed(ActionEvent evt) {
            if (progressMonitor.isCanceled() || task.done()) {

                boolean canceled = progressMonitor.isCanceled();
                progressMonitor.close();
                task.stop();
                timer.stop();
                if (canceled) {
                    listener.taskCanceled();
                } else {
                    listener.taskFinished();
                    //Toolkit.getDefaultToolkit().beep();
                }
            } else {
                progressMonitor.setNote(task.getMessage());
                progressMonitor.setProgress(task.getCurrent());
            }
        }
    }
}
