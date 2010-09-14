/*
 * TaskListener.java
 *
 * Copyright (c) 2009 JAM Development Team
 *
 * This package is distributed under the Lesser Gnu Public Licence (LGPL)
 *
 */

package dr.app.gui.util;

public interface TaskListener {

    void taskFinished();

    void taskCanceled();
}
