/*
 * LogColumn.java
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

package dr.inference.loggers;

import java.io.Serializable;

/**
 * An interface for a column in a log.
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: LogColumn.java,v 1.5 2005/05/24 20:25:59 rambaut Exp $
 */

public interface LogColumn extends Serializable {

    /**
     * Set the label (heading) for this column
     *
     * @param label the column label
     */
    void setLabel(String label);

    /**
     * @return the label (heading) for this column
     */
    String getLabel();

    /**
     * Set the minimum width in characters for this column
     *
     * @param minimumWidth the minimum width in characters
     */
    void setMinimumWidth(int minimumWidth);

    /**
     * @return the minimum width in characters for this column
     */
    int getMinimumWidth();

    /**
     * Returns a string containing the current value for this column with
     * appropriate formatting.
     *
     * @return the formatted string.
     */
    String getFormatted();

    public abstract class Abstract implements LogColumn {

        private String label;
        private int minimumWidth;

        public Abstract(String label) {

            setLabel(label);
            minimumWidth = -1;
        }

        public void setLabel(String label) {
            if (label == null) throw new IllegalArgumentException("column label is null");
            this.label = label;
        }

        public String getLabel() {
            StringBuffer buffer = new StringBuffer(label);

            if (minimumWidth > 0) {
                while (buffer.length() < minimumWidth) {
                    buffer.append(' ');
                }
            }

            return buffer.toString();
        }

        public void setMinimumWidth(int minimumWidth) {
            this.minimumWidth = minimumWidth;
        }

        public int getMinimumWidth() {
            return minimumWidth;
        }

        public final String getFormatted() {
            StringBuffer buffer = new StringBuffer(getFormattedValue());

            if (minimumWidth > 0) {
                while (buffer.length() < minimumWidth) {
                    buffer.append(' ');
                }
            }

            return buffer.toString();
        }

        protected abstract String getFormattedValue();
    }

    public class Default extends Abstract {

        private Object object;

        public Default(String label, Object object) {
            super(label);
            this.object = object;
        }

        protected String getFormattedValue() {
            return object.toString();
        }
    }

}
