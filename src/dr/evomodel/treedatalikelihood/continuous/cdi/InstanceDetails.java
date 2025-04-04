/*
 * InstanceDetails.java
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

package dr.evomodel.treedatalikelihood.continuous.cdi;

import java.io.Serializable;

/**
 * @author Marc A. Suchard
 */
public class InstanceDetails implements Serializable {
    private int resourceNumber;
    private long flags;
    private String resourceName;
    private String implementationName;

    public InstanceDetails() {
    }

    public int getResourceNumber() {
        return resourceNumber;
    }

    public void setResourceNumber(int i) {
        resourceNumber = i;
    }

    public String getResourceName() {
        return resourceName;
    }

    public void setResourceName(String str) {
        resourceName = str;
    }

    public String getImplementationName() {
        return implementationName;
    }

    public void setImplementationName(String str) {
        implementationName = str;
    }

    public long getFlags() {
        return flags;
    }

    public void setFlags(long flags) {
        this.flags = flags;
    }

    public String toString() {
        return CDIFlag.toString(this.getFlags());
    }
}
