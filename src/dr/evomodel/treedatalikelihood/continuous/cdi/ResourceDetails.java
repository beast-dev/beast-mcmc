/*
 * ResourceDetail.java
 *
 * Copyright (c) 2002-2016 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodel.treedatalikelihood.continuous.cdi;

public class ResourceDetails {
    private final int number;
    private String name;
    private String description;
    private long flags;

    public ResourceDetails(int var1) {
        this.number = var1;
    }

    public int getNumber() {
        return this.number;
    }

    public String getName() {
        return this.name;
    }

    public void setName(String var1) {
        this.name = var1;
    }

    public String getDescription() {
        return this.description;
    }

    public void setDescription(String var1) {
        this.description = var1;
    }

    public long getFlags() {
        return this.flags;
    }

    public void setFlags(long var1) {
        this.flags = var1;
    }

    public String toString() {
        StringBuilder var1 = new StringBuilder();
        var1.append("").append(this.getNumber()).append(" : ").append(this.getName()).append("\n");
        if(this.getDescription() != null) {
            String[] var2 = this.getDescription().split("\\|");
            int var4 = var2.length;

            for(int var5 = 0; var5 < var4; ++var5) {
                String var6 = var2[var5];
                if(var6.trim().length() > 0) {
                    var1.append("    ").append(var6.trim()).append("\n");
                }
            }
        }

        var1.append("    Flags:");
        var1.append(CDIFlag.toString(this.getFlags()));
        var1.append("\n");
        return var1.toString();
    }
}
