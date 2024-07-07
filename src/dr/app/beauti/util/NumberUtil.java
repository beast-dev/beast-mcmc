/*
 * NumberUtil.java
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

package dr.app.beauti.util;

import java.text.DecimalFormat;

/**
 * @author Walter Xie
 */
public class NumberUtil {

    protected static DecimalFormat formatter = new DecimalFormat("0.####E0");
    protected static DecimalFormat formatter2 = new DecimalFormat("####0.####");

    public static String formatDecimal(double value, int maxFractionDigits1, int maxFractionDigits2) {
        formatter.setMaximumFractionDigits(maxFractionDigits1);
        formatter2.setMaximumFractionDigits(maxFractionDigits2);

        if (value > 0 && (Math.abs(value) < 0.001 || Math.abs(value) >= 100000.0)) {
            return formatter.format(value);
        } else {
            return formatter2.format(value);
        }
    }
}