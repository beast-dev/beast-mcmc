/*
 * NotBooleanStatistic.java
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

package dr.inference.model;

import java.util.List;

/**
 * @author Marc Suchard
 */
public class NotBooleanStatistic extends Statistic.Abstract implements BooleanStatistic {

    public NotBooleanStatistic(BooleanStatistic originalStatistic) {
        this(originalStatistic, null);
    }

    public NotBooleanStatistic(BooleanStatistic originalStatistic, List<Integer> mark) {
        super(originalStatistic.getStatisticName());
        this.originalStatistic = originalStatistic;
        this.mark = mark;
    }

    /**
     * @return boolean result of test.
     */
    public double getStatisticValue(int dim) {
        return getBoolean(dim) ? 1.0 : 0.0;
    }

    @Override
    public boolean getBoolean(int dim) {
        boolean rtnValue = originalStatistic.getBoolean(dim);
        if (mark == null || mark.contains(dim)) {
            rtnValue = !rtnValue;
        }
        return rtnValue;
    }

    public int getDimension() {
        return originalStatistic.getDimension();
    }

    private final BooleanStatistic originalStatistic;
    private final List<Integer> mark;
}
