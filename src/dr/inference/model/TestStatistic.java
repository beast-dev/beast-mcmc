/*
 * TestStatistic.java
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

package dr.inference.model;

import dr.util.Attribute;

/**
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @version $Id: TestStatistic.java,v 1.12 2005/07/11 14:06:25 rambaut Exp $
 */
public class TestStatistic extends Statistic.Abstract implements BooleanStatistic {

    private Attribute attribute = null;
    private Attribute attribute2 = null;
    private double testValue1, testValue2;
    private final int mode;

    public static final int EQUALS = 0;
    public static final int GREATER_THAN = 1;
    public static final int LESS_THAN = 2;
    public static final int INSIDE = 3;
    public static final int OUTSIDE = 4;

    public TestStatistic(String name, Attribute attr, double value, int mode) {
        super(name);

        attribute = attr;
        testValue1 = value;
        this.mode = mode;
    }

    /**
     * A constructor for modes that require 2 values (i.e., INSIDE and OUTSIDE).
     */
    public TestStatistic(String name, Attribute attr, double value1, double value2, int mode) {
        super(name);

        attribute = attr;
        testValue1 = value1;
        testValue2 = value2;
        this.mode = mode;
    }

    /**
     * A constructor that compares 2 attributes (i.e., GREATER_THAN or LESS_THAN).
     */
    public TestStatistic(String name, Attribute attr1, Attribute attr2, int mode) {
        super(name);

        attribute = attr1;
        attribute2 = attr2;
        this.mode = mode;
    }

    public int getDimension() {
        return 1;
    }

    /**
     * @return boolean result of test.
     */
    public double getStatisticValue(int dim) {
        return getBoolean(dim) ? 1.0 : 0.0;
    }

    /**
     * @return boolean result of test.
     */
    public boolean getBoolean(int i) {

        double num;
        if (attribute instanceof Statistic) {
            num = ((Statistic) attribute).getStatisticValue(0);
        } else {
            num = ((Number) attribute.getAttributeValue()).doubleValue();
        }

        if (attribute2 != null) {
            if (attribute2 instanceof Statistic) {
                testValue1 = ((Statistic) attribute2).getStatisticValue(0);
            } else {
                testValue1 = ((Number) attribute2.getAttributeValue()).doubleValue();
            }
        }

        switch (mode) {
            case EQUALS:
                if (num == testValue1) return true;
                break;
            case GREATER_THAN:
                if (num > testValue1) return true;
                break;
            case LESS_THAN:
                if (num < testValue1) return true;
                break;
            case INSIDE:
                if (num > testValue1 && num < testValue2) return true;
                break;
            case OUTSIDE:
                if (num < testValue1 && num > testValue2) return true;
                break;
        }
        return false;
    }

}
