/*
 * UniformIntegerOperatorTest.java
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

package test.dr.evomodel.operators;

import dr.inference.model.Bounds;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.inference.operators.UniformIntegerOperator;
import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

/**
 * @author Walter Xie
 */
public class UniformIntegerOperatorTest extends TestCase {
    private final int dimension = 3;
    private int[][] count;

    public static Test suite() {
        return new TestSuite(UniformIntegerOperatorTest.class);
    }

    public void testParameterBound() {
        count = new int[dimension][4]; // 4 vaules {0, 1, 2, 3}
        Parameter parameter = new Parameter.Default(new double[]{1.0, 0.0, 3.0});
        parameter.addBounds(new Parameter.DefaultBounds(new double[]{3.0, 3.0, 3.0}, new double[dimension]));
        UniformIntegerOperator uniformIntegerOperator = new UniformIntegerOperator(parameter, 0, 3, 10, 1);
        for (int i = 0; i < 300; i++) {
            uniformIntegerOperator.doOperation();
            countParaValueFrequency(uniformIntegerOperator.getVariable());
        }
        printCount("Parameter (Double) lower = 0, upper = 3");

        assertTrue("Expected count[0][0-3] > 0", count[0][0] > 0 && count[0][1] > 0 && count[0][2] > 0 && count[0][3] > 0);
        assertTrue("Expected count[1][0-3] > 0", count[1][0] > 0 && count[1][1] > 0 && count[1][2] > 0 && count[1][3] > 0);
        assertTrue("Expected count[2][0-3] > 0", count[2][0] > 0 && count[2][1] > 0 && count[2][2] > 0 && count[2][3] > 0);
    }

    public void testIntegerParameterStaircaseBound() {
        count = new int[dimension][3]; // 3 vaules
        Variable<Integer> parameterInt = new Variable.I(new int[dimension-1]); // dimension = 3
        parameterInt.addBounds(new Bounds.Staircase(parameterInt)); // integer index parameter size = real size - 1
        UniformIntegerOperator uniformIntegerOperator = new UniformIntegerOperator(parameterInt, 10, 1);
        for (int i = 0; i < 300; i++) {
            uniformIntegerOperator.doOperation();
            countParaValueFrequency(uniformIntegerOperator.getVariable());
        }
        printCount("Integer Parameter using Staircase Bound");

//        assertTrue("Expected count[0][0] > 0", count[0][0] > 0);
        assertTrue("Expected count[1][0] && [1][1] > 0", count[1][0] > 0 && count[1][1] > 0);
        assertTrue("Expected count[2][0] && [2][1] && [2][2] > 0", count[2][0] > 0 && count[2][1] > 0 && count[2][2] > 0);
    }

    private void countParaValueFrequency(Variable para) {
        for (int i = 0; i < para.getSize(); i++) {
            int j;
            if (para.getValue(i) instanceof Double) {
                j = (int) (double) (Double) para.getValue(i);
                count[i][j] += 1;
            } else {
                j = (Integer) para.getValue(i);
                count[i+1][j] += 1; // integer index parameter size = real size - 1
            }
        }
    }

    private void printCount(String m) {
        System.out.println("\n-----------------------\n");
        System.out.println(m);
        int i = 0;
        for (int[] row : count) {
            System.out.print(i + " : ");
            i++;
            for (int col : row) {
                System.out.printf("\t%10d", col);
            }
            System.out.print("\n");
        }
    }
}
