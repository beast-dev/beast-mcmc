/*
 * HeapSort.java
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

package dr.util;

import java.util.AbstractList;
import java.util.Comparator;
import java.util.Random;
import java.util.Vector;

/**
 * sorts numbers and comparable objects by treating contents of array as a binary tree.
 * KNOWN BUGS: There is a horrible amount of code duplication here!
 *
 * @author Alexei Drummond
 * @author Korbinian Strimmer
 */
public class HeapSort {

    //
    // Public stuff
    //

    /**
     * Sorts an array of indices to vector of comparable objects
     * into increasing order.
     */
    public static void sort(AbstractList array, int[] indices) {

        // ensures we are starting with valid indices
        for (int i = 0; i < indices.length; i++) {
            indices[i] = i;
        }

        int temp;
        int j, n = array.size();

        // turn input array into a heap
        for (j = n / 2; j > 0; j--) {
            adjust(array, indices, j, n);
        }

        // remove largest elements and put them at the end
        // of the unsorted region until you are finished
        for (j = n - 1; j > 0; j--) {
            temp = indices[0];
            indices[0] = indices[j];
            indices[j] = temp;
            adjust(array, indices, 1, j);
        }
    }

    /**
     * Sorts a vector of comparable objects into increasing order.
     */
    public static void sort(AbstractList array) {

        Object temp;
        int j, n = array.size();

        // turn input array into a heap
        for (j = n / 2; j > 0; j--) {
            adjust(array, j, n);
        }

        // remove largest elements and put them at the end
        // of the unsorted region until you are finished
        for (j = n - 1; j > 0; j--) {
            temp = array.get(0);
            array.set(0, array.get(j));
            array.set(j, temp);
            adjust(array, 1, j);
        }
    }

    /**
     * Sorts an array of comparable objects into increasing order.
     */
    public static void sort(Comparable[] array) {

        Comparable temp;
        int j, n = array.length;

        // turn input array into a heap
        for (j = n / 2; j > 0; j--) {
            adjust(array, j, n);
        }

        // remove largest elements and put them at the end
        // of the unsorted region until you are finished
        for (j = n - 1; j > 0; j--) {
            temp = array[0];
            array[0] = array[j];
            array[j] = temp;
            adjust(array, 1, j);
        }
    }

    /**
     * Sorts an array of objects into increasing order given a comparator.
     */
    public static void sort(Object[] array, Comparator c) {

        Object temp;
        int j, n = array.length;

        // turn input array into a heap
        for (j = n / 2; j > 0; j--) {
            adjust(array, c, j, n);
        }

        // remove largest elements and put them at the end
        // of the unsorted region until you are finished
        for (j = n - 1; j > 0; j--) {
            temp = array[0];
            array[0] = array[j];
            array[j] = temp;
            adjust(array, c, 1, j);
        }
    }

    /**
     * Sorts an array of doubles into increasing order.
     */
    public static void sort(double[] array) {

        double temp;
        int j, n = array.length;

        // turn input array into a heap
        for (j = n / 2; j > 0; j--) {
            adjust(array, j, n);
        }

        // remove largest elements and put them at the end
        // of the unsorted region until you are finished
        for (j = n - 1; j > 0; j--) {
            temp = array[0];
            array[0] = array[j];
            array[j] = temp;
            adjust(array, 1, j);
        }
    }

    /**
     * Sorts an array of doubles into increasing order, ingoring sign.
     */
    public static void sortAbs(double[] array) {

        double temp;
        int j, n = array.length;

        // turn input array into a heap
        for (j = n / 2; j > 0; j--) {
            adjustAbs(array, j, n);
        }

        // remove largest elements and put them at the end
        // of the unsorted region until you are finished
        for (j = n - 1; j > 0; j--) {
            temp = array[0];
            array[0] = array[j];
            array[j] = temp;
            adjustAbs(array, 1, j);
        }
    }

    /**
     * Sorts an array of indices into an array of doubles
     * into increasing order.
     */
    public static void sort(double[] array, int[] indices) {

        // ensures we are starting with valid indices
        for (int i = 0; i < indices.length; i++) {
            indices[i] = i;
        }

        int temp;
        int j, n = indices.length;

        // turn input array into a heap
        for (j = n / 2; j > 0; j--) {
            adjust(array, indices, j, n);
        }

        // remove largest elements and put them at the end
        // of the unsorted region until you are finished
        for (j = n - 1; j > 0; j--) {
            temp = indices[0];
            indices[0] = indices[j];
            indices[j] = temp;
            adjust(array, indices, 1, j);
        }
    }


    /**
     * test harness for heapsort algorithm
     */
    public static void main(String[] args) {

        int testSize = 100;

        // test array of Comparable objects

        ComparableDouble[] test = new ComparableDouble[testSize];

        Random random = new Random();

        for (int i = 0; i < test.length; i++) {
            test[i] = new ComparableDouble(random.nextInt(testSize * 10));
        }

        sort(test);
        for(ComparableDouble aTest : test) {
            System.out.print(aTest + " ");
        }
        System.out.println();

        // test index to Vector of Comparable objects

        Vector testv = new Vector();
        int[] indices = new int[testSize];

        for (int i = 0; i < testSize; i++) {
            testv.addElement(new ComparableDouble(random.nextInt(testSize * 10)));
        }

        sort(testv, indices);
        for (int i = 0; i < test.length; i++) {
            System.out.print(testv.elementAt(indices[i]) + " ");
        }
        System.out.println();


        // test index to array of doubles

        double[] testd = new double[testSize];
        //int[] indices = new int[testSize];

        for (int i = 0; i < testSize; i++) {
            testd[i] = random.nextInt(testSize * 10);
        }

        sort(testd, indices);
        for (int i = 0; i < test.length; i++) {
            System.out.print(testd[indices[i]] + " ");
        }
        System.out.println();

    }

    // PRIVATE STUFF

    /**
     * helps sort an array of indices into a vector of comparable objects.
     * Assumes that array[lower+1] through to array[upper] is
     * already in heap form and then puts array[lower] to
     * array[upper] in heap form.
     */
    private static void adjust(AbstractList<Comparable> array, int[] indices, int lower, int upper) {

        int j, k;
        int temp;

        j = lower;
        k = lower * 2;

        while (k <= upper) {
            if ((k < upper) && (array.get(indices[k - 1]).compareTo(array.get(indices[k])) < 0)) {
                k += 1;
            }
            if (array.get(indices[j - 1]).compareTo(array.get(indices[k - 1])) < 0) {
                temp = indices[j - 1];
                indices[j - 1] = indices[k - 1];
                indices[k - 1] = temp;
            }
            j = k;
            k *= 2;
        }
    }

    /**
     * helps sort an vector of comparable objects.
     * Assumes that array[lower+1] through to array[upper] is
     * already in heap form and then puts array[lower] to
     * array[upper] in heap form.
     */
    private static void adjust(AbstractList array, int lower, int upper) {

        int j, k;
        Object temp;

        j = lower;
        k = lower * 2;

        while (k <= upper) {
            if ((k < upper) && (((Comparable) array.get(k - 1)).compareTo(array.get(k)) < 0)) {
                k += 1;
            }
            if (((Comparable) array.get(j - 1)).compareTo(array.get(k - 1)) < 0) {
                temp = array.get(j - 1);
                array.set(j - 1, array.get(k - 1));
                array.set(k - 1, temp);
            }
            j = k;
            k *= 2;
        }
    }

    /**
     * Assumes that array[lower+1] through to array[upper] is
     * already in heap form and then puts array[lower] to
     * array[upper] in heap form.
     */
    private static void adjust(Comparable[] array, int lower, int upper) {

        int j, k;
        Comparable temp;

        j = lower;
        k = lower * 2;

        while (k <= upper) {
            if ((k < upper) && (array[k - 1].compareTo(array[k]) < 0)) {
                k += 1;
            }
            if (array[j - 1].compareTo(array[k - 1]) < 0) {
                temp = array[j - 1];
                array[j - 1] = array[k - 1];
                array[k - 1] = temp;
            }
            j = k;
            k *= 2;
        }
    }

    /**
     * Assumes that array[lower+1] through to array[upper] is
     * already in heap form and then puts array[lower] to
     * array[upper] in heap form.
     */
    private static void adjust(Object[] array, Comparator c, int lower, int upper) {

        int j, k;
        Object temp;

        j = lower;
        k = lower * 2;

        while (k <= upper) {
            if ((k < upper) && (c.compare(array[k - 1], array[k]) < 0)) {
                k += 1;
            }
            if (c.compare(array[j - 1], array[k - 1]) < 0) {
                temp = array[j - 1];
                array[j - 1] = array[k - 1];
                array[k - 1] = temp;
            }
            j = k;
            k *= 2;
        }
    }

    /**
     * helps sort an array of doubles.
     * Assumes that array[lower+1] through to array[upper] is
     * already in heap form and then puts array[lower] to
     * array[upper] in heap form.
     */
    private static void adjust(double[] array, int lower, int upper) {

        int j, k;
        double temp;

        j = lower;
        k = lower * 2;

        while (k <= upper) {
            if ((k < upper) && (array[k - 1] < array[k])) {
                k += 1;
            }
            if (array[j - 1] < array[k - 1]) {
                temp = array[j - 1];
                array[j - 1] = array[k - 1];
                array[k - 1] = temp;
            }
            j = k;
            k *= 2;
        }
    }

    /**
     * helps sort an array of doubles.
     * Assumes that array[lower+1] through to array[upper] is
     * already in heap form and then puts array[lower] to
     * array[upper] in heap form.
     */
    private static void adjustAbs(double[] array, int lower, int upper) {

        int j, k;
        double temp;

        j = lower;
        k = lower * 2;

        while (k <= upper) {
            if ((k < upper) && (Math.abs(array[k - 1]) < Math.abs(array[k]))) {
                k += 1;
            }
            if (Math.abs(array[j - 1]) < Math.abs(array[k - 1])) {
                temp = array[j - 1];
                array[j - 1] = array[k - 1];
                array[k - 1] = temp;
            }
            j = k;
            k *= 2;
        }
    }

    /**
     * helps sort an array of indices into an array of doubles.
     * Assumes that array[lower+1] through to array[upper] is
     * already in heap form and then puts array[lower] to
     * array[upper] in heap form.
     */
    private static void adjust(double[] array, int[] indices, int lower, int upper) {

        int j, k;
        int temp;

        j = lower;
        k = lower * 2;

        while (k <= upper) {
            if ((k < upper) && (array[indices[k - 1]] < array[indices[k]])) {
                k += 1;
            }
            if (array[indices[j - 1]] < array[indices[k - 1]]) {
                temp = indices[j - 1];
                indices[j - 1] = indices[k - 1];
                indices[k - 1] = temp;
            }
            j = k;
            k *= 2;
		}
	} 
}

