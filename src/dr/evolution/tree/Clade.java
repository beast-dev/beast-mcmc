/*
 * Clade.java
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

/**
 * 
 */
package dr.evolution.tree;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;

import dr.math.MathUtils;
import dr.math.distributions.Distribution;
import dr.math.distributions.NormalDistribution;

/**
 * @author Sebastian Hoehna
 *
 */
public class Clade implements Comparable<Clade> {

    private final BitSet bits;
    private double height;
    private List<Double> heights;
    private final int size;
    
    public Clade(final BitSet bits, final double height) {
        this.bits = bits;
        this.height = height;
        size = bits.cardinality();
        heights = new ArrayList<Double>();
    }

    public BitSet getBits() {
        return bits;
    }

    public double getHeight() {
        return height;
    }
    
    public double getHeight(int index) {
        return heights.get(index);
    }
    
    public int getSampleCount(){
    	return heights.size();
    }

    public int getSize() {
        return size;
    }

    public int compareTo(Clade clade) {
		int setBitIndexI = -1;
		int setBitIndexJ = -1;
		BitSet otherBits = clade.getBits();
		do {
			setBitIndexI = bits.nextSetBit(setBitIndexI + 1);
			setBitIndexJ = otherBits.nextSetBit(setBitIndexJ + 1);
		} while (setBitIndexI == setBitIndexJ && setBitIndexI != -1);

		return (setBitIndexI < setBitIndexJ ? -1
				: (setBitIndexI > setBitIndexJ ? 1 : 0));
	}
    
    public boolean equals(Object o){
    	if (o instanceof Clade){
    		return equals((Clade)o);
    	}
    	return false;
    }
    
    public boolean equals(Clade c){
    	if (size != c.getSize()){
    		return false;
    	}
    	return compareTo(c) == 0;
    }
    
    public String toString(){
    	return bits.toString();
    }
    
    public int hashCode(){
    	return bits.hashCode();
    }
    
    public void addHeight(double height){
    	heights.add(height);
    	setHeight(height);
    }
    
    public void setHeight(double height){
    	this.height = height;
    }

    public static Set<Clade> getCladeSet(Tree tree) {
        return new HashSet<Clade>(getCladeList(tree));
    }

    public static List<Clade> getCladeList(Tree tree) {
        List<Clade> clades = new ArrayList<Clade>();
        getClades(clades, tree, tree.getRoot());
        return clades;
    }

    private static BitSet getClades(List<Clade> clades, Tree tree, NodeRef node) {

        BitSet bits = new BitSet();

        if (tree.isExternal(node)) {

            bits.set(node.getNumber());

        } else {
            bits.or(getClades(clades, tree, tree.getChild(node, 0)));
            bits.or(getClades(clades, tree, tree.getChild(node, 1)));

            clades.add(new Clade(bits, tree.getNodeHeight(node)));
        }

        return bits;
    }

}
