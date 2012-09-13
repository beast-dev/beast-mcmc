/*
 * CompositePrior.java
 *
 * Copyright (C) 2002-2006 Alexei Drummond and Andrew Rambaut
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

package dr.inference.prior;

/**
 * This prior provides a composition of a number of simple priors.
 */
//public class CompositePrior implements Prior {
//
//	Vector priors = null;
//	String name = null;
//
//	public CompositePrior(String name) {
//		priors = new Vector();
//		this.name = name;
//	}
//
//	public CompositePrior(String name, Prior[] priorArray) {
//
//		this.name = name;
//		priors = new Vector(priorArray.length);
//		for (int i =0; i < priorArray.length; i++) {
//			priors.add(priorArray[i]);
//		}
//	}
//
//	public double getLogPrior(dr.inference.model.Model model) {
//
//		if (priors == null) return 0.0;
//
//		double logPrior = 0.0;
//
//		for (int i = 0; i < priors.size(); i++) {
//			Prior p = (Prior)priors.elementAt(i);
//			double l = p.getLogPrior(model);
//			logPrior += l;
//		}
//
//		return logPrior;
//	}
//
//	public void addPrior(Prior p) {
//		if (p != null) priors.add(p);
//	}
//
//	public Iterator getPriorIterator() {
//		return priors.iterator();
//	}
//
//	public String getPriorName() {
//		return name;
//	}
//
//	public Element createElement(Document d) {
//		throw new RuntimeException("Not implemented!");
//	}
//
//	public String toString() {
//		StringBuffer buffer = new StringBuffer();
//		for (int i = 0; i < priors.size(); i++) {
//			Prior p = (Prior)priors.elementAt(i);
//			buffer.append(p.toString());
//		}
//		return buffer.toString();
//	}
//}
