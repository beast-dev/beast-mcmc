/*
 * BEASTInputFile.java
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

package test.dr.integration;

/**
 * GTR Parameter Estimation Tester.
 *
 * @author Walter Xie
 * @since <pre>08/06/2009</pre>
 */
public class BEASTInputFile {

	public static final double ct = 1; // fixed	
	public static final double birthRate = 50;
	
	private double ac;
	private double ag;
	private double at;
	private double cg;
	private double gt;
    private String fileNamePrefix;    

    
	public double getAc() {
		return ac;
	}
	public void setAc(double ac) {
		this.ac = ac;
	}
	public double getAg() {
		return ag;
	}
	public void setAg(double ag) {
		this.ag = ag;
	}
	public double getAt() {
		return at;
	}
	public void setAt(double at) {
		this.at = at;
	}
	public double getCg() {
		return cg;
	}
	public void setCg(double cg) {
		this.cg = cg;
	}
	public double getGt() {
		return gt;
	}
	public void setGt(double gt) {
		this.gt = gt;
	}
	public String getFileNamePrefix() {
		return fileNamePrefix;
	}
	public void setFileNamePrefix (String fileNamePrefix) {
		this.fileNamePrefix = fileNamePrefix;
	}

 
}
