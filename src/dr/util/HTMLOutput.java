/*
 * HTMLOutput.java
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

package dr.util;

public class HTMLOutput {

	
	public static final String HTML_OPEN = "<html>\n";
	
	public static final String HTML_HEAD =
		"  <head>\n" +			
		"    <style type=\"text/css\">\n" +
		"      body { \n" +
		"        margin-left: 10%; \n" +
		"        margin-right: 10%; \n" +
		"        font-family: sans-serif\n" +
		"      }\n" +
		"      h2 { margin-top: 1em }\n" +
		"      h2,h3,h4,h5,h6 { \n" +
		"        margin-left: -3%;\n" +
		"        font-family:Optima;\n" +
		"      }\n" +
		"      pre {\n" +
		"         color: green; font-weight: bold;\n" +
		"         white-space: pre; font-family: \"Courier New\", monospace;\n" +
		"      }\n" +
		"      tt { color: green }\n" +
		"      em { font-style: normal; font-weight: bold }\n" +
		"      strong { text-transform: uppercase; font-weight: bold }\n" +
		"      table {\n" +
		"        margin-left: -4%;\n" +
		"        font-family: sans-serif;\n" +
		"        background: white;\n" +
		"        border-width: 1px;\n" +
		"        border-color: white;\n" +
		"      }\n" +
		"    </style>\n" +
		"  </head>\n";
}