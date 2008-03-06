/*
 * BeastVersion.java
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

package dr.app.beast;

import dr.util.Version;

/**
 * This class provides a mechanism for returning the version number of the
 * dr software. It relies on the administrator of the dr source using the
 * module tagging system in CVS. The method getVersionString() will return
 * the version of dr under the following condition: <BR>
 * 1. the dr source has been checked out *by tag* before being packaged for
 * distribution.
 *
 * Version last changed 2007/07/31 by AER
 *
 * @author Alexei Drummond
 * @author Andrew Rambaut
 */
public class BeastVersion implements Version {


	/**
	 * Version string: assumed to be in format x.x.x
	 */
	private static String VERSION = "1.5a1";
	private static String DATE_STRING = "2002-2007";

	/**
     * this used to parse the CVS ID string but there is no equivalent
     * for SVN. We must increment this manually.
	 */
	private static String BUILD_ID = "Build r616";

	public String getVersionString() {
		return "v" + VERSION;
	}

	public String getDateString() {
		return DATE_STRING;
	}

    public String getBuildString() {
        return BUILD_ID;
	}
}
