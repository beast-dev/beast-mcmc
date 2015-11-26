/*
 * SimpleFileFilter.java
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

package dr.app.bss;

import java.io.File;

import javax.swing.filechooser.FileFilter;

/**
 * @author Filip Bielejec
 * @version $Id$
 */
public class SimpleFileFilter extends FileFilter {

	private String[] extensions;
	private String description;

	public SimpleFileFilter(String ext) {
		this(new String[] { ext }, null);
	}

	public SimpleFileFilter(String[] exts, String descr) {

		// clone and lowercase the extensions
		extensions = new String[exts.length];
		for (int i = exts.length - 1; i >= 0; i--) {
			extensions[i] = exts[i].toLowerCase();
		}

		// make sure we have a valid (if simplistic) description
		description = (descr == null ? exts[0] + " files" : descr);
	}

	public boolean accept(File f) {

		// we always allow directories, regardless of their extension
		if (f.isDirectory()) {
			return true;
		}

		// ok, it's a regular file so check the extension
		String name = f.getName().toLowerCase();
		for (int i = extensions.length - 1; i >= 0; i--) {
			if (name.endsWith(extensions[i])) {
				return true;
			}
		}
		return false;
	}

	public String getDescription() {
		return description;
	}

}// END: class
