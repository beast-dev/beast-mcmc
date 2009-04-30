/*
 * Utils.java
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

package dr.app.util;

import java.awt.*;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.lang.reflect.Constructor;
import java.util.*;
import java.util.List;

/**
 *
 * @author  adru001
 * @version
 */
public class Utils {
	
	public static String getLoadFileName(String message) {
		java.io.File file = getLoadFile(message);
		if (file == null) return null;
		return file.getAbsolutePath();
	}

	public static String getSaveFileName(String message) {
		java.io.File file = getSaveFile(message);
		if (file == null) return null;
		return file.getAbsolutePath();
	}

	public static File getLoadFile(String message) {
		// No file name in the arguments so throw up a dialog box...
		java.awt.Frame frame = new java.awt.Frame();
		java.awt.FileDialog chooser = new java.awt.FileDialog(frame, message,
															java.awt.FileDialog.LOAD);
		chooser.show();
		if (chooser.getFile() == null) return null;
		java.io.File file = new java.io.File(chooser.getDirectory(), chooser.getFile());
		chooser.dispose();
		frame.dispose();

		return file;
	}

	public static File getSaveFile(String message) {
		// No file name in the arguments so throw up a dialog box...
		java.awt.Frame frame = new java.awt.Frame();
		java.awt.FileDialog chooser = new java.awt.FileDialog(frame, message,
															java.awt.FileDialog.SAVE);
		chooser.show();
		java.io.File file = new java.io.File(chooser.getDirectory(), chooser.getFile());
		chooser.dispose();
		frame.dispose();

		return file;
	}

    // detect type of text value - return class of type
    public static Class detectTYpe(final String valueString) {
        if ( valueString.equalsIgnoreCase("TRUE") || valueString.equalsIgnoreCase("FALSE")) {
            return Boolean.class;
        }

        try {
            final double number = Double.parseDouble(valueString);
            if( Math.round(number) == number ) {
                return Integer.class;
            }
            return Double.class;
        } catch (NumberFormatException pe) {
            return String.class;
        }
    }

    // Ner object of type cl from text.
    // return null if can't be done of value can't be converted

    public static Object constructFromString(Class cl, String value) {
        for( Constructor c : cl.getConstructors() ) {
            final Class[] classes = c.getParameterTypes();
            if (classes.length == 1 && classes[0].equals(String.class)) {
                try {
                    return c.newInstance(value);
                } catch( Exception e ) {
                    return null;
                }
            }
        }
        return null;
    }

    // skip over comment and empty lines

	private static String nextNonCommentLine(BufferedReader reader) throws IOException {
        String line;
        do {
            line = reader.readLine();
            // ignore empty or comment lines
        } while( line != null && (line.trim().length() == 0 || line.trim().charAt(0) == '#') );
        return line;
    }

    /**
     * Load traits from file.
     * @param file
     * @param delimiter
     * @return  A map whose key is the trait. The value is a list of <taxa, value> as a string array of size 2.
     * @throws IOException
     */
    public static Map<String, List<String[]>> importTraitsFromFile(File file, final String delimiter)
            throws IOException, Arguments.ArgumentException {
        final BufferedReader reader = new BufferedReader(new FileReader(file));

        String line = nextNonCommentLine(reader);

        // first column is label for the redundant "taxa" name
        final int startAt = 1;

        final String[] labels = line.split(delimiter);
        for(int k = 0; k < labels.length; ++k) {
            labels[k] = labels[k].trim();
        }

        Map<String, java.util.List<String[]>> traits = new HashMap<String, java.util.List<String[]>>();
        for(int i = startAt; i < labels.length; i++) {
            traits.put(labels[i], new ArrayList<String[]>());
        }

        line = nextNonCommentLine(reader);
        while (line != null) {
            String[] values = line.split(delimiter);

            assert ( values.length > 0 );

            final String first = values[0].trim();
            int k = Arrays.asList(labels).indexOf(first);
            if( k >= 0 ) {
                java.util.List<String[]> trait = traits.get(first);
                if( trait == null ) {
                    throw new Arguments.ArgumentException("undefined trait " + first);
                }
                final String traitVal = values[1].trim();
                for(int i = 2; i < values.length; i++) {
                    trait.add(new String[]{values[i], traitVal});
                }
            } else {
                for(int i = startAt; i < values.length; i++) {
                    if( i < labels.length ) {
                        java.util.List<String[]> column = traits.get(labels[i]);
                        column.add(new String[]{first, values[i].trim()});
                    }
                }
            }

            line = nextNonCommentLine(reader);
        }
        return traits;
    }

	/**
	 * This function takes a file name and an array of extensions (specified
	 * without the leading '.'). If the file name ends with one of the extensions
	 * then it is returned with this trimmed off. Otherwise the file name is
	 * return as it is.
	 * @return the trimmed filename
	 */
	public static String trimExtensions(String fileName, String[] extensions) {

		String newName = null;

        for( String extension : extensions ) {
            final String ext = "." + extension;
            if( fileName.toUpperCase().endsWith(ext.toUpperCase()) ) {
                newName = fileName.substring(0, fileName.length() - ext.length());
            }
        }

		return ( newName != null ) ? newName : fileName;
	}

	/**
	 * @return a named image from file or resource bundle.
	 */
	public static Image getImage(Object caller, String name) {

		java.net.URL url = caller.getClass().getResource(name);
		if (url != null) {
			return Toolkit.getDefaultToolkit().createImage(url);
		} else {
			if (caller instanceof Component) {
				Component c = (Component)caller;
				Image i = c.createImage(100,20);
				Graphics g = c.getGraphics();
				g.drawString("Not found!", 1, 15);
				return i;
			} else return null;
		}
	}
}
