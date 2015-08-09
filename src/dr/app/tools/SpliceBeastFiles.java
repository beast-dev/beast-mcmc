/*
 * SpliceBeastFiles.java
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

package dr.app.tools;

import java.io.*;
import java.util.*;

/**
 * Simple program to splice a set of (probably simulated) beast files containting sequences
 * with a fixed second half containing some analysis setup. The first half, up to but not
 * including the line '\t<patterns id="patterns" from="1">' comes from a set of XML numbered
 * from 1 to 100, the second half is read from a specified file and a new set numbered 1 to
 * 100 is created.
 *
 * @author Andrew Rambaut
 * @version $Id$
 */
public class SpliceBeastFiles {

    public SpliceBeastFiles(String templateFileName, String inputFileStem, String outputFileStem) {
        try {
            List<String> template = new ArrayList<String>();

            BufferedReader reader = new BufferedReader(new FileReader(templateFileName));
            String line = reader.readLine();
            while (line != null) {
                template.add(line);
                line = reader.readLine();
            }

            reader.close();

            for (int i = 1; i < 201; i++) {
                String number = (i < 10? "00" : (i < 100? "0" : "")) + i;
                reader = new BufferedReader(new FileReader(inputFileStem + number + ".xml"));

                PrintWriter writer = new PrintWriter(outputFileStem + number + ".xml");

                line = reader.readLine();
                while (line != null && !line.equals("\t<patterns id=\"patterns\" from=\"1\">")) {
                    writer.println(line);
                    line = reader.readLine();
                }

                for (String line1 : template) {
                    if (line1.contains("$stem")) {
                        line1 = line1.replace("$stem", "sim" + number);
                    }
                    writer.println(line1);
                }

                writer.close();

                reader.close();
            }
        } catch (FileNotFoundException e) {
            e.printStackTrace();
            return;
        } catch (IOException e) {
            e.printStackTrace();
            return;
        }
    }

    public static void main(String[] argv) {
        new SpliceBeastFiles(argv[0], argv[1], argv[2]);
    }
}