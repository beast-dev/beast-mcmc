/*
 * XMLWriter.java
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

package dr.app.beauti.util;

import dr.util.Attribute;

import java.io.Writer;
import java.util.List;

/**
 * @author Alexei Drummond
 * @version $Id: XMLWriter.java,v 1.3 2005/06/27 21:18:40 rambaut Exp $
 */
public class XMLWriter extends java.io.PrintWriter {

    int level = 0;

    public XMLWriter(Writer writer) {
        super(writer);
    }

    public void increaseLevel() {
        level += 1;
    }

    public void decreaseLevel() {
        level -= 1;
    }

    public void writeComment(String... comments) {
        writeComment(80, comments);
    }

    public void writeComment(int length, String... comments) {
        writeBlankLine();

        for (String comment : comments) {
            StringBuffer buffer = new StringBuffer("<!-- ");
            buffer.append(comment);
            for (int i = buffer.length(); i < (length - 3); i++) {
                buffer.append(' ');
            }
            buffer.append("-->");
            writeText(buffer.toString());
        }
    }

    public void writeOpenTag(String tagname) {
        writeText("<" + tagname + ">");
        increaseLevel();
    }

    public void writeOpenTag(String tagname, Attribute attribute) {
        writeTag(tagname, new Attribute[]{attribute}, false);
    }

    public void writeOpenTag(String tagname, Attribute[] attributes) {
        writeTag(tagname, attributes, false);
    }

    public void writeOpenTag(String tagname, List<Attribute> attributes) {
        writeTag(tagname, attributes.toArray(new Attribute[attributes.size()]), false);
    }

    public void writeTag(String tagname, Attribute attribute, boolean close) {
        writeTag(tagname, new Attribute[]{attribute}, close);
    }

    public void writeTag(String tagname, Attribute[] attributes, boolean close) {
        StringBuffer buffer = new StringBuffer("<");
        buffer.append(tagname);
        for (Attribute attribute : attributes) {
            if (attribute != null) {
                buffer.append(' ');
                buffer.append(attribute.getAttributeName());
                buffer.append("=\"");
                buffer.append(attribute.getAttributeValue());
                buffer.append("\"");
            }
        }
        if (close) {
            buffer.append("/");
        }
        buffer.append(">");
        writeText(buffer.toString());
        if (!close) {
            increaseLevel();
        }
    }

    public void writeTag(String tagname, Attribute[] attributes, String content, boolean close) {
        StringBuffer buffer = new StringBuffer("<");
        buffer.append(tagname);
        for (Attribute attribute : attributes) {
            if (attribute != null) {
                buffer.append(' ');
                buffer.append(attribute.getAttributeName());
                buffer.append("=\"");
                buffer.append(attribute.getAttributeValue());
                buffer.append("\"");
            }
        }
        if (content != null) {
            buffer.append(">");
            buffer.append(content);
            if (close) {
                buffer.append("</");
                buffer.append(tagname);
                //buffer.append("/");
            }
        } else if (close) {
            buffer.append("/");
        }
        buffer.append(">");
        writeText(buffer.toString());
        if (!close) {
            increaseLevel();
        }
    }

    public void writeCloseTag(String tagname) {
        decreaseLevel();
        writeText("</" + tagname + ">");
    }

    public void writeBlankLine() {
        println();
        flush();
    }

    public void writeText(String string) {
        for (int i = 0; i < level; i++) {
            write('\t');
        }
        println(string);

        // wrapping lines at a fixed length can only do bad things...
//        int subLength = 2000;
//
//        for (int l = 0; l < string.length(); l += subLength) {
//            int e = l + subLength;
//            if (e > string.length()) e = string.length();
//            for (int i = 0; i < level; i++) {
//                write('\t');
//            }
//            println(string.substring(l, e));
//            flush();
//        }

    }

    public void writeIDref(String tagName, String paramName) {
        writeTag(tagName, new Attribute[]{new Attribute.Default<String>(dr.xml.XMLParser.IDREF, paramName)}, true);
    }
}

