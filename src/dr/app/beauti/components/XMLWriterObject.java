/*
 * XMLWriterObject.java
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

package dr.app.beauti.components;

import dr.app.beauti.util.XMLWriter;
import dr.util.Attribute;

public class XMLWriterObject {

    private final String tag;
    private final String id;
    private final String content;
    private final XMLWriterObject[] children;
    private final Attribute[] attributes;
    private boolean hasBeenWritten = false;

    public XMLWriterObject(String tag, String id, XMLWriterObject[] children, Attribute[] attributes, String content) {
        this.tag = tag;
        this.id = id;
        this.content = content;

        if (children == null) children = new XMLWriterObject[0];
        this.children = children;

        if (attributes == null) attributes = new Attribute[0];
        this.attributes = attributes;

    }


    public XMLWriterObject(String tag, XMLWriterObject child) {
        this(tag, null, new XMLWriterObject[]{child}, null, null);
    }

    public XMLWriterObject(String tag, String id) {
        this(tag, id, null, null, null);
    }


    public XMLWriterObject(String tag, String id, XMLWriterObject[] children, Attribute[] attributes) {
        this(tag, id, children, attributes, null);
    }

    public XMLWriterObject(String tag, String id, XMLWriterObject[] children) {
        this(tag, id, children, null, null);
    }

    public static XMLWriterObject refObj(String tag, String id) {
        XMLWriterObject obj = new XMLWriterObject(tag, id);
        obj.setAlreadyWritten(true);
        return obj;
    }

    public void setAlreadyWritten(boolean b) {
        this.hasBeenWritten = b;
    }

    public void writeOrReference(XMLWriter writer) {
        if (hasBeenWritten) {
            reference(writer);
        } else {
            Attribute[] finalAttributes;
            if (id == null) {
                finalAttributes = attributes;
            } else {
                finalAttributes = new Attribute[attributes.length + 1];
                finalAttributes[0] = new Attribute.Default("id", id);
                System.arraycopy(attributes, 0, finalAttributes, 1, attributes.length);
            }

            if (children.length == 0) {
                writer.writeTag(tag, finalAttributes, true);
            } else {
                writer.writeOpenTag(tag, finalAttributes);
                for (XMLWriterObject child : children) {
                    child.writeOrReference(writer);
                }
                writer.writeCloseTag(tag);
            }

            hasBeenWritten = true;
        }
    }

    public void reference(XMLWriter writer) {
        writer.writeTag(tag, new Attribute.Default("idref", id), true);
    }
}


