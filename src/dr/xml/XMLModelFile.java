/*
 * XMLModelFile.java
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

package dr.xml;

import static dr.util.HeapSort.sort;
import org.jdom.Attribute;
import org.jdom.Document;
import org.jdom.Element;
import org.jdom.JDOMException;
import org.jdom.filter.Filter;
import org.jdom.input.SAXBuilder;
import org.jdom.output.XMLOutputter;

import javax.swing.*;
import javax.swing.event.ListDataListener;
import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.util.*;

/**
 * Package: XMLModelFile
 * Description:
 * <p/>
 * <p/>
 * Created by
 * Alexander V. Alekseyenko (alexander.alekseyenko@gmail.com)
 * Date: Apr 10, 2009
 * Time: 12:35:42 PM
 */
public class XMLModelFile implements ListModel {
    private final Element xmlModel;
    private final Map<String, XMLModelIdentifiable> identifiedElements;

    public XMLModelFile(Element root) {
        xmlModel = root; //doc.getRootElement();

        identifiedElements = new HashMap<String, XMLModelIdentifiable>();
        /* Parse out all the elements with id's and idref's
          Create a mapping of parameter names to the corresponding xml elements
          and a mapping of other identified elements to the corresponding xml elements
         */

        Iterator identifiedIter = xmlModel.getDescendants(new IdentifiedFilter());

        while (identifiedIter.hasNext()) {
            Element child = (Element) identifiedIter.next();
            Attribute idAttr = child.getAttribute("id");
            if (idAttr != null) {
                identifiedElements.put(idAttr.getValue(), new XMLModelIdentifiable(idAttr.getValue(), child));
            } else {
                idAttr = child.getAttribute("idref");
                identifiedElements.get(idAttr.getValue()).addReference(child);
            }
        }

        Iterator mcmcIter = xmlModel.getDescendants(new MCMCFilter());
        Element mcmc = (Element) mcmcIter.next();
        if (mcmc == null) {
            System.err.println("Error: cannot find mcmc element");
        } else {
            for (Object child : mcmc.getChildren("posterior")) {
                if ((child instanceof Element)) {
                    Element posteriorElement = (Element) child;
                }
            }
            mcmc.detach();
        }

    }

    public Set<String> getIdentifiedElementNames() {
        return identifiedElements.keySet();
    }

    public void prefixIdentifiedNames(String prefix, Map<String, String> special, boolean keepID) {
        if (special == null) {
            special = new HashMap<String, String>();
        }

        for (Map.Entry<String, XMLModelIdentifiable> stringXMLModelIdentifiableEntry : identifiedElements.entrySet()) {
            String newName;
            if (special.containsKey(stringXMLModelIdentifiableEntry.getKey())) {//assign a special name
                newName = special.get(stringXMLModelIdentifiableEntry.getKey());
                if (keepID) {
                    stringXMLModelIdentifiableEntry.getValue().restoreDefinition();
                } else {
                    stringXMLModelIdentifiableEntry.getValue().removeDefinition();
                }
            } else {
                newName = prefix + stringXMLModelIdentifiableEntry.getKey();
            }
            stringXMLModelIdentifiableEntry.getValue().rename(newName);
        }
    }

    public void printIdentified() {
        System.out.println("Identified elements follow");
        XMLOutputter outputter = new XMLOutputter();
        for (Map.Entry<String, XMLModelIdentifiable> stringXMLModelIdentifiableEntry : identifiedElements.entrySet()) {
            System.out.println("Original name: " + stringXMLModelIdentifiableEntry.getKey());
            stringXMLModelIdentifiableEntry.getValue().print(outputter, System.out);
        }
    }

    public void print(XMLOutputter outputter, OutputStream ostream) {
        try {
            outputter.output(xmlModel.getContent(), ostream);
        }
        catch (IOException e) {
            System.err.println("Error writing model!");
        }
    }


    public class MCMCFilter implements Filter {
        public boolean matches(Object o) {
            return (o instanceof Element) && ((Element) o).getName().equals("mcmc");
        }
    }

    public class IdentifiedFilter implements Filter {

        public boolean matches(Object o) {
            if (!(o instanceof Element)) {
                return false;
            }
            return ((Element) o).getAttribute("id") != null || ((Element) o).getAttribute("idref") != null;
        }
    }


    // ListModel implementation
    public int getSize() {
        return identifiedElements.size();
    }

    public Object getElementAt(int i) {
        if (i < identifiedElements.size()) {
            Object[] elements = identifiedElements.keySet().toArray();
            sort(elements, new Comparator() {
                public int compare(Object o, Object o1) {
                    return o.toString().compareTo(o1.toString());
                }
            });
            return elements[i];
        }
        return null;
    }

    public void addListDataListener(ListDataListener listDataListener) {
        //AUTOGENERATED METHOD IMPLEMENTATION
    }

    public void removeListDataListener(ListDataListener listDataListener) {
        //AUTOGENERATED METHOD IMPLEMENTATION
    }


    public static void main(String[] args) {  //main method for debugging
        SAXBuilder parser = new SAXBuilder();

        Document doc;
        try {
            doc = parser.build(new File("testSimplePathSampling.xml"));

            XMLModelFile z = new XMLModelFile(doc.getRootElement());

            z.printIdentified();
            HashMap<String, String> hm = new HashMap<String, String>();
            hm.put("samplingMean", "samplingMean");
            z.prefixIdentifiedNames("model1.", hm, false);
            z.printIdentified();
        } catch (IOException e) {
          //
        }
        catch (JDOMException e) {
         //
        }
    }

}
