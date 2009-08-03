package dr.xml;

import org.jdom.Element;
import org.jdom.output.XMLOutputter;

import java.io.IOException;
import java.io.OutputStream;
import java.util.LinkedList;

/**
 * Package: XMLModelIdentifiable
 * Description:
 * <p/>
 * <p/>
 * Created by
 * Alexander V. Alekseyenko (alexander.alekseyenko@gmail.com)
 * Date: Apr 13, 2009
 * Time: 10:38:05 AM
 */
public class XMLModelIdentifiable {
    private String name;
    private Element definition;
    private Element definitionIdref;
    private boolean defined;
    private LinkedList<Element> references;

    public XMLModelIdentifiable(String name, Element definition) {
        this.name = name;
        this.definition = definition;
        definitionIdref = new Element(definition.getName());

        definitionIdref.setAttribute("idref", name);

        defined = true;
        references = new LinkedList<Element>();
    }

    public void removeDefinition() {
        if (defined) {
            definition.getParentElement().addContent(definitionIdref);
            definition.detach();
            defined = false;
        }
    }

    public void restoreDefinition() {
        if (!defined) {
            definitionIdref.getParentElement().addContent(definition);
            definitionIdref.detach();
            defined = true;
        }

    }

    public void rename(String newName) {
        name = newName;
        definition.setAttribute("id", newName);
        definitionIdref.setAttribute("idref", newName);
        for (Element ref : references) {
            ref.setAttribute("idref", newName);
        }

    }

    public void addReference(Element newRef) {
        if (newRef.getAttribute("idref").getValue().equals(name)) {
            references.addLast(newRef);
        }
    }

    public void print(XMLOutputter outputter, OutputStream ostream) {
        if (ostream == null) {
            ostream = System.out;
        }
        try {
            if (defined) {
                outputter.output(definition, ostream);
            } else {
                outputter.output(definitionIdref, ostream);
            }
            for (Element ref : references) {
                outputter.output(ref, ostream);
            }

        } catch (IOException e) {
            System.err.println(e);
        }
    }
}
