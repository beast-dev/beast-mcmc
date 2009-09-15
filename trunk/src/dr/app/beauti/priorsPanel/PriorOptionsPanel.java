package dr.app.beauti.priorsPanel;

import dr.app.beauti.options.Parameter;
import dr.math.distributions.Distribution;
import org.virion.jam.components.RealNumberField;
import org.virion.jam.panels.OptionsPanel;

import javax.swing.*;
import java.util.ArrayList;
import java.util.List;

/**
 * @author Alexei Drummond
 */
public abstract class PriorOptionsPanel extends OptionsPanel {

    private List<JComponent> argumentFields = new ArrayList<JComponent>();
    private List<String> argumentNames = new ArrayList<String>();

    public PriorOptionsPanel() {

        super(12, 12);
    }

    public PriorOptionsPanel(String[] names, RealNumberField[] fields) {
        for (int i = 0; i < names.length; i++) {
            addField(names[i], fields[i]);
        }
        setupComponents();
    }

    void addField(String name, double initialValue, double min, double max) {

        RealNumberField field = new RealNumberField(min, max);
        field.setValue(initialValue);
        addField(name, field);
    }

    void addField(String name, RealNumberField field) {
        argumentNames.add(name);

        field.setColumns(8);
        argumentFields.add(field);
        setupComponents();
    }

    void addCheckBox(String name, JCheckBox jCheckBox) {
        argumentNames.add(name);

        argumentFields.add(jCheckBox);
        setupComponents();
    }

    void replaceFieldName(int i, String name) {
        argumentNames.set(i, name);
        setupComponents();
    }

    double getValue(int i) {
        return ((RealNumberField) argumentFields.get(i)).getValue();
    }

    void setupComponents() {
        removeAll();

        for (int i = 0; i < argumentFields.size(); i++) {
            addComponentWithLabel(argumentNames.get(i) + ":", argumentFields.get(i));
        }
    }

    public abstract Distribution getDistribution();

    /**
     * Set the prior on the given parameter
     *
     * @param parameter the parameter to set this prior on
     */
    public abstract void setParameterPrior(Parameter parameter);

    RealNumberField getField(int i) {
        return (RealNumberField) argumentFields.get(i);
    }

    public List<JComponent> getJComponents() {
        return argumentFields;
    }
}
