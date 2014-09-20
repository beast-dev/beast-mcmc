package dr.app.bss;

import java.awt.Component;
import java.util.HashSet;
import java.util.Set;

import javax.swing.JComboBox;
import javax.swing.JList;
import javax.swing.UIManager;
import javax.swing.plaf.basic.BasicComboBoxRenderer;

@SuppressWarnings("serial")
public class DisabledItemsComboBox extends JComboBox {

	private Set<Object> disabledItems = new HashSet<Object>();

	@SuppressWarnings("unchecked")
	public DisabledItemsComboBox() {
		super();
		super.setRenderer(new DisabledItemsRenderer());
	}// END: Constructor

	public void addItem(Object object, boolean isDisabled) {
		super.addItem(object);
		if (isDisabled) {
			disabledItems.add(getItemCount() - 1);
		}
	}// END: addItem

	@Override
	public void removeAllItems() {
		super.removeAllItems();
		disabledItems = new HashSet<Object>();
	}// END: removeAllItems

	@Override
	public void removeItemAt(final int index) {
		super.removeItemAt(index);
		disabledItems.remove(index);
	}// END: removeItemAt

	@Override
	public void removeItem(final Object anObject) {
		for (int i = 0; i < getItemCount(); i++) {

			if (getItemAt(i) == anObject) {
				disabledItems.remove(i);
			}
		}

		super.removeItem(anObject);
	}// END: removeItem

	@Override
	public void setSelectedIndex(int index) {

		if (!disabledItems.contains(index)) {
			super.setSelectedIndex(index);
		}

	}// END: setSelectedIndex

	private class DisabledItemsRenderer extends BasicComboBoxRenderer {

		@Override
		public Component getListCellRendererComponent(@SuppressWarnings("rawtypes") JList list, Object value,
				int index, boolean isSelected, boolean cellHasFocus) {

			if (isSelected) {

				setBackground(list.getSelectionBackground());
				setForeground(list.getSelectionForeground());

			} else {

				setBackground(list.getBackground());
				setForeground(list.getForeground());
			}

			if (disabledItems.contains(index)) {

				setFocusable(false);
				setEnabled(false);
				setBackground(list.getBackground());
				setForeground(UIManager.getColor("Label.disabledForeground"));

			} else {

				setFocusable(true);
				setEnabled(true);

			}

			setText((value == null) ? "" : value.toString());

			return this;
		}// END: getListCellRendererComponent

	}// END: DisabledItemsRenderer

}// END: class
