package edu.washington.gs.noble.crux.gui;

import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.util.logging.Logger;

import javax.swing.BorderFactory;
import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JToggleButton;
import javax.swing.border.EtchedBorder;

@SuppressWarnings("serial")
class QRankerPanel extends JPanel implements ItemListener {

	private static Logger logger = 
		Logger.getLogger("edu.washington.gs.noble.crux.gui");
	private final JCheckBox runToolCheckBox = new JCheckBox("Run this tool");
	private final JCheckBox showAdvancedParameters = new JCheckBox("Show advanced parameters");
	private final JButton saveButton = new JButton("Save");
	private final JButton cancelButton = new JButton("Cancel");
	private final JButton loadDefaults = new JButton("Load Defaults");
	private CruxAnalysisModel model = null;
	private CruxComponentButton button;
	private JToggleButton dummyButton;

	public QRankerPanel(CruxAnalysisModel model, final CruxComponentButton button, final JToggleButton dummy) {
		super();
		this.model = model;
		this.button = button;
		this.dummyButton = dummy;
		setBorder(
			BorderFactory.createTitledBorder(
				BorderFactory.createEtchedBorder(EtchedBorder.RAISED), 
				"q-ranker parameters"
			)
		);
		setBackground(Color.white);
		setLayout(new BoxLayout(this, BoxLayout.PAGE_AXIS));
		runToolCheckBox.setBackground(Color.white);
		showAdvancedParameters.setBackground(Color.white);
		showAdvancedParameters.setEnabled(false);
		JPanel buttonPanel = new JPanel();
		buttonPanel.setLayout(new BoxLayout(buttonPanel, BoxLayout.LINE_AXIS));
		buttonPanel.setBackground(Color.white);
		buttonPanel.setAlignmentX(Component.LEFT_ALIGNMENT);
		//saveButton.setEnabled(false);
		saveButton.addActionListener(new SaveButtonListener());
		buttonPanel.add(Box.createRigidArea(new Dimension(36,0)));
		buttonPanel.add(saveButton);
		cancelButton.addActionListener(new CancelButtonListener());
		buttonPanel.add(Box.createRigidArea(new Dimension(12,0)));
		buttonPanel.add(cancelButton);
		buttonPanel.add(Box.createRigidArea(new Dimension(12,0)));
		loadDefaults.setEnabled(false);
		buttonPanel.add(loadDefaults);
		add(runToolCheckBox);
		add(Box.createRigidArea(new Dimension(0,12)));
		add(showAdvancedParameters);
		add(Box.createRigidArea(new Dimension(0,12)));
		add(buttonPanel);
		runToolCheckBox.addItemListener(new RunComponentChangeListener());
		updateFromModel();
		setVisible(false);
	}
	
	private void updateFromModel() {
		runToolCheckBox.setSelected(model.getRunSearchForMatches());
	}
	
	class RunComponentChangeListener implements ItemListener {
		public void itemStateChanged(final ItemEvent event) {
			// saveButton.setEnabled(runToolCheckBox.isSelected());
		}
	}
	
	public void itemStateChanged(final ItemEvent event) {
		logger.info("User selected compute-q-values component.");
		updateFromModel();
		if (((CruxComponentButton) event.getSource()).isSelected() == true) {
			setVisible(true);
		}
		else {
			setVisible(false);
		}
	}
	
	class SaveButtonListener implements ActionListener {
		public void actionPerformed(ActionEvent e) {
			boolean checked = runToolCheckBox.isSelected();
		}
	}
	
	class CancelButtonListener implements ActionListener {
		public void actionPerformed(ActionEvent e) {
			setVisible(false);			
			dummyButton.setSelected(true);
		}
	}
}