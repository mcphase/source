//package mcphas;

import java.awt.*;
import java.awt.event.*;
import java.applet.*;

//import com.sun.java.swing.*;
import javax.swing.*;
import borland.jbcl.control.*;
import borland.jbcl.layout.*;

/******************************************
 * class declaration
 ******************************************/
public class IniConfigurator extends JPanel 
{

  String CONST_MCPHASE_RUNTIME_CONTROL = "MCPHASE RUNTIME CONTROL";
  String CONST_XY_PHASEDIAGRAM_PARAMETERS = "XY PHASEDIAGRAM PARAMETERS";
  String CONST_GENERATION_OF_SPIN_CONFIGURATIONS = "GENERATION OF SPIN CONFIGURATIONS";
  String CONST_PARAMETERS_FOR_SUB_FECALC_SELFCONSISTENCY_PROCESS = "PARAMETERS FOR SUB FECALC SELFCONSISTENCY PROCESS";
  String CONST_OUTPUT_OF_PHYSICAL_PROPERTIES = "OUTPUT OF PHYSICAL PROPERTIES";
	
  // Member variables		
  String strWorkingDir;					// contains the working directory of the program
  String strFileSeparator;
  IniFile m_IniFile = new IniFile();	// class for handling file transactions
  boolean isStandalone = false;

  /*************************************************
   * elements of the form
   *************************************************/
  
  // Button-Panel
  JPanel panOkCancel = new JPanel();
  
  // Buttons
  JButton btnRead = new JButton();  
  JButton btnWrite = new JButton();
  
  // Checkboxes
  JCheckBox chkExit = new JCheckBox();
  JCheckBox chkPause = new JCheckBox();
  JCheckBox chkDisplayAll = new JCheckBox();
  JCheckBox chkLogFevsQ = new JCheckBox();
  
  // Textfield for Status-display
  JTextField txStatus = new JTextField();
  
  // Tab-Panel
  JTabbedPane panMain = new JTabbedPane();
  // 1. Tab
  JPanel panXYPhaseParams = new JPanel();	// [XY PHASEDIAGRAM PARAMETERS]
	// Labels
	JLabel labXT = new JLabel();
	JLabel labXHa = new JLabel();
	JLabel labXHb = new JLabel();
	JLabel labXHc = new JLabel();
	JLabel labXMin = new JLabel();
	JLabel labXMax = new JLabel();
	JLabel labXStep = new JLabel();
	JLabel labXMax1 = new JLabel();
	JLabel labXMin1 = new JLabel();
	JLabel labXT1 = new JLabel();
	JLabel labXHa1 = new JLabel();
	JLabel labXHb1 = new JLabel();
	JLabel labXHc1 = new JLabel();
	JLabel labXStep1 = new JLabel();
        JLabel labT0 = new JLabel();
        JLabel labHa0 = new JLabel();
        JLabel labHb0 = new JLabel();
        JLabel labHc0 = new JLabel();
	
	// Textfields
	JTextField txXT = new JTextField();
	JTextField txXHa = new JTextField();
	JTextField txXHb = new JTextField();
	JTextField txXHc = new JTextField();
	JTextField txXMin = new JTextField();
	JTextField txXMax = new JTextField();
	JTextField txXStep = new JTextField();
	JTextField txYMin = new JTextField();
	JTextField txYT = new JTextField();
	JTextField txYHa = new JTextField();
	JTextField txYHb = new JTextField();
	JTextField txYHc = new JTextField();
	JTextField txYStep = new JTextField();
	JTextField txYMax = new JTextField();
	JTextField txT0 = new JTextField();
	JTextField txHa0 = new JTextField();
	JTextField txHb0 = new JTextField();
	JTextField txHc0 = new JTextField();
  
  // 2. Tab
  JPanel panQVector = new JPanel();			// [GENERATION OF SPIN-CONFIGURATIONS]
    // Labels
	JLabel labDeltaH = new JLabel();
	JLabel labHMax = new JLabel();
	JLabel labHMin = new JLabel();
	JLabel labDeltaL = new JLabel();
	JLabel labLMax = new JLabel();
	JLabel labLMin = new JLabel();
	JLabel labDeltaK = new JLabel();
	JLabel labKMax = new JLabel();
	JLabel labKMin = new JLabel();
	JLabel labMaxNofSpins = new JLabel();
	JLabel labNofRndTries = new JLabel();
	JLabel labMaxQPeriod = new JLabel();
  	JLabel labmaxnoftestspincf = new JLabel();
  	// Textfields
	JTextField txDeltaH = new JTextField();
	JTextField txHMax = new JTextField();
	JTextField txHMin = new JTextField();
	JTextField txDeltaK = new JTextField();
	JTextField txKMax = new JTextField();
	JTextField txKMin = new JTextField();
	JTextField txDeltaL = new JTextField();
	JTextField txLMax = new JTextField();
	JTextField txLMin = new JTextField();
	JTextField txMaxNofSpins = new JTextField();
	JTextField txNofRndTries = new JTextField();
	JTextField txMaxQPeriod = new JTextField();
	JTextField txmaxnoftestspincf = new JTextField();
  // 3. Tab		
  JPanel panAdditionalParams = new JPanel();	// [PARAMETERS FOR SUB FECALC - SELFCONSISTENCY PROCESS]
	// Labels
	JLabel labMaxSpinChange = new JLabel();
	//JLabel labSmallStep = new JLabel();
	JLabel labBigStep = new JLabel();
	JLabel labMaxStaMf = new JLabel();
	JLabel labMaxNofMfLoops = new JLabel();
	// Textfields
	JTextField txMaxSpinChange = new JTextField();
	//JTextField txSmallStep = new JTextField();
	JTextField txBigStep = new JTextField();
	JTextField txMaxStaMf = new JTextField();
	JTextField txMaxNofMfLoops = new JTextField();

  // 4. Tab		
  JPanel panOutputOfProperties = new JPanel();	// [PARAMETERS FOR SUB FECALC - SELFCONSISTENCY PROCESS]
	// Labels
	JLabel labMaxNofHkls = new JLabel();
	JLabel labNofSpinCorrs = new JLabel();
	JLabel labMaxQ = new JLabel();
	// Textfields
	JTextField txMaxNofHkls = new JTextField();
	JTextField txNofSpinCorrs = new JTextField();
	JTextField txMaxQ = new JTextField();

  // Helper Panels for display on the Tab-panles
  JPanel jPanel1 = new JPanel();
  JPanel jPanel2 = new JPanel();
  JPanel jPanel3 = new JPanel();
  JPanel jPanel4 = new JPanel();
  JPanel jPanel4a = new JPanel();
  JPanel jPanel4b = new JPanel();
  JPanel jPanel5 = new JPanel();
  JPanel jPanel6 = new JPanel();
  JPanel jPanel7 = new JPanel();
  JPanel jPanel8 = new JPanel();
  JPanel jPanel9 = new JPanel();
  JPanel jPanel10 = new JPanel();
  JPanel jPanel11 = new JPanel();
  JPanel jPanel12 = new JPanel();
  JPanel jPanel13 = new JPanel();
  JPanel jPanel14 = new JPanel();
  JPanel jPanel15 = new JPanel();
  JPanel jPanel16 = new JPanel();
  
  JPanel jPanel17 = new JPanel();
  JPanel jPanel18 = new JPanel();
  
  // Layouts
  FlowLayout flowLayout1 = new FlowLayout();
  VerticalFlowLayout verticalFlowLayout1 = new VerticalFlowLayout();
  VerticalFlowLayout verticalFlowLayout2 = new VerticalFlowLayout();
  VerticalFlowLayout verticalFlowLayout3 = new VerticalFlowLayout();
  VerticalFlowLayout verticalFlowLayout4 = new VerticalFlowLayout();
  VerticalFlowLayout verticalFlowLayout4a = new VerticalFlowLayout();
  VerticalFlowLayout verticalFlowLayout4b = new VerticalFlowLayout();
  FlowLayout flowLayout2 = new FlowLayout();
  VerticalFlowLayout verticalFlowLayout5 = new VerticalFlowLayout();
  VerticalFlowLayout verticalFlowLayout6 = new VerticalFlowLayout();
  VerticalFlowLayout verticalFlowLayout7 = new VerticalFlowLayout();
  VerticalFlowLayout verticalFlowLayout8 = new VerticalFlowLayout();
  VerticalFlowLayout verticalFlowLayout9 = new VerticalFlowLayout();
  VerticalFlowLayout verticalFlowLayout10 = new VerticalFlowLayout();
  FlowLayout flowLayout3 = new FlowLayout();
  VerticalFlowLayout verticalFlowLayout11 = new VerticalFlowLayout();
  VerticalFlowLayout verticalFlowLayout12 = new VerticalFlowLayout();
  VerticalFlowLayout verticalFlowLayout13 = new VerticalFlowLayout();
  VerticalFlowLayout verticalFlowLayout14 = new VerticalFlowLayout();
  VerticalFlowLayout verticalFlowLayout15 = new VerticalFlowLayout();
  FlowLayout flowLayout4 = new FlowLayout();
  FlowLayout flowLayout5 = new FlowLayout();
  VerticalFlowLayout verticalFlowLayout16 = new VerticalFlowLayout();
  VerticalFlowLayout verticalFlowLayout17 = new VerticalFlowLayout();


  /***************************************************************
   * Initializing function
   * All elements of the form are initialized
   ***************************************************************/
  private void jbInit() throws Exception 
  {
    //this.setSize(2000,2000);		// overall size of the form
	//this.setPreferredSize(new Dimension(2000,2000));	// preferred size
	//this.setMinimumSize(new Dimension(2000,2000));		// minimum size
	
	/**************************************
	 * Buttons (Read / Write)
	 **************************************/
	
	// Read-Button
    btnRead.setText("Read file mcphas.ini");		// display-text of the button
	btnRead.setToolTipText("Reads Data from file mcphas.ini");	// tooltiptext of button
    // action-listener for button
	btnRead.addActionListener(new java.awt.event.ActionListener() 
		{
		public void actionPerformed(ActionEvent e) 
			{
				btnRead_actionPerformed(e);
			}
		});
	// Write-Button
    btnWrite.setText("Write file mcphas.ini");
	btnWrite.setToolTipText("Writes Data to file mcphas.ini; Resets the runtime-parameters!");
    btnWrite.addActionListener(new java.awt.event.ActionListener() 
		{
		public void actionPerformed(ActionEvent e) 
			{
				btnWrite_actionPerformed(e);
			}
		});
	
	/******************************************************
	 * Checkboxes
	 ******************************************************/
	// Exit-Checkbox
	chkExit.setText("Exit");
	chkExit.setToolTipText("Stop McPhase");
    chkExit.addActionListener(new java.awt.event.ActionListener() 
		{
		public void actionPerformed(ActionEvent e) 
			{
				btnExit_actionPerformed(e);
			}
		});
	
	// Pause-Checkbox
	chkPause.setText("Pause");
	chkPause.setToolTipText("Pause McPhase");
    chkPause.addActionListener(new java.awt.event.ActionListener() 
		{
		public void actionPerformed(ActionEvent e) 
			{
				btnPause_actionPerformed(e);
			}
		});
	
	// Display-All-Checkbox
	chkDisplayAll.setText("Display All");
	chkDisplayAll.setToolTipText("Display all magnetic structures while iterating (low speed!)");
    chkDisplayAll.addActionListener(new java.awt.event.ActionListener() 
		{
		public void actionPerformed(ActionEvent e) 
			{
				btnDisplayAll_actionPerformed(e);
			}
		});
	
	// Log FEVSQ-Checkbox
	chkLogFevsQ.setText("Log FEVSQ");
	chkLogFevsQ.setToolTipText("Create a logfile of Propagation vs Free Energy (needs a lot disc space!)");
    chkLogFevsQ.addActionListener(new java.awt.event.ActionListener() 
		{
		public void actionPerformed(ActionEvent e) 
			{
				btnLogFevsQ_actionPerformed(e);
			}
		});
	
	/**********************************
	 * Statustextfield
	 **********************************/
	// Statustext - Textfield
	InitTextfield(txStatus, "Statustext", 700, 19);
	txStatus.setHorizontalAlignment(JTextField.CENTER);
	txStatus.setBackground(Color.lightGray);
	
	/***********************************
	 * First tab-panel
	 ***********************************/
	// Labels
	InitLabel(labXT, "xT:", 70, 19);
	InitLabel(labXHa, "xHa:", 70, 19);
	InitLabel(labXHb, "xHb:", 70, 19);    
	InitLabel(labXHc, "xHc:", 70, 19);
	InitLabel(labXMin, "xmin:", 70, 19);
	InitLabel(labXMax, "xmax:", 70, 19);
	InitLabel(labXStep, "xstep:", 70, 19);
	InitLabel(labXMax1, "ymax:", 70, 19);
	InitLabel(labXMin1, "ymin:", 70, 19);
	InitLabel(labXT1, "yT:", 70, 19);
	InitLabel(labXHa1, "yHa:", 70, 19);
	InitLabel(labXHb1, "yHb:", 70, 19);
	InitLabel(labXHc1, "yHc:", 70, 19);
	InitLabel(labXStep1, "ystep:", 70, 19);

	InitLabel(labT0, "T0:", 70, 19);
	InitLabel(labHa0, "Ha0:", 70, 19);
	InitLabel(labHb0, "Hb0:", 70, 19);
	InitLabel(labHc0, "Hc0:", 70, 19);
	
	// Textfields	
	InitTextfield(txXT, "Vector in (H-T) space corresponding to phasediagram x axis (xT [K] xHa [T] xHb [T] xHc [T])", 50, 19);
	InitTextfield(txXHa, "Vector in (H-T) space corresponding to phasediagram x axis (xT [K] xHa [T] xHb [T] xHc [T])", 50, 19);
	InitTextfield(txXHb, "Vector in (H-T) space corresponding to phasediagram x axis (xT [K] xHa [T] xHb [T] xHc [T])", 50, 19);
	InitTextfield(txXHc, "Vector in (H-T) space corresponding to phasediagram x axis (xT [K] xHa [T] xHb [T] xHc [T])", 50, 19);
	InitTextfield(txXMin, "Minimum of x axis", 50, 19);
	InitTextfield(txXMax, "Maximum of x axis", 50, 19);
	InitTextfield(txXStep, "Stepwidth of x axis ", 50, 19);
	InitTextfield(txYMin, "Minimum of y axis", 50, 19);
	InitTextfield(txYT, "Vector in (H-T) space corresponding to phasediagram y axis (yT [K] yHa [T] yHb [T] yHc [T])", 50, 19);
	InitTextfield(txYHa, "Vector in (H-T) space corresponding to phasediagram y axis (yT [K] yHa [T] yHb [T] yHc [T])", 50, 19);
	InitTextfield(txYHb, "Vector in (H-T) space corresponding to phasediagram y axis (yT [K] yHa [T] yHb [T] yHc [T])", 50, 19);
	InitTextfield(txYHc, "Vector in (H-T) space corresponding to phasediagram y axis (yT [K] yHa [T] yHb [T] yHc [T])", 50, 19);
	InitTextfield(txYStep, "Stepwidth of y axis", 50, 19);
	InitTextfield(txYMax, "Maximum of y axis", 50, 19);

	InitTextfield(txT0, "Temperature offset", 50, 19);
	InitTextfield(txHa0, "Field  offset along a", 50, 19);
	InitTextfield(txHb0, "Field  offset along b", 50, 19);
	InitTextfield(txHc0, "Field  offset along c", 50, 19);
    
	
	/***********************************
	 * Second tab-panel
	 ***********************************/
	// Labels
	InitLabel(labDeltaH, "deltah:", 70, 19);
	InitLabel(labHMax, "hmax:", 70, 19);
	InitLabel(labHMin, "hmin:", 70, 19);
	InitLabel(labDeltaL, "deltal:", 70, 19);
	InitLabel(labLMax, "lmax:", 70, 19);
	InitLabel(labLMin, "lmin:", 70, 19);
	InitLabel(labDeltaK, "deltak:", 70, 19);
	InitLabel(labKMax, "kmax:", 70, 19);
	InitLabel(labKMin, "kmin:", 70, 19);
	InitLabel(labNofRndTries, "nofrndtries:", 150, 19);
	InitLabel(labMaxQPeriod, "maxqperiod:", 150, 19);
	InitLabel(labmaxnoftestspincf, "maxnoftestspincf:", 150, 19);
	InitLabel(labMaxNofSpins, "MaxNofSpins:", 150, 19);
	// Textfields
	InitTextfield(txDeltaH, "test q vector (qmin qmax deltaq)", 50, 19);
 	InitTextfield(txHMax, "test q vector (qmin qmax deltaq)", 50, 19);
 	InitTextfield(txHMin, "test q vector (qmin qmax deltaq)", 50, 19);
	InitTextfield(txDeltaK, "test q vector (qmin qmax deltaq)", 50, 19);
	InitTextfield(txKMax, "test q vector (qmin qmax deltaq)", 50, 19);
	InitTextfield(txKMin, "test q vector (qmin qmax deltaq)", 50, 19);
	InitTextfield(txDeltaL, "test q vector (qmin qmax deltaq)", 50, 19);
	InitTextfield(txLMax, "test q vector (qmin qmax deltaq)", 50, 19);
	InitTextfield(txLMin, "test q vector (qmin qmax deltaq)", 50, 19);
 	InitTextfield(txNofRndTries, "number of random (Monte Carlo) spin inversions  to try for each initial spinconfiguration", 50, 19);
 	InitTextfield(txMaxQPeriod, "maximal periodicity of spinconfigurations generated by q vectors", 50, 19);
 	InitTextfield(txmaxnoftestspincf, "maximal number of test spinconfigurations to be stored", 50, 19);
 	InitTextfield(txMaxNofSpins, "maximal number of spins in spinconfigurations generated by q vectors", 50, 19);

	/***********************************
	 * Third tab-panel
	 ***********************************/
	// Labels
	InitLabel(labMaxSpinChange, "maxspinchange:", 190, 19);
	//InitLabel(labSmallStep, "smallstep:", 150, 19);
	InitLabel(labBigStep, "bigstep:", 190, 19);
	InitLabel(labMaxStaMf, "maxstamf:", 190, 19);
	InitLabel(labMaxNofMfLoops, "maxnofmfloops:", 190, 19);    
	// Textfields
	InitTextfield(txMaxSpinChange, "Sum_{i=1}^{n} abs(change of <Ji> with respect to initial  configuration) - limit to abandon selfconsistency process", 50, 19);
	//InitTextfield(txSmallStep, "Tooltip", 50, 19);
	InitTextfield(txBigStep, "Mean field - step ratio (=actual step/calculated step)", 50, 19);
	InitTextfield(txMaxStaMf, "Standard deviation - limit to end selfconsistency process", 50, 19);
	InitTextfield(txMaxNofMfLoops, "Maximum number of loops to abandon selfconsistency process", 50, 19);

	/***********************************
	 * Fourth tab-panel
	 ***********************************/
	// Labels
	InitLabel(labNofSpinCorrs, "nofspincorrs:", 150, 19);
	InitLabel(labMaxNofHkls, "maxnofhkls:", 150, 19);
	InitLabel(labMaxQ, "maxQ:", 150, 19);
	// Textfields
	InitTextfield(txNofSpinCorrs, "For thermal expansion and magnetostriction - how many spinspin correlation functions should be calculated", 50, 19);
	InitTextfield(txMaxNofHkls, "For Neutron Diffraction - calculation of mxnofhkl strongest reflections", 50, 19);
	InitTextfield(txMaxQ, "For Neutron Diffraction - maximum scattering vector |Q|[1/A] for calculated hkl's", 50, 19);
	 	
	this.add(panMain, BorderLayout.NORTH);
		panMain.setPreferredSize(new Dimension(700,300));
		panMain.addTab("XY Phasediagram Parameters", panXYPhaseParams);
		    panXYPhaseParams.setLayout(flowLayout1);
			panXYPhaseParams.add(jPanel1, null);
				//jPanel1.setPreferredSize(new Dimension(50, 173));
				//jPanel1.setMinimumSize(new Dimension(50, 173));
			    jPanel1.setLayout(verticalFlowLayout1);
					jPanel1.add(labXT, null);
					jPanel1.add(labXHa, null);
					jPanel1.add(labXHb, null);
					jPanel1.add(labXHc, null);
					jPanel1.add(labXMin, null);
					jPanel1.add(labXMax, null);
					jPanel1.add(labXStep, null);
			panXYPhaseParams.add(jPanel2, null);
			    jPanel2.setLayout(verticalFlowLayout2);
					jPanel2.add(txXT, null);
					jPanel2.add(txXHa, null);
					jPanel2.add(txXHb, null);
					jPanel2.add(txXHc, null);
					jPanel2.add(txXMin, null);
					jPanel2.add(txXMax, null);
					jPanel2.add(txXStep, null);
			panXYPhaseParams.add(jPanel3, null);
				//jPanel3.setPreferredSize(new Dimension(50, 173));
				//jPanel3.setMinimumSize(new Dimension(50, 173));
				jPanel3.setLayout(verticalFlowLayout3);
					jPanel3.add(labXT1, null);
					jPanel3.add(labXHa1, null);
					jPanel3.add(labXHb1, null);
					jPanel3.add(labXHc1, null);
					jPanel3.add(labXMin1, null);
					jPanel3.add(labXMax1, null);
					jPanel3.add(labXStep1, null);
			panXYPhaseParams.add(jPanel4, null);
			    jPanel4.setLayout(verticalFlowLayout4);
					jPanel4.add(txYT, null);
					jPanel4.add(txYHa, null);
					jPanel4.add(txYHb, null);
					jPanel4.add(txYHc, null);
					jPanel4.add(txYMin, null);
					jPanel4.add(txYMax, null);
					jPanel4.add(txYStep, null);
			panXYPhaseParams.add(jPanel4a, null);
			    jPanel4a.setLayout(verticalFlowLayout4a);
					jPanel4a.add(labT0, null);
					jPanel4a.add(labHa0, null);
					jPanel4a.add(labHb0, null);
					jPanel4a.add(labHc0, null);
                        panXYPhaseParams.add(jPanel4b, null);
			    jPanel4b.setLayout(verticalFlowLayout4b);
					jPanel4b.add(txT0, null);
					jPanel4b.add(txHa0, null);
					jPanel4b.add(txHb0, null);
					jPanel4b.add(txHc0, null);

		panMain.addTab("Generation of Spin-Configurations", jPanel13);
			jPanel13.setLayout(verticalFlowLayout13);
			jPanel13.add(panQVector,null);
				panQVector.setLayout(flowLayout2);
				panQVector.add(jPanel5, null);
					jPanel5.setLayout(verticalFlowLayout5);
					//jPanel5.setPreferredSize(new Dimension(70, 77));
					//jPanel5.setMinimumSize(new Dimension(70, 77));
						jPanel5.add(labHMin, null);
						jPanel5.add(labHMax, null);
						jPanel5.add(labDeltaH, null);
				panQVector.add(jPanel8, null);
					jPanel8.setLayout(verticalFlowLayout8);
					//jPanel8.setPreferredSize(new Dimension(70, 77));
					//jPanel8.setMinimumSize(new Dimension(70, 77));
						jPanel8.add(txHMin, null);
						jPanel8.add(txHMax, null);
						jPanel8.add(txDeltaH, null);
				panQVector.add(jPanel7, null);
				    jPanel7.setLayout(verticalFlowLayout7);
					//jPanel7.setPreferredSize(new Dimension(70, 77));
					//jPanel7.setMinimumSize(new Dimension(70, 77));
						jPanel7.add(labKMin, null);
						jPanel7.add(labKMax, null);
						jPanel7.add(labDeltaK, null);
				panQVector.add(jPanel9, null);
					jPanel9.setLayout(verticalFlowLayout9);
					//jPanel9.setPreferredSize(new Dimension(70, 77));
					//jPanel9.setMinimumSize(new Dimension(70, 77));
						jPanel9.add(txKMin, null);
						jPanel9.add(txKMax, null);
						jPanel9.add(txDeltaK, null);
				panQVector.add(jPanel6, null);
				    jPanel6.setLayout(verticalFlowLayout6);
					//jPanel6.setPreferredSize(new Dimension(70, 77));
					//jPanel6.setMinimumSize(new Dimension(70, 77));
						jPanel6.add(labLMin, null);
						jPanel6.add(labLMax, null);
						jPanel6.add(labDeltaL, null);
				panQVector.add(jPanel10, null);
					jPanel10.setLayout(verticalFlowLayout10);
					//jPanel10.setPreferredSize(new Dimension(70, 77));
					//jPanel10.setMinimumSize(new Dimension(70, 77));			
						jPanel10.add(txLMin, null);
						jPanel10.add(txLMax, null);
						jPanel10.add(txDeltaL, null);
			jPanel13.add(jPanel14, null);
				jPanel14.setLayout(flowLayout4);
				jPanel14.add(jPanel15, null);
					jPanel15.setLayout(verticalFlowLayout14);
					//jPanel15.setPreferredSize(new Dimension(120, 77));
					//jPanel15.setMinimumSize(new Dimension(120, 77));			
						jPanel15.add(labNofRndTries, null);
						jPanel15.add(labMaxQPeriod, null);
						jPanel15.add(labmaxnoftestspincf, null);
						jPanel15.add(labMaxNofSpins, null);
				jPanel14.add(jPanel16, null);
					jPanel16.setLayout(verticalFlowLayout15);
					//jPanel16.setPreferredSize(new Dimension(70, 77));
					//jPanel16.setMinimumSize(new Dimension(70, 77));			
						jPanel16.add(txNofRndTries, null);
						jPanel16.add(txMaxQPeriod, null);	
						jPanel16.add(txmaxnoftestspincf, null);	
						jPanel16.add(txMaxNofSpins, null);
					
			panMain.addTab("Parameters for sub fecalc - Selfconsistency process", panAdditionalParams);
		    panAdditionalParams.setLayout(flowLayout3);
			panAdditionalParams.add(jPanel12, null);
				//jPanel12.setPreferredSize(new Dimension(140, 240));
				//jPanel12.setMinimumSize(new Dimension(140, 240));
				jPanel12.setLayout(verticalFlowLayout11);
					jPanel12.add(labMaxNofMfLoops, null);
					jPanel12.add(labMaxStaMf, null);
					jPanel12.add(labBigStep, null);
					//jPanel12.add(labSmallStep, null);
					jPanel12.add(labMaxSpinChange, null);
			panAdditionalParams.add(jPanel11, null);
				jPanel11.setLayout(verticalFlowLayout12);
				//jPanel11.setPreferredSize(new Dimension(200, 240));
				//jPanel11.setMinimumSize(new Dimension(200, 240));
					jPanel11.add(txMaxNofMfLoops, null);
					jPanel11.add(txMaxStaMf, null);
					jPanel11.add(txBigStep, null);
					//jPanel11.add(txSmallStep, null);
					jPanel11.add(txMaxSpinChange, null);

			panMain.addTab("Output of physical properties", panOutputOfProperties);
		    panOutputOfProperties.setLayout(flowLayout5);
			panOutputOfProperties.add(jPanel17, null);
				//jPanel17.setPreferredSize(new Dimension(140, 240));
				//jPanel17.setMinimumSize(new Dimension(140, 240));
				jPanel17.setLayout(verticalFlowLayout16);
					jPanel17.add(labNofSpinCorrs, null);
					jPanel17.add(labMaxNofHkls, null);
					jPanel17.add(labMaxQ, null);
			panOutputOfProperties.add(jPanel18, null);
				jPanel18.setLayout(verticalFlowLayout17);
				//jPanel18.setPreferredSize(new Dimension(200, 240));
				//jPanel18.setMinimumSize(new Dimension(200, 240));
					jPanel18.add(txNofSpinCorrs, null);
					jPanel18.add(txMaxNofHkls, null);
					jPanel18.add(txMaxQ, null);
			
	this.add(txStatus, BorderLayout.CENTER);
	
    this.add(panOkCancel, BorderLayout.SOUTH);
		panOkCancel.add(chkExit, null);
		panOkCancel.add(chkPause, null);
		panOkCancel.add(chkDisplayAll, null);
		panOkCancel.add(chkLogFevsQ, null);
		panOkCancel.add(btnWrite, null);
		panOkCancel.add(btnRead, null);

	// Get system variables and save as members
	strWorkingDir = System.getProperty("user.dir");
	strFileSeparator = System.getProperty("file.separator");
	if (!strWorkingDir.endsWith(strFileSeparator))
	{
		strWorkingDir = strWorkingDir + strFileSeparator;
	}
	// Add Filename to working-directory
	strWorkingDir = strWorkingDir + "mcphas.ini";

  }
//Start the applet

  public void start() 
  {
  }
//Stop the applet

  public void stop() 
  {
  }
//Destroy the applet

  public void destroy() 
  {
  }
//Get Applet information

  public String getAppletInfo() {
    return "Applet Information";
  }
//Get parameter info

  public String[][] getParameterInfo() {
    return null;
  }

  private void InitLabel(JLabel jLab, String strText, int iLength, int iHeight)
  {
    jLab.setText(strText);
    jLab.setMaximumSize(new Dimension(iLength, iHeight));
    jLab.setPreferredSize(new Dimension(iLength, iHeight));
    jLab.setMinimumSize(new Dimension(iLength, iHeight));
	//jLab.setHorizontalAlignment(JLabel.RIGHT);
  }

  private void InitTextfield(JTextField jText, String strToolTipText, int iLength, int iHeight)
  {
	jText.setToolTipText(strToolTipText);
  	jText.setPreferredSize(new Dimension(iLength, iHeight));
    jText.setMinimumSize(new Dimension(iLength, iHeight));
	jText.setMaximumSize(new Dimension(iLength, iHeight));
    jText.setHorizontalAlignment(JTextField.RIGHT);
  }
  
  /************************************************
   * Main method
   ************************************************/
  
  public static void main(String[] args) {

	  try
	  {		
		IniConfigurator panel = new IniConfigurator();
				
		JFrame frame = new JFrame("IniConfigurator");
		frame.addWindowListener(new WindowAdapter() 
			{
				public void windowClosing(WindowEvent e) 
				{
					System.exit(0);
				}
			});
		
		frame.getContentPane().add("Center", panel);
		frame.pack();
		frame.setSize(800,400);
		// Initialize Form
		panel.jbInit();
		// Start Reading of INI
		panel.ReadFromIni(panel.strWorkingDir);
		// Set the form visible
		frame.setVisible(true);
	  }
	  catch(Exception e)
	  {
		System.out.println(e.getMessage());
	  }
  }

  void btnRead_actionPerformed(ActionEvent e) 
  {
	  ReadFromIni(strWorkingDir);
  }

  void ReadFromIni(String strFileName)
  {	
	m_IniFile.SetFileName(strFileName);
	DisplayStatus("Reading INI-File !");
    int iRet = m_IniFile.Read();

	if (iRet != 0)
    {
      // An error has occurred !!  (MsgBox ??)
      return;
    }

	String strPause = m_IniFile.GetValue(CONST_MCPHASE_RUNTIME_CONTROL, "pause");
	String strDisplayAll = m_IniFile.GetValue(CONST_MCPHASE_RUNTIME_CONTROL, "displayall");
	String strLogFevQS = m_IniFile.GetValue(CONST_MCPHASE_RUNTIME_CONTROL, "logfevsQ");

	chkPause.setSelected(strPause.equalsIgnoreCase("1"));	
	chkDisplayAll.setSelected(strDisplayAll.equalsIgnoreCase("1"));
	chkLogFevsQ.setSelected(strLogFevQS.equalsIgnoreCase("1"));
	
    txXT.setText(m_IniFile.GetValue(CONST_XY_PHASEDIAGRAM_PARAMETERS, "xT"));
    txXHa.setText(m_IniFile.GetValue(CONST_XY_PHASEDIAGRAM_PARAMETERS, "xHa"));
    txXHb.setText(m_IniFile.GetValue(CONST_XY_PHASEDIAGRAM_PARAMETERS, "xHb"));
    txXHc.setText(m_IniFile.GetValue(CONST_XY_PHASEDIAGRAM_PARAMETERS, "xHc"));
    txXMin.setText(m_IniFile.GetValue(CONST_XY_PHASEDIAGRAM_PARAMETERS, "xmin"));
    txXMax.setText(m_IniFile.GetValue(CONST_XY_PHASEDIAGRAM_PARAMETERS, "xmax"));
    txXStep.setText(m_IniFile.GetValue(CONST_XY_PHASEDIAGRAM_PARAMETERS, "xstep"));
    txYT.setText(m_IniFile.GetValue(CONST_XY_PHASEDIAGRAM_PARAMETERS, "yT"));
    txYHa.setText(m_IniFile.GetValue(CONST_XY_PHASEDIAGRAM_PARAMETERS, "yHa"));
    txYHb.setText(m_IniFile.GetValue(CONST_XY_PHASEDIAGRAM_PARAMETERS, "yHb"));
    txYHc.setText(m_IniFile.GetValue(CONST_XY_PHASEDIAGRAM_PARAMETERS, "yHc"));
    txYMin.setText(m_IniFile.GetValue(CONST_XY_PHASEDIAGRAM_PARAMETERS, "ymin"));
    txYMax.setText(m_IniFile.GetValue(CONST_XY_PHASEDIAGRAM_PARAMETERS, "ymax"));
    txYStep.setText(m_IniFile.GetValue(CONST_XY_PHASEDIAGRAM_PARAMETERS, "ystep"));

    txT0.setText(m_IniFile.GetValue(CONST_XY_PHASEDIAGRAM_PARAMETERS, "T0"));
    txHa0.setText(m_IniFile.GetValue(CONST_XY_PHASEDIAGRAM_PARAMETERS, "Ha0"));
    txHb0.setText(m_IniFile.GetValue(CONST_XY_PHASEDIAGRAM_PARAMETERS, "Hb0"));
    txHc0.setText(m_IniFile.GetValue(CONST_XY_PHASEDIAGRAM_PARAMETERS, "Hc0"));


    txHMin.setText(m_IniFile.GetValue(CONST_GENERATION_OF_SPIN_CONFIGURATIONS, "hmin"));
    txHMax.setText(m_IniFile.GetValue(CONST_GENERATION_OF_SPIN_CONFIGURATIONS, "hmax"));
    txDeltaH.setText(m_IniFile.GetValue(CONST_GENERATION_OF_SPIN_CONFIGURATIONS, "deltah"));
    txKMin.setText(m_IniFile.GetValue(CONST_GENERATION_OF_SPIN_CONFIGURATIONS, "kmin"));
    txKMax.setText(m_IniFile.GetValue(CONST_GENERATION_OF_SPIN_CONFIGURATIONS, "kmax"));
    txDeltaK.setText(m_IniFile.GetValue(CONST_GENERATION_OF_SPIN_CONFIGURATIONS, "deltak"));
    txLMin.setText(m_IniFile.GetValue(CONST_GENERATION_OF_SPIN_CONFIGURATIONS, "lmin"));
    txLMax.setText(m_IniFile.GetValue(CONST_GENERATION_OF_SPIN_CONFIGURATIONS, "lmax"));
    txDeltaL.setText(m_IniFile.GetValue(CONST_GENERATION_OF_SPIN_CONFIGURATIONS, "deltal"));
    txNofRndTries.setText(m_IniFile.GetValue(CONST_GENERATION_OF_SPIN_CONFIGURATIONS, "nofrndtries"));
    txMaxQPeriod.setText(m_IniFile.GetValue(CONST_GENERATION_OF_SPIN_CONFIGURATIONS, "maxqperiod"));	  
    txmaxnoftestspincf.setText(m_IniFile.GetValue(CONST_GENERATION_OF_SPIN_CONFIGURATIONS, "maxnoftestspincf"));	  
	txMaxNofSpins.setText(m_IniFile.GetValue(CONST_GENERATION_OF_SPIN_CONFIGURATIONS, "maxnofspins"));

	txMaxNofMfLoops.setText(m_IniFile.GetValue(CONST_PARAMETERS_FOR_SUB_FECALC_SELFCONSISTENCY_PROCESS, "maxnofmfloops"));
    txMaxStaMf.setText(m_IniFile.GetValue(CONST_PARAMETERS_FOR_SUB_FECALC_SELFCONSISTENCY_PROCESS, "maxstamf"));
    txBigStep.setText(m_IniFile.GetValue(CONST_PARAMETERS_FOR_SUB_FECALC_SELFCONSISTENCY_PROCESS, "bigstep"));
    //txSmallStep.setText(m_IniFile.GetValue(CONST_PARAMETERS_FOR_SUB_FECALC_SELFCONSISTENCY_PROCESS, "smallstep"));
    txMaxSpinChange.setText(m_IniFile.GetValue(CONST_PARAMETERS_FOR_SUB_FECALC_SELFCONSISTENCY_PROCESS, "maxspinchange"));

	txNofSpinCorrs.setText(m_IniFile.GetValue(CONST_OUTPUT_OF_PHYSICAL_PROPERTIES, "nofspincorrs"));
    txMaxNofHkls.setText(m_IniFile.GetValue(CONST_OUTPUT_OF_PHYSICAL_PROPERTIES, "maxnofhkls"));
    txMaxQ.setText(m_IniFile.GetValue(CONST_OUTPUT_OF_PHYSICAL_PROPERTIES, "maxQ"));
  }

  void WriteCommand()
  {	
	String strExit;
	String strPause;
	String strDispAll;
	String strLogFevsQ;
	
	if (chkExit.isSelected() == true)
	{
		strExit = "1";
	}
	else
	{
		strExit = "0";		
	}
	
	if (chkPause.isSelected() == true)
	{
		strPause = "1";
	}
	else
	{
		strPause = "0";
	}
	
	if (chkDisplayAll.isSelected() == true)
	{
		strDispAll = "1";
	}
	else
	{
		strDispAll = "0";
	}
	
	if (chkLogFevsQ.isSelected() == true)
	{
		strLogFevsQ = "1";
	}
	else
	{
		strLogFevsQ = "0";
	}
	m_IniFile.SetValue(CONST_MCPHASE_RUNTIME_CONTROL, "exit", strExit, "to stop program set exit to 1");
	m_IniFile.SetValue(CONST_MCPHASE_RUNTIME_CONTROL, "pause", strPause, "to hold program set pause to 1");
	m_IniFile.SetValue(CONST_MCPHASE_RUNTIME_CONTROL, "displayall", strDispAll, "to display all structures while iterating set displayall to 1\n# (mind that by using this option mcphas gets very slow)");
	m_IniFile.SetValue(CONST_MCPHASE_RUNTIME_CONTROL, "logfevsQ", strLogFevsQ, "to create a logfile of the propagation versus free energy  set logfevsQ to 1\n# (mind this uses a lot of disc space)");
	WriteToIni(strWorkingDir);	  
	DisplayStatus("Sending Runtime-Control to INI-File");
  }
  
  void btnExit_actionPerformed(ActionEvent e) {
	  WriteCommand();
  }							
  
  void btnPause_actionPerformed(ActionEvent e) {
	  WriteCommand();
  }
  
  void btnDisplayAll_actionPerformed(ActionEvent e) {
	  WriteCommand();
  }
  
  void btnLogFevsQ_actionPerformed(ActionEvent e) {
	  WriteCommand();
  }
  
  void btnWrite_actionPerformed(ActionEvent e) {
    //System.out.println("Write !!!");
	m_IniFile.SetValue(CONST_MCPHASE_RUNTIME_CONTROL, "exit", "0", "to stop program set exit to 1");
	m_IniFile.SetValue(CONST_MCPHASE_RUNTIME_CONTROL, "pause", "0", "to hold program set pause to 1");
	m_IniFile.SetValue(CONST_MCPHASE_RUNTIME_CONTROL, "displayall", "0", "to display all structures while iterating set displayall to 1\n# (mind that by using this option mcphas gets very slow)");
	m_IniFile.SetValue(CONST_MCPHASE_RUNTIME_CONTROL, "logfevsQ", "0", "to create a logfile of the propagation versus free energy  set logfevsQ to 1\n# (mind this uses a lot of disc space)");

	WriteToIni(strWorkingDir);
  }
  
  
  void WriteToIni(String strFileName)
  {	  
	m_IniFile.SetFileName(strFileName);
	DisplayStatus("Updating INI-File");
	
	m_IniFile.SetValue(CONST_XY_PHASEDIAGRAM_PARAMETERS, "xT", txXT.getText(), "xy phasediagram axes - parameters\n# structures are calculated in the xy - phasediagram\n# the direction of x and y can be chosen:\n# vector in (H-T) space corresponding to x axis (xT [K] xHa [T] xHb [T] xHc [T])");
    m_IniFile.SetValue(CONST_XY_PHASEDIAGRAM_PARAMETERS, "xHa", txXHa.getText());
    m_IniFile.SetValue(CONST_XY_PHASEDIAGRAM_PARAMETERS, "xHb", txXHb.getText());
    m_IniFile.SetValue(CONST_XY_PHASEDIAGRAM_PARAMETERS, "xHc", txXHc.getText());
    m_IniFile.SetValue(CONST_XY_PHASEDIAGRAM_PARAMETERS, "xmin", txXMin.getText(), "range of x");
    m_IniFile.SetValue(CONST_XY_PHASEDIAGRAM_PARAMETERS, "xmax", txXMax.getText());
    m_IniFile.SetValue(CONST_XY_PHASEDIAGRAM_PARAMETERS, "xstep", txXStep.getText());
    m_IniFile.SetValue(CONST_XY_PHASEDIAGRAM_PARAMETERS, "yT", txYT.getText(), "vector in (H-T) space corresponding to y axis (yT [K] yHa [T] yHb [T] yHc [T])");
    m_IniFile.SetValue(CONST_XY_PHASEDIAGRAM_PARAMETERS, "yHa", txYHa.getText());
    m_IniFile.SetValue(CONST_XY_PHASEDIAGRAM_PARAMETERS, "yHb", txYHb.getText());
    m_IniFile.SetValue(CONST_XY_PHASEDIAGRAM_PARAMETERS, "yHc", txYHc.getText());
    m_IniFile.SetValue(CONST_XY_PHASEDIAGRAM_PARAMETERS, "ymin", txYMin.getText(), "range of y");
    m_IniFile.SetValue(CONST_XY_PHASEDIAGRAM_PARAMETERS, "ymax", txYMax.getText());
    m_IniFile.SetValue(CONST_XY_PHASEDIAGRAM_PARAMETERS, "ystep", txYStep.getText());
    m_IniFile.SetValue(CONST_XY_PHASEDIAGRAM_PARAMETERS, "T0", txT0.getText(), "offset");
    m_IniFile.SetValue(CONST_XY_PHASEDIAGRAM_PARAMETERS, "Ha0", txHa0.getText());
    m_IniFile.SetValue(CONST_XY_PHASEDIAGRAM_PARAMETERS, "Hb0", txHb0.getText());
    m_IniFile.SetValue(CONST_XY_PHASEDIAGRAM_PARAMETERS, "Hc0", txHc0.getText());

    m_IniFile.SetValue(CONST_GENERATION_OF_SPIN_CONFIGURATIONS, "hmin", txHMin.getText(), "test q vector (qmin qmax deltaq)");
    m_IniFile.SetValue(CONST_GENERATION_OF_SPIN_CONFIGURATIONS, "hmax", txHMax.getText());
    m_IniFile.SetValue(CONST_GENERATION_OF_SPIN_CONFIGURATIONS, "deltah", txDeltaH.getText());
    m_IniFile.SetValue(CONST_GENERATION_OF_SPIN_CONFIGURATIONS, "kmin", txKMin.getText());
    m_IniFile.SetValue(CONST_GENERATION_OF_SPIN_CONFIGURATIONS, "kmax", txKMax.getText());
    m_IniFile.SetValue(CONST_GENERATION_OF_SPIN_CONFIGURATIONS, "deltak", txDeltaK.getText());
    m_IniFile.SetValue(CONST_GENERATION_OF_SPIN_CONFIGURATIONS, "lmin", txLMin.getText());
    m_IniFile.SetValue(CONST_GENERATION_OF_SPIN_CONFIGURATIONS, "lmax", txLMax.getText());
    m_IniFile.SetValue(CONST_GENERATION_OF_SPIN_CONFIGURATIONS, "deltal", txDeltaL.getText());
    m_IniFile.SetValue(CONST_GENERATION_OF_SPIN_CONFIGURATIONS, "maxqperiod", txMaxQPeriod.getText(), "maximal periodicity of spinconfigurations generated by q vectors");
    m_IniFile.SetValue(CONST_GENERATION_OF_SPIN_CONFIGURATIONS, "maxnoftestspincf", txmaxnoftestspincf.getText(), "maximal number of test spinconfigurations");
    m_IniFile.SetValue(CONST_GENERATION_OF_SPIN_CONFIGURATIONS, "maxnofspins", txMaxNofSpins.getText(), "maximal number of spins in spinconfigurations generated by q vectors");
    m_IniFile.SetValue(CONST_GENERATION_OF_SPIN_CONFIGURATIONS, "nofrndtries", txNofRndTries.getText(), "number of random (Monte Carlo) spin inversions  to try for each initial spinconfiguration");
	
    m_IniFile.SetValue(CONST_PARAMETERS_FOR_SUB_FECALC_SELFCONSISTENCY_PROCESS, "maxnofmfloops", txMaxNofMfLoops.getText(), "maximum number of selfconsistency loops");
    m_IniFile.SetValue(CONST_PARAMETERS_FOR_SUB_FECALC_SELFCONSISTENCY_PROCESS, "maxstamf", txMaxStaMf.getText(), "standard deviation - limit to end selfconsistency process\n# standard deviation is defined by ...sta=sqrt(sum_{i=1}^{n} (newmf-old mf)i^2/n)\n# the meanfield is given by mf=gj mb H [meV] (gj...lande factor, mb... bohr magneton)");
    m_IniFile.SetValue(CONST_PARAMETERS_FOR_SUB_FECALC_SELFCONSISTENCY_PROCESS, "bigstep", txBigStep.getText(), "mean field step ratio (bigstep=actual step/calculated step<1) to perform actually\n# note: if sta increases - then for 10 iterations set step ratio to smallstep=bigstep/n\n# by default n=5. However, if bigstep>1 then n=integervalue(bigstep) and step ratio=bigstep-n");
    //m_IniFile.SetValue(CONST_PARAMETERS_FOR_SUB_FECALC_SELFCONSISTENCY_PROCESS, "smallstep", txSmallStep.getText(), "a small step (=step/calculated step) to perform actually when sta rises");
    m_IniFile.SetValue(CONST_PARAMETERS_FOR_SUB_FECALC_SELFCONSISTENCY_PROCESS, "maxspinchange", txMaxSpinChange.getText(), "sum_{i=1}^{n} abs(actual change of angular momentum <Ji> with respect to\n# initial  configuration) > maxspinchange will  end selfconsistency process");

	m_IniFile.SetValue(CONST_OUTPUT_OF_PHYSICAL_PROPERTIES, "nofspincorrs", txNofSpinCorrs.getText(), "output of physical properties to compare with experiment\n# 1. For thermal expansion and magnetostriction\n#  how many spinspin correlation functions\n#  should be calculated");
    m_IniFile.SetValue(CONST_OUTPUT_OF_PHYSICAL_PROPERTIES, "maxnofhkls", txMaxNofHkls.getText(), " 2. For Neutron Diffraction\n#  calculation of mxnofhkl strongest reflections");
    m_IniFile.SetValue(CONST_OUTPUT_OF_PHYSICAL_PROPERTIES, "maxQ", txMaxQ.getText(), "maximum scattering vector |Q|[1/A] for calculated hkl's");

	m_IniFile.Write();
  }
  
  void DisplayStatus(String strStatus)
  {
	  txStatus.setText(strStatus);
  }
}



