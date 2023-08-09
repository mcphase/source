
//package mcphas;

import java.io.*;
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
public class graphic_parameters_configurator extends JPanel
{

  String CONST_POINTCHARGES = "POINTCHARGES";
  String CONST_UNIT_CELL_VIEW = "UNIT CELL VIEW";
  String CONST_SPINS_MOMENTS = "SPINS MOMENTS";
  String CONST_DENSITY = "DENSITY";
  String CONST_PHONON = "PHONON";	
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
//  JButton btnRead = new JButton();
  JButton btnWrite = new JButton();
  
  // Checkboxes
  JCheckBox chkExit = new JCheckBox();
  
  // textfield for Status-display
 // JTextField txStatus = new JTextField();
  
  // Tab-Panel
  JTabbedPane panMain = new JTabbedPane();
  // 1. Tab
  JPanel panPointcharges = new JPanel();	// [UNIT CELL VIEW]
	// Labels
	JLabel labShow_Pointcharges = new JLabel();
	JLabel labScale_Pointcharges = new JLabel();
	JLabel labShow_Pointcharges1 = new JLabel();
	JLabel labScale_Pointcharges1 = new JLabel();
	
	// textfields
	JTextField txShow_Pointcharges = new JTextField();
	JTextField txScale_Pointcharges = new JTextField();
	JTextField tShow_Pointcharges0 = new JTextField();
	JTextField tScale_Pointcharges0 = new JTextField();
  
  // 2. Tab
  JPanel panQVector = new JPanel();			// UNITCELL VIEW
    // Labels
	JLabel labshow_abc_unitcell = new JLabel();
	JLabel labshow_primitive_crystal_unitcell = new JLabel();
	JLabel labshow_magnetic_unitcell = new JLabel();
	JLabel labshow_atoms = new JLabel();
	JLabel labscale_view_1 = new JLabel();
	JLabel labscale_view_2 = new JLabel();
	JLabel labscale_view_3 = new JLabel();
	JLabel labshowprim = new JLabel();
  	// textfields
	JTextField txshow_abc_unitcell = new JTextField();
	JTextField txshow_primitive_crystal_unitcell = new JTextField();
	JTextField txshow_magnetic_unitcell = new JTextField();
	JTextField txshowprim = new JTextField();
	JTextField txshow_atoms = new JTextField();
	JTextField txscale_view_1 = new JTextField();
	JTextField txscale_view_2 = new JTextField();
	JTextField txscale_view_3 = new JTextField();

  // 3. Tab		
  JPanel panAdditionalParams = new JPanel();	// MOMENTS
	// Labels
	JLabel labspins_scale_moment = new JLabel();
	JLabel labspins_show_oscillation = new JLabel();
	JLabel labspins_wave_amplitude = new JLabel();
	JLabel labspins_show_ellipses = new JLabel();
	JLabel labspins_show_static_moment_direction = new JLabel();
	// textfields
	JTextField txspins_scale_moment = new JTextField();
	JTextField txspins_show_oscillation = new JTextField();
	JTextField txspins_wave_amplitude = new JTextField();
	JTextField txspins_show_ellipses = new JTextField();
	JTextField txspins_show_static_moment_direction = new JTextField();

  // 4. Tab		
  JPanel panOutputOfProperties = new JPanel();	// DENSITY
	// Labels
	JLabel labshow_density = new JLabel();
	JLabel labscale_density_vectors = new JLabel();
	JLabel labdensity_dtheta = new JLabel();
	JLabel labdensity_dfi = new JLabel();
	JLabel labdensity_threshhold = new JLabel();
	JLabel labdensity_gridi = new JLabel();
	JLabel labdensity_gridj = new JLabel();
	JLabel labdensity_gridk = new JLabel();
	// textfields
	JTextField txshow_density = new JTextField();
	JTextField txscale_density_vectors = new JTextField();
	JTextField txdensity_dtheta = new JTextField();
	JTextField txdensity_dfi = new JTextField();
	JTextField txdensity_threshhold = new JTextField();
	JTextField txdensity_gridi = new JTextField();
	JTextField txdensity_gridj = new JTextField();
	JTextField txdensity_gridk = new JTextField();

  // 5. Tab		
  JPanel panPhonon = new JPanel();	// DENSITY
	// Labels
	JLabel labphonon_scale_static_displacements = new JLabel();
	JLabel labphonon_wave_amplitude = new JLabel();
	// textfields
	JTextField txphonon_scale_static_displacements = new JTextField();
	JTextField txphonon_wave_amplitude = new JTextField();


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
  JPanel jPanel11a = new JPanel();
  JPanel jPanel12a = new JPanel();
  JPanel jPanel13 = new JPanel();
  JPanel jPanel14 = new JPanel();
  JPanel jPanel15 = new JPanel();
  JPanel jPanel16 = new JPanel();
  
  JPanel jPanel17 = new JPanel();
  JPanel jPanel18 = new JPanel();
  JPanel jPanel17a = new JPanel();
  JPanel jPanel17b = new JPanel();
  JPanel jPanel18a = new JPanel();
  JPanel jPanel18b = new JPanel();
  JPanel jPanel19 = new JPanel();
  JPanel jPanel20 = new JPanel();
  
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
  VerticalFlowLayout verticalFlowLayout18 = new VerticalFlowLayout();
  FlowLayout flowLayout6 = new FlowLayout();
  VerticalFlowLayout verticalFlowLayout19 = new VerticalFlowLayout();
  VerticalFlowLayout verticalFlowLayout20 = new VerticalFlowLayout();


  /***************************************************************
   * Initializing function
   * All elements of the form are initialized
   ***************************************************************/
  private void jbInit() throws Exception 
  {
   
	/**************************************
	 * Buttons (Read / Write)
	 **************************************/
	
	// Read-Button
//    btnRead.setText("Read file results/graphic_parameters.set");		// display-text of the button
//	btnRead.setToolTipText("Reads Data from file results/graphic_parameters.set");	// tooltiptext of button
    // action-listener for button
/*	btnRead.addActionListener(new java.awt.event.ActionListener()
		{
		public void actionPerformed(ActionEvent e) 
			{
				btnRead_actionPerformed(e);
			}
		});
*/
	// Write-Button
    btnWrite.setText("Write file results/graphic_parameters.set and recalculate");
	btnWrite.setToolTipText("Writes Data to file results/graphic_parameters.set; Resets the runtime-parameters!");
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
	chkExit.setToolTipText("Stop Graphics");
    chkExit.addActionListener(new java.awt.event.ActionListener() 
		{
		public void actionPerformed(ActionEvent e) 
			{
				btnExit_actionPerformed(e);
			}
		});
	
	/**********************************
	 * Statustext
	 **********************************/
	// Statustext - textfield
	//InitTextfield(txStatus, "Statustext", 700, 19);
	//txStatus.setHorizontalAlignment(JTextField.CENTER);
	//txStatus.setBackground(Color.lightGray);
	
	/***********************************
	 * First tab-panel
	 ***********************************/
	// Labels
	InitLabel(labShow_Pointcharges, "Show_Pointcharges (1=yes,0=no):", 400, 19);
	InitLabel(labScale_Pointcharges, "Scale_Pointcharges (common radius scale factor):", 400, 19);

	// textfields
	InitTextfield(txShow_Pointcharges, "Show Pointcharges as Spheres (1=yes,0=no)", 50, 19);
	InitTextfield(txScale_Pointcharges, "Common Scale factor for Pointcharge Sphere Radius", 50, 19);
    
	
	/***********************************
	 * Second tab-panel
	 ***********************************/
	// Labels
	InitLabel(labshow_abc_unitcell, "show_abc_unitcell (1=yes,0=no):", 220, 19);
	InitLabel(labshow_primitive_crystal_unitcell, "show_primitive_crystal_unitcell:", 200, 19);
	InitLabel(labshow_magnetic_unitcell, "show_magnetic_unitcell:", 170, 19);
	InitLabel(labshow_atoms, "show_atoms (1=yes,0=no):", 100, 19);
	InitLabel(labscale_view_1, "scale_view_1:", 100, 19);
	InitLabel(labscale_view_2, "scale_view_2:", 100, 19);
	InitLabel(labscale_view_3, "scale_view_3:", 100, 19);
	InitLabel(labshowprim, "showprim (atoms in prim unit):", 170, 19);
	// textfields
	InitTextfield(txshow_abc_unitcell, "show abc unitcell (1=yes,0=no)", 50, 19);
 	InitTextfield(txshow_primitive_crystal_unitcell, "show prim crystal unitcell (1=yes,0=no)", 50, 19);
 	InitTextfield(txshow_magnetic_unitcell, "show magnetic unit cell (1=yes,0=no)", 50, 19);
	InitTextfield(txshowprim, "show atoms nuclei in primitive unit cell (1=yes,0=no)", 50, 19);
	InitTextfield(txshow_atoms, "show atoms nuclei (1=yes,0=no)", 50, 19);
	InitTextfield(txscale_view_1, "scale view of crystal in 1 direction (default=1)", 50, 19);
	InitTextfield(txscale_view_2, "scale view of crystal in 2 direction (default=1)", 50, 19);
	InitTextfield(txscale_view_3, "scale view of crystal in 3 direction (default=1)", 50, 19);
 	
	/***********************************
	 * Third tab-panel
	 ***********************************/
	// Labels
	InitLabel(labspins_scale_moment, "scale moment:", 190, 19);
	InitLabel(labspins_show_oscillation, "spins_show_oscillation  (1=yes,0=no):", 210, 19);
	InitLabel(labspins_wave_amplitude, "densitywave&spins_wave_amplitude:", 190, 19);
	InitLabel(labspins_show_ellipses, "show ellipses (1=yes,0=no):", 190, 19);
	InitLabel(labspins_show_static_moment_direction, "show static moment direction:", 190, 19);
	// textfields
	InitTextfield(txspins_scale_moment, "scale moments by common scale factor", 50, 19);
	InitTextfield(txspins_show_oscillation, "show oscillation of spins (1=yes,0=no)", 50, 19);
	InitTextfield(txspins_wave_amplitude, "spin&charge density wave amplitude common scale factor", 50, 19);
	InitTextfield(txspins_show_ellipses, "show ellipses indicating precession of the magnetic moment (1=yes,0=no)", 50, 19);
	InitTextfield(txspins_show_static_moment_direction, "show static moment direction (1=yes,0=no)", 50, 19);

	/***********************************
	 * Fourth tab-panel
	 ***********************************/
	// Labels
	InitLabel(labscale_density_vectors, "scale_density_vectors:", 150, 19);
	InitLabel(labshow_density, "show_density (1=yes,0=no):", 150, 19);
	InitLabel(labdensity_dtheta, "dtheta (rad):", 70, 19);
	InitLabel(labdensity_dfi, "dfi (rad):", 70, 19);
	InitLabel(labdensity_threshhold, "threshhold (in |e|/A^3,<0 for electrons):", 300, 19);
	InitLabel(labdensity_gridi, "gridi:", 50, 19);
	InitLabel(labdensity_gridj, "gridj:", 50, 19);
	InitLabel(labdensity_gridk, "gridk:", 50, 19);
	// textfields
	InitTextfield(txscale_density_vectors, "common scale factor for vector density plots e.g. current density", 50, 19);
	InitTextfield(txshow_density, "show the density", 50, 19);
	InitTextfield(txdensity_dtheta, "dtheta for const density surface", 50, 19);
	InitTextfield(txdensity_dfi, "dfi for const density surface", 50, 19);
	InitTextfield(txdensity_threshhold, "threshhold=constant for calculation of constant density surface", 50, 19);
	InitTextfield(txdensity_gridi, "number of gridpoints in i direction for file *.grid", 50, 19);
	InitTextfield(txdensity_gridj, "number of gridpoints in j direction for file *.grid", 50, 19);
	InitTextfield(txdensity_gridk, "number of gridpoints in k direction for file *.grid", 50, 19);

	/***********************************
	 * Fifth tab-panel
	 ***********************************/
	// Labels
	InitLabel(labphonon_scale_static_displacements, "scale_static_displacements:", 170, 19);
	InitLabel(labphonon_wave_amplitude, "phonon_wave_amplitude:", 170, 19);
	// textfields
	InitTextfield(txphonon_scale_static_displacements, "scale factor for static lattice displacements (by mean fields)", 50, 19);
	InitTextfield(txphonon_wave_amplitude, "scale factor for oscillation of lattice positions", 50, 19);


	this.add(panMain, BorderLayout.NORTH);
		panMain.setPreferredSize(new Dimension(790,120));
		panMain.addTab("POINTCHARGES", panPointcharges);
		    panPointcharges.setLayout(flowLayout1);
			panPointcharges.add(jPanel1, null);
				jPanel1.setLayout(verticalFlowLayout1);
					jPanel1.add(labShow_Pointcharges, null);
					jPanel1.add(labScale_Pointcharges, null);
			panPointcharges.add(jPanel2, null);
			    jPanel2.setLayout(verticalFlowLayout2);
					jPanel2.add(txShow_Pointcharges, null);
					jPanel2.add(txScale_Pointcharges, null);
			panPointcharges.add(jPanel3, null);
				jPanel3.setLayout(verticalFlowLayout3);
			panPointcharges.add(jPanel4, null);
			    jPanel4.setLayout(verticalFlowLayout4);
			panPointcharges.add(jPanel4a, null);
			    jPanel4a.setLayout(verticalFlowLayout4a);
                        panPointcharges.add(jPanel4b, null);
			    jPanel4b.setLayout(verticalFlowLayout4b);

		panMain.addTab("UNIT CELL VIEW", jPanel13);
			jPanel13.setLayout(verticalFlowLayout13);
			jPanel13.add(panQVector,null);
				panQVector.setLayout(flowLayout2);
				panQVector.add(jPanel5, null);
					jPanel5.setLayout(verticalFlowLayout5);
						jPanel5.add(labshow_abc_unitcell, null);
						jPanel5.add(labshow_primitive_crystal_unitcell, null);
						jPanel5.add(labshow_atoms, null);
				panQVector.add(jPanel8, null);
					jPanel8.setLayout(verticalFlowLayout8);
						jPanel8.add(txshow_abc_unitcell, null);
						jPanel8.add(txshow_primitive_crystal_unitcell, null);
						jPanel8.add(txshow_atoms, null);
				panQVector.add(jPanel7, null);
				    jPanel7.setLayout(verticalFlowLayout7);
						jPanel7.add(labshow_magnetic_unitcell, null);
						jPanel7.add(labshowprim, null);
				panQVector.add(jPanel9, null);
					jPanel9.setLayout(verticalFlowLayout9);
						jPanel9.add(txshow_magnetic_unitcell, null);
						jPanel9.add(txshowprim, null);
				panQVector.add(jPanel6, null);
				    jPanel6.setLayout(verticalFlowLayout6);
						jPanel6.add(labscale_view_1, null);
						jPanel6.add(labscale_view_2, null);
						jPanel6.add(labscale_view_3, null);
				panQVector.add(jPanel10, null);
					jPanel10.setLayout(verticalFlowLayout10);
						jPanel10.add(txscale_view_1, null);
						jPanel10.add(txscale_view_2, null);
						jPanel10.add(txscale_view_3, null);
			jPanel13.add(jPanel14, null);
				jPanel14.setLayout(flowLayout4);
				jPanel14.add(jPanel15, null);
					jPanel15.setLayout(verticalFlowLayout14);
				jPanel14.add(jPanel16, null);
					jPanel16.setLayout(verticalFlowLayout15);
					
			panMain.addTab("SPINS/MOMENTS", panAdditionalParams);
		    panAdditionalParams.setLayout(flowLayout3);
			panAdditionalParams.add(jPanel12, null);
				jPanel12.setLayout(verticalFlowLayout12);
					jPanel12.add(labspins_scale_moment, null);
					jPanel12.add(labspins_show_oscillation, null);
					jPanel12.add(labspins_wave_amplitude, null);
			panAdditionalParams.add(jPanel11, null);
				jPanel11.setLayout(verticalFlowLayout11);
					jPanel11.add(txspins_scale_moment, null);
					jPanel11.add(txspins_show_oscillation, null);
					jPanel11.add(txspins_wave_amplitude, null);

			panAdditionalParams.add(jPanel12a, null);
				jPanel12a.setLayout(verticalFlowLayout11);
					jPanel12a.add(labspins_show_static_moment_direction, null);
					jPanel12a.add(labspins_show_ellipses, null);
			panAdditionalParams.add(jPanel11a, null);
				jPanel11a.setLayout(verticalFlowLayout12);
					jPanel11a.add(txspins_show_static_moment_direction, null);
                                        jPanel11a.add(txspins_show_ellipses, null);

			panMain.addTab("DENSITY", panOutputOfProperties);
		    panOutputOfProperties.setLayout(flowLayout5);
			panOutputOfProperties.add(jPanel17, null);
				jPanel17.setLayout(verticalFlowLayout17);
					jPanel17.add(labscale_density_vectors, null);
					jPanel17.add(labshow_density, null);
			panOutputOfProperties.add(jPanel18, null);
				jPanel18.setLayout(verticalFlowLayout18);
					jPanel18.add(txscale_density_vectors, null);
					jPanel18.add(txshow_density, null);
			panOutputOfProperties.add(jPanel17a, null);
				jPanel17a.setLayout(verticalFlowLayout17);
					jPanel17a.add(labdensity_dtheta, null);
					jPanel17a.add(labdensity_dfi, null);
					jPanel17a.add(labdensity_threshhold, null);
			panOutputOfProperties.add(jPanel18a, null);
				jPanel18a.setLayout(verticalFlowLayout18);
					jPanel18a.add(txdensity_dtheta, null);
					jPanel18a.add(txdensity_dfi, null);
					jPanel18a.add(txdensity_threshhold, null);
                        panOutputOfProperties.add(jPanel17b, null);
                                jPanel17b.setLayout(verticalFlowLayout17);
					jPanel17b.add(labdensity_gridi, null);
					jPanel17b.add(labdensity_gridj, null);
					jPanel17b.add(labdensity_gridk, null);
			panOutputOfProperties.add(jPanel18b, null);
				jPanel18b.setLayout(verticalFlowLayout18);
					jPanel18b.add(txdensity_gridi, null);
					jPanel18b.add(txdensity_gridj, null);
					jPanel18b.add(txdensity_gridk, null);

			panMain.addTab("PHONON", panPhonon);
		    panPhonon.setLayout(flowLayout6);
			panPhonon.add(jPanel19, null);
				jPanel19.setLayout(verticalFlowLayout19);
					jPanel19.add(labphonon_scale_static_displacements, null);
					jPanel19.add(labphonon_wave_amplitude, null);
			panPhonon.add(jPanel20, null);
				jPanel20.setLayout(verticalFlowLayout20);
					jPanel20.add(txphonon_scale_static_displacements, null);
					jPanel20.add(txphonon_wave_amplitude, null);


			
	//this.add(txStatus, BorderLayout.CENTER);
	
    this.add(panOkCancel, BorderLayout.SOUTH);
		panOkCancel.add(chkExit, null);
		panOkCancel.add(btnWrite, null);
//		panOkCancel.add(btnRead, null);

	// Get system variables and save as members
	strWorkingDir = System.getProperty("user.dir");
	strFileSeparator = System.getProperty("file.separator");
	if (!strWorkingDir.endsWith(strFileSeparator))
	{
		strWorkingDir = strWorkingDir + strFileSeparator;
	}
	// Add Filename to working-directory
	strWorkingDir = strWorkingDir + "results/graphic_parameters.set";

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

  private void InitLabel(JLabel jLab, String strtext, int iLength, int iHeight)
  {
    jLab.setText(strtext);
    jLab.setMaximumSize(new Dimension(iLength, iHeight));
    jLab.setPreferredSize(new Dimension(iLength, iHeight));
    jLab.setMinimumSize(new Dimension(iLength, iHeight));
	//jLab.setHorizontalAlignment(JLabel.RIGHT);
  }

  private void InitTextfield(JTextField jtext, String strToolTiptext, int iLength, int iHeight)
  {
	jtext.setToolTipText(strToolTiptext);
  	jtext.setPreferredSize(new Dimension(iLength, iHeight));
    jtext.setMinimumSize(new Dimension(iLength, iHeight));
	jtext.setMaximumSize(new Dimension(iLength, iHeight));
    jtext.setHorizontalAlignment(JTextField.RIGHT);
  }
  
  /************************************************
   * Main method
   ************************************************/
  
  public static void main(String[] args) {

	  try
	  {		
		graphic_parameters_configurator panel = new graphic_parameters_configurator();
				
		JFrame frame = new JFrame("graphic_parameters_configurator");
		frame.addWindowListener(new WindowAdapter() 
			{
				public void windowClosing(WindowEvent e) 
				{
                                         System.exit(0);
				}
			});
		
		frame.getContentPane().add("Center", panel);
		frame.pack();
		frame.setSize(800,200);
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
//	DisplayStatus("Reading results/graphic_parameters.set-File !");
    int iRet = m_IniFile.Read();

	if (iRet != 0)
    {
      // An error has occurred !!  (MsgBox ??)
      return;
    }


    txShow_Pointcharges.setText(m_IniFile.GetValue(CONST_POINTCHARGES, "show_pointcharges"));
    txScale_Pointcharges.setText(m_IniFile.GetValue(CONST_POINTCHARGES, "scale_pointcharges"));
    
    txshow_magnetic_unitcell.setText(m_IniFile.GetValue(CONST_UNIT_CELL_VIEW, "show_magnetic_unitcell"));
    txshow_primitive_crystal_unitcell.setText(m_IniFile.GetValue(CONST_UNIT_CELL_VIEW, "show_primitive_crystal_unitcell"));
    txshow_abc_unitcell.setText(m_IniFile.GetValue(CONST_UNIT_CELL_VIEW, "show_abc_unitcell"));
    txshowprim.setText(m_IniFile.GetValue(CONST_UNIT_CELL_VIEW, "showprim"));
    txscale_view_3.setText(m_IniFile.GetValue(CONST_UNIT_CELL_VIEW, "scale_view_3"));
    txscale_view_2.setText(m_IniFile.GetValue(CONST_UNIT_CELL_VIEW, "scale_view_2"));
    txscale_view_1.setText(m_IniFile.GetValue(CONST_UNIT_CELL_VIEW, "scale_view_1"));
    txshow_atoms.setText(m_IniFile.GetValue(CONST_UNIT_CELL_VIEW, "show_atoms"));
   
    txspins_show_static_moment_direction.setText(m_IniFile.GetValue(CONST_SPINS_MOMENTS, "spins_show_static_moment_direction"));
    txspins_show_ellipses.setText(m_IniFile.GetValue(CONST_SPINS_MOMENTS, "spins_show_ellipses"));
    txspins_wave_amplitude.setText(m_IniFile.GetValue(CONST_SPINS_MOMENTS, "spins_wave_amplitude"));
    txspins_show_oscillation.setText(m_IniFile.GetValue(CONST_SPINS_MOMENTS, "spins_show_oscillation"));
    txspins_scale_moment.setText(m_IniFile.GetValue(CONST_SPINS_MOMENTS, "spins_scale_moment"));

    txscale_density_vectors.setText(m_IniFile.GetValue(CONST_DENSITY, "scale_density_vectors"));
    txshow_density.setText(m_IniFile.GetValue(CONST_DENSITY, "show_density"));
    txdensity_dtheta.setText(m_IniFile.GetValue(CONST_DENSITY, "density_dtheta"));
    txdensity_dfi.setText(m_IniFile.GetValue(CONST_DENSITY, "density_dfi"));
    txdensity_threshhold.setText(m_IniFile.GetValue(CONST_DENSITY, "density_threshhold"));
    txdensity_gridi.setText(m_IniFile.GetValue(CONST_DENSITY, "gridi"));
    txdensity_gridj.setText(m_IniFile.GetValue(CONST_DENSITY, "gridj"));
    txdensity_gridk.setText(m_IniFile.GetValue(CONST_DENSITY, "gridk"));

    txphonon_wave_amplitude.setText(m_IniFile.GetValue(CONST_PHONON, "phonon_wave_amplitude"));
    txphonon_scale_static_displacements.setText(m_IniFile.GetValue(CONST_PHONON, "phonon_scale_static_displacements"));

  }

  void WriteCommand()
  {	
	String strExit;
	
	if (chkExit.isSelected() == true)
	{
		strExit = "1";
	}
	else
	{
		strExit = "0";		
	}
	
//	m_IniFile.SetValue(CONST_POINTCHARGES, "exit", strExit, "to stop program set exit to 1");
	WriteToIni(strWorkingDir);	  
	//DisplayStatus("Sending Runtime-Control to INI-File");
  }
  
  void btnExit_actionPerformed(ActionEvent e) {
//	  WriteCommand();
  File f1 = new File(strWorkingDir);
  f1.delete();
  System.exit(0);
  }
							
  
  
  void btnWrite_actionPerformed(ActionEvent e) {
    //System.out.println("Write !!!");
//	m_IniFile.SetValue(CONST_POINTCHARGES, "exit", "0", "to stop program set exit to 1");

	WriteToIni(strWorkingDir);
  }
  
  
  void WriteToIni(String strFileName)
  {	  
	m_IniFile.SetFileName(strFileName);
	//DisplayStatus("Updating graphic parameters -File");
//System.out.println("Writing"+txShow_Pointcharges.getText());
	
    m_IniFile.SetValue(CONST_POINTCHARGES, "show_pointcharges", txShow_Pointcharges.getText(), " ");
    m_IniFile.SetValue(CONST_POINTCHARGES, "scale_pointcharges", txScale_Pointcharges.getText());
  
    m_IniFile.SetValue(CONST_UNIT_CELL_VIEW, "show_magnetic_unitcell", txshow_magnetic_unitcell.getText(), " ");
    m_IniFile.SetValue(CONST_UNIT_CELL_VIEW, "show_primitive_crystal_unitcell", txshow_primitive_crystal_unitcell.getText());
    m_IniFile.SetValue(CONST_UNIT_CELL_VIEW, "show_abc_unitcell", txshow_abc_unitcell.getText());
    m_IniFile.SetValue(CONST_UNIT_CELL_VIEW, "showprim", txshowprim.getText());
    m_IniFile.SetValue(CONST_UNIT_CELL_VIEW, "scale_view_3", txscale_view_3.getText());
    m_IniFile.SetValue(CONST_UNIT_CELL_VIEW, "scale_view_2", txscale_view_2.getText());
    m_IniFile.SetValue(CONST_UNIT_CELL_VIEW, "scale_view_1", txscale_view_1.getText());
    m_IniFile.SetValue(CONST_UNIT_CELL_VIEW, "show_atoms", txshow_atoms.getText());
    	
    m_IniFile.SetValue(CONST_SPINS_MOMENTS, "spins_show_static_moment_direction", txspins_show_static_moment_direction.getText(), " ");
    m_IniFile.SetValue(CONST_SPINS_MOMENTS, "spins_show_ellipses", txspins_show_ellipses.getText(), " ");
    m_IniFile.SetValue(CONST_SPINS_MOMENTS, "spins_wave_amplitude", txspins_wave_amplitude.getText(), " ");
    m_IniFile.SetValue(CONST_SPINS_MOMENTS, "spins_show_oscillation", txspins_show_oscillation.getText(), " ");
    m_IniFile.SetValue(CONST_SPINS_MOMENTS, "spins_scale_moment", txspins_scale_moment.getText(), " ");

    m_IniFile.SetValue(CONST_DENSITY, "scale_density_vectors", txscale_density_vectors.getText(), " ");
    m_IniFile.SetValue(CONST_DENSITY, "show_density", txshow_density.getText(), " ");
    m_IniFile.SetValue(CONST_DENSITY, "density_dtheta", txdensity_dtheta.getText(), " ");
    m_IniFile.SetValue(CONST_DENSITY, "density_dfi", txdensity_dfi.getText(), " ");
    m_IniFile.SetValue(CONST_DENSITY, "density_threshhold", txdensity_threshhold.getText(), " ");
    m_IniFile.SetValue(CONST_DENSITY, "gridi", txdensity_gridi.getText(), " ");
    m_IniFile.SetValue(CONST_DENSITY, "gridj", txdensity_gridj.getText(), " ");
    m_IniFile.SetValue(CONST_DENSITY, "gridk", txdensity_gridk.getText(), " ");

    m_IniFile.SetValue(CONST_PHONON, "phonon_wave_amplitude", txphonon_wave_amplitude.getText(), " ");
    m_IniFile.SetValue(CONST_PHONON, "phonon_scale_static_displacements", txphonon_scale_static_displacements.getText(), " ");

	m_IniFile.Write();
        System.exit(0);
  }
  
  //void DisplayStatus(String strStatus)
 // {
//	  txStatus.setText(strStatus);
 // }
}



