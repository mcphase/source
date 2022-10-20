
//Title:        McPhas INI-Configurator
//Version:
//Copyright:    Copyright (c) 1998
//Author:       Stefan Rotter
//Company:      Stefan Rotter
//Description:  Configurator for mcphas INI-File


//package mcphas;
import java.util.Hashtable;
import java.util.Enumeration;
import java.io.*;

public class IniFile {

private String m_strFileName;
private Hashtable m_Keys;
private int iMaxKeys;

  public IniFile() {
    m_Keys = new Hashtable();
    m_strFileName = "";
  }

  public IniFile(String strFileName) {
    m_Keys = new Hashtable();
    m_strFileName = strFileName;
  }

  public void SetFileName(String strFileName) {
    m_strFileName = strFileName;
  }

  public String GetFileName() {
    return(m_strFileName);
  }

  public int Read() {
    File fileIni;
    int iRet;
    String strLine;
    String strKey = "";
    String strItem = "";
    String strValue = "";

	int iSequKey = 0;
	int iSequItem = 0;
	
    //System.out.println("Read File: " + m_strFileName);
    try{
    fileIni = new File(m_strFileName);

    //Oeffnen der Ini-Datei
    DataInputStream inStream = new DataInputStream(new FileInputStream(fileIni));


     //Auslesen der Ini-Datei
    while (inStream.available() > 0)
    {
      strLine = inStream.readLine();
      if ((strLine.length() == 0)
        ||(strLine.substring(0, 1).equalsIgnoreCase("#")))
      {
        continue;
      }

      if (IsKey(strLine))
      {
        strKey = GetKeyFromLine(strLine);
		iSequKey++;
		iSequItem = 0; // now we start a new Key-session
        //System.out.println("New Key: " + strKey);
      }

      IniKey m_Key = (IniKey)m_Keys.get(strKey);

      if (m_Key == null)
      {
        m_Key = new IniKey(strKey);
		m_Key.iSequence = iSequKey;
        m_Keys.put(strKey, m_Key);
      }

      int iPos = strLine.indexOf("=");
      if (iPos < 0)
      {
        continue;
      }

	  iSequItem++;
	  
      strItem = TrimString(strLine.substring(0, iPos));
      strValue = TrimString(strLine.substring(iPos + 1, strLine.length()));

      IniItem m_Item = new IniItem(strItem, strValue, 1);

	  m_Item.iSequence = iSequItem;
	  
      int i = m_Key.AddItem(m_Item);

      //System.out.println("LINE: '" + strLine + "'");
      //System.out.println("Item: '" + strItem + "', Value: '" + strValue + "'");

    }
	iMaxKeys = iSequKey;
    inStream.close();
    iRet = 0;
    }
    catch(EOFException e)
    {
      System.out.println("EOF: " + e.getLocalizedMessage());
      iRet = 1;
      //EntSession.CWatch("Unplanmaessiges 'End Of File' in DSN-Konfigurationsdatei!");
    }

    catch (FileNotFoundException e)
    {
      System.out.println("File not found: " + m_strFileName + ", " + e.getLocalizedMessage());
      iRet = 2;
      //EntSession.CWatch("Konfigurationsdatei cti_listener.ini nicht gefunden!");
    }

    //Sonstiger Dateifehler
    catch (IOException e)
    {
      System.out.println("Dateifehler: " + e.getLocalizedMessage());
      iRet = 3;
      //EntSession.CWatch("Fehler beim Zugriff auf Datei cti_listener.ini!");
    }

    
    return(iRet);

  }

  public int Write() {
    File fileIni;
    int iRet = 0;
    String strLine;
    String strKey = "";
    String strItem = "";
    String strValue = "";

	
    //System.out.println("Read File: " + m_strFileName);
    try
	{
		fileIni = new File(m_strFileName);
		DataOutputStream out_Stream = new DataOutputStream(new FileOutputStream(fileIni));
	
					out_Stream.writeBytes("#\n#<!--mcphase.mcphas.ini-->\n");
		
		for (int iSequKey = 1; iSequKey <= iMaxKeys;++iSequKey)
		{
		
			IniKey m_Key = GetNextKey(iSequKey);
			if (m_Key != null)
			{
				if (m_Key.GetKey().length() > 0)
				{
					out_Stream.writeBytes("[" + m_Key.GetKey() + "]\n");
					out_Stream.writeBytes(m_Key.GetItemsForWrite()+ "\n");		
				}
			}
			else
			{
				continue;
			}
		}
		out_Stream.close();
		
	}
    catch(EOFException e)
    {
      System.out.println("EOF: " + e.getLocalizedMessage());
      iRet = 1;
      //EntSession.CWatch("Unplanm��iges 'End Of File' in DSN-Konfigurationsdatei!");
    }

    catch (FileNotFoundException e)
    {
      System.out.println("File not found: " + m_strFileName + ", " + e.getLocalizedMessage());
      iRet = 2;
      //EntSession.CWatch("Konfigurationsdatei cti_listener.ini nicht gefunden!");
    }

    //Sonstiger Dateifehler
    catch (IOException e)
    {
      System.out.println("Dateifehler: " + e.getLocalizedMessage());
      iRet = 3;
      //EntSession.CWatch("Fehler beim Zugriff auf Datei cti_listener.ini!");
    }
	return(iRet);
	
  }

  public String GetValue(String strKey, String strItem)
  {
    IniKey m_Key = (IniKey)m_Keys.get(strKey);
    if (m_Key != null)
    {
      IniItem m_Item = m_Key.GetItem(strItem);
      if (m_Item != null)
      {
        return(m_Item.GetValue());
      }
    }
    return("");
  }

  public int SetValue(String strKey, String strItem, String strValue)
  {
    return(SetValue(strKey, strItem, strValue, ""));
  }

  public int SetValue(String strKey, String strItem, String strValue, String strComment)
  {
	IniKey m_Key = (IniKey)m_Keys.get(strKey);
	if (m_Key != null)
	{
		IniItem m_Item = m_Key.GetItem(strItem);
		if (m_Item != null)
		{
			m_Item.SetValues(strItem, strValue,1, strComment);
		}
		else
		{ if(!("").equals(TrimString(strValue))){
			m_Key.AddItem(new IniItem(strItem, strValue, 1, strComment));						  
                                  }
		}
	}
	else
	{ if(!("").equals(TrimString(strValue))){
		m_Key = new IniKey(strKey);

		iMaxKeys++;
		m_Key.iSequence = iMaxKeys;
		m_Key.AddItem(new IniItem(strItem, strValue, 1, strComment));
	        m_Keys.put(strKey, m_Key);
                          }
	}
    return(0);
  }

  private boolean IsKey(String strLine)
  {
    int iPos1 = strLine.indexOf("[");
    int iPos2 = strLine.indexOf("]");

    //System.out.println("Pos1: " + iPos1 + ", Pos2: " + iPos2);

    if ((iPos1 == 0) && (iPos2 >= 0) && (iPos2 > iPos1))
    {
      return(true);
    }
    else
    {
      return(false);
    }

  }

  private String GetKeyFromLine(String strLine)
  {
    int iPos1 = strLine.indexOf("[");
    int iPos2 = strLine.indexOf("]");

    //System.out.println("Pos1: " + iPos1 + ", Pos2: " + iPos2);

    if ((iPos1 >= 0) && (iPos2 >= 0) && (iPos2 > iPos1))
    {
      return(strLine.substring(iPos1 + 1, iPos2));
    }
    else
    {
      return("");
    }
  }

  private String TrimString(String strSource)
  {
    while ((strSource.startsWith(" "))
      && (strSource.length() > 0))
      {
        strSource = strSource.substring(1, strSource.length());
      }

    while ((strSource.endsWith(" "))
        && (strSource.length() > 0))
      {
        strSource = strSource.substring(0, strSource.length() - 1);
      }

    return(strSource);
  }

  public void DumpParams()
  {
    for (Enumeration enumKeys = m_Keys.elements(); enumKeys.hasMoreElements();)
		{
      IniKey m_Key = (IniKey)enumKeys.nextElement();
      System.out.println("KEY: " + m_Key.GetKey());
      m_Key.DumpItems();
    }
  }
  
  private IniKey GetNextKey(int iSequ)
  {
	for (Enumeration enumKey = m_Keys.elements(); enumKey.hasMoreElements();)
	{
		IniKey m_Key = (IniKey)enumKey.nextElement();
		if (m_Key != null)
		{
			if (m_Key.iSequence == iSequ)
			{
				return(m_Key);
			}
		}
		else
		{
			return(null);
		}
	}
	return(null);
  }
}