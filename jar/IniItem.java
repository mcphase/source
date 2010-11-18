
//Title:        McPhas INI-Configurator
//Version:      
//Copyright:    Copyright (c) 1998
//Author:       Stefan Rotter
//Company:      Stefan Rotter
//Description:  Configurator for mcphas INI-File


//package mcphas;

public class IniItem {

public final int TYPE_IS_STRING = 1;
public final int TYPE_IS_NUMBER = 2;

public int iSequence;

private String m_strItem;
private String m_strValue;
private int m_iType;
private String m_strComment;

  public IniItem() {
    m_strItem = "";
    m_strValue = "";
	m_strComment = "";
    m_iType = TYPE_IS_STRING;
  }

  public IniItem(String strItem, String strValue, int iType) {
    m_strItem = strItem;
    m_strValue = strValue;
	m_strComment = "";
    m_iType = iType;
  }

  public IniItem(String strItem, String strValue, int iType, String strComment) {
    m_strItem = strItem;
    m_strValue = strValue;
	if (strComment.length() > 0)
	{
		m_strComment = "\n# " + strComment + "\n";
	}
	else
	{
		m_strComment = "";
	}
    m_iType = iType;
  }

  public void SetValues(String strItem, String strValue, int iType) {
    m_strItem = strItem;
    m_strValue = strValue;
    m_iType = iType;
	m_strComment = "";
  }

  public void SetValues(String strItem, String strValue, int iType, String strComment) {
    m_strItem = strItem;
    m_strValue = strValue;
    m_iType = iType;
	if (strComment.length() > 0)
	{
		m_strComment = "\n# " + strComment + "\n";
	}
	else
	{
		m_strComment = "";
	}
  }

  public String GetItem() {
    return(m_strItem);
  }

  public String GetValue() {
    return(m_strValue);
  }

  public int GetType() {
    return(m_iType);
  }

  public String GetComment() {
    return(m_strComment);
  }
}