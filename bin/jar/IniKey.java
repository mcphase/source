
//Title:        McPhas INI-Configurator
//Version:
//Copyright:    Copyright (c) 1998
//Author:       Stefan Rotter
//Company:      Stefan Rotter
//Description:  Configurator for mcphas INI-File



//package mcphas;

import java.util.Hashtable;
import java.util.Enumeration;

public class IniKey {

public int iSequence;

public int iMaxItems;

private Hashtable m_Items;
private String m_Key;

  public IniKey() {
    m_Items = new Hashtable();
    m_Key = "";
	iMaxItems = 0;
  }

  public IniKey(String strKey) {
    m_Items = new Hashtable();
    m_Key = strKey;
	iMaxItems = 0;
  }

  public int AddItem(IniItem Item)
  { Item.iSequence=iMaxItems+1; // introduced by Martin
    m_Items.put(Item.GetItem(), Item);

	iMaxItems++;
    return(0);
  }

  public int RemoveItem(String strKey)
  {
    m_Items.remove(strKey);
	iMaxItems--;
    return(0);
  }

  public int RefreshItem(IniItem Item)
  {

    return(0);
  }

  public IniItem GetItem(String strKey)
  {
    return((IniItem)m_Items.get(strKey));
  }

  public String GetKey()
  {
    return(m_Key);
  }

  public void DumpItems()
  {
    for (Enumeration enumItems = m_Items.elements(); enumItems.hasMoreElements();)
		{
      IniItem m_Item = (IniItem)enumItems.nextElement();
      System.out.println("ITEM: " + m_Item.GetItem() + ": " + m_Item.GetValue());

    }

  }

  public String GetItemsForWrite()
  {
	String strReturn = "";
	
	for (int iSequItem = 1; iSequItem <= iMaxItems;++iSequItem)
	{
		IniItem m_Item = GetNextItem(iSequItem);
		if (m_Item != null)
		{
			strReturn = strReturn + m_Item.GetComment() + m_Item.GetItem() + "=" + m_Item.GetValue() + "\n";
		}
	}
	return(strReturn);
  }

  private IniItem GetNextItem(int iSequ)
  {
	for (Enumeration enumItems = m_Items.elements(); enumItems.hasMoreElements();)
	{
		IniItem m_Item = (IniItem)enumItems.nextElement();
		if (m_Item != null)
		{
			if (m_Item.iSequence == iSequ)
			{
				return(m_Item);
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