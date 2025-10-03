class myStringfunc {

public myStringfunc(){}

public int NofCols (String strLine)
{String sx;
 // select N
      sx=TrimString(strLine);
      int cx =0;
      if(sx.length() > 0){++cx;}   
      int iPos = sx.indexOf(" ");
      while (iPos >= 0)
      {++cx;
       sx=sx.substring(iPos);
       sx=TrimString(sx);iPos = sx.indexOf(" ");}
       return (cx);
}   

public String NthWord (String strLine,int N)
 {String sx;
 // select N
      sx=TrimString(strLine);
      int cx =N-1;
      int colxnotfound=0;
      while (cx>0)
      {--cx;
       int iPos = sx.indexOf(" ");
       if (iPos < 0)
       {
         colxnotfound=1;
        continue;
       }
       sx=sx.substring(iPos);
       sx=TrimString(sx);
      }
       cx=sx.indexOf(" ");
       if (cx>0) {sx=sx.substring(0,cx);}
  if (colxnotfound==1){sx=sx.substring(0,0);}
  return (sx);
 }

public String FirstWord(String strSource)
 { String fw;
  fw=TrimString(strSource);
       int iPos = fw.indexOf(" ");
       if (iPos >= 0)
       {
       fw=strSource.substring(0,iPos);
       fw=TrimString(fw);
       }
 return(fw);
 }

public String DataCol(String strSource)
 { String fw;String f;
  fw=TrimString(strSource);f=fw;
       int iPos = fw.indexOf("e");int tipos=iPos;
       // try to exclude "exp" function from being detected es error separator
       while(iPos>=0&&fw.substring(iPos+1).startsWith("x")){fw=fw.substring(iPos+1); iPos = fw.indexOf("e");tipos+=iPos+1;}
     
       if (iPos >= 0)
       {
       fw=f.substring(0,tipos);
       fw=TrimString(fw);
       }
       else
       {iPos = f.indexOf("b");
      // try to exclude "abs" function from being detected es error separator
       while(iPos>=0&&fw.substring(iPos+1).startsWith("s")){fw=fw.substring(iPos+1); iPos = fw.indexOf("b");tipos+=iPos+1;}
     
       if (iPos >= 0)
       {
       fw=f.substring(0,iPos);
       fw=TrimString(fw);
       } else
       {fw=f;}
       }
 return(fw);
 }

public String ErrorCol(String strSource)
 { String fw;
  fw=TrimString(strSource);
       int iPos = fw.indexOf("e"); // try to exclude "exp" function from being detected es error separator
       while(iPos>=0&&fw.substring(iPos+1).startsWith("x")){fw=fw.substring(iPos+1); iPos = fw.indexOf("e");}
       if (iPos >= 0)
       {
       fw=fw.substring(iPos+1);
       fw=TrimString(fw);
       }
       else
       {fw="0";}
 return(fw);
 }

public String BubbleCol(String strSource)
 { String fw;
  fw=TrimString(strSource);
       int iPos = fw.indexOf("b");// try to exclude "abs" function from being detected es bubble separator
       while(iPos>=0&&fw.substring(iPos+1).startsWith("s")){fw=fw.substring(iPos+1); iPos = fw.indexOf("b");}
       if (iPos >= 0)
       {
       fw="b"+fw.substring(iPos+1);
       fw=TrimString(fw);
       }
       else
       {fw="0";}
 return(fw);
 }


public String DropWord(String strSource)
 {String fw;
  fw=TrimString(strSource);
       int iPos = fw.indexOf(" ");
       if (iPos >= 0)
       {
       fw=strSource.substring(iPos);
       fw=TrimString(fw);
       }
       else
       {fw="";}
 return(fw);
 }

public String TrimString(String strSource)
 {while ((strSource.startsWith(" ")) || (strSource.startsWith("\""))
      && (strSource.length() > 0))
      {
        strSource = strSource.substring(1, strSource.length());
      }
    while ((strSource.endsWith(" "))|| (strSource.startsWith("\""))
        && (strSource.length() > 0))
      {
        strSource = strSource.substring(0, strSource.length() - 1);
      }
    return(strSource);
 }// TrimString

} // myStringfunc