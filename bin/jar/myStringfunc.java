class myStringfunc {

public myStringfunc(){}

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