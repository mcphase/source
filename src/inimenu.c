// this file contains rountines needed for the 
// main menu

#include "inipar.hpp"
#include "EZ.h"
#include "vector.h"

void checkquit();
inipar ini("mcphas.ini");



static void ExitCallBack(EZ_Widget *button, void *data);
static void DisplayCallBack(EZ_Widget *button, void *data);
static void floatCallback(EZ_Widget *entry, void *data);
static void doubleCallback(EZ_Widget *entry, void *data);
static void intCallBack(EZ_Widget *entry, void *data);


/*******************************************************************/
// necessary for EZ graphics server - program control
int EZ_initialized=0;
/**********************************************************************/

#define LINEWIDTH 25
#define TAB       55

// this sub checks if the user has pressed a button or changed any parameters
void checkquit()
{ char args[15]="spincf.cpp";  char *argv=&args[0];  
  static  EZ_Widget *button,*displayallbutton,*pausebutton,
                    *logqbutton,*frame,*xlabel,*ylabel,*xl1,*yl1,
                    *xT,*xHa,*xHb,*xHc,*xmin,*xmax,*xstep,
		    *yT,*yHa,*yHb,*yHc,*ymin,*ymax,*ystep,
		    *nofrndtries,*mclabel,
		    *maxnofmfloops,*maxnofmfloopslabel,
		    *maxstamf,*maxstamflabel,
		    *bigstep,*bigsteplabel,
		    *maxspinchange,*maxspinchangelabel;
		    
  static int washere=0;
  char xTstr[20],xHastr[20],xHbstr[20],xHcstr[20],xminstr[20],xmaxstr[20],xstepstr[20];
  char yTstr[20],yHastr[20],yHbstr[20],yHcstr[20],yminstr[20],ymaxstr[20],ystepstr[20];
  char nofrndtriesstr[20],maxnofmfloopsstr[20],maxstamfstr[20],bigstepstr[20],maxspinchangestr[20];
  
  snprintf(xTstr,sizeof(xTstr),"%g",ini.xv(0));
  snprintf(xHastr,sizeof(xHastr),"%g",ini.xv(1));
  snprintf(xHbstr,sizeof(xHbstr),"%g",ini.xv(2));
  snprintf(xHcstr,sizeof(xHcstr),"%g",ini.xv(3));
  snprintf(xminstr,sizeof(xminstr),"%g",ini.xmin);
  snprintf(xmaxstr,sizeof(xmaxstr),"%g",ini.xmax);
  snprintf(xstepstr,sizeof(xstepstr),"%g",ini.xstep);

  snprintf(yTstr,sizeof(yTstr),"%g",ini.yv(0));
  snprintf(yHastr,sizeof(yHastr),"%g",ini.yv(1));
  snprintf(yHbstr,sizeof(yHbstr),"%g",ini.yv(2));
  snprintf(yHcstr,sizeof(yHcstr),"%g",ini.yv(3));
  snprintf(yminstr,sizeof(yminstr),"%g",ini.ymin);
  snprintf(ymaxstr,sizeof(ymaxstr),"%g",ini.ymax);
  snprintf(ystepstr,sizeof(ystepstr),"%g",ini.ystep);
  
  snprintf(nofrndtriesstr,sizeof(nofrndtriesstr),"%i",ini.nofrndtries);
  snprintf(maxnofmfloopsstr,sizeof(maxnofmfloopsstr),"%i",ini.maxnofmfloops);
  snprintf(maxstamfstr,sizeof(maxstamfstr),"%g",ini.maxstamf);
  snprintf(bigstepstr,sizeof(bigstepstr),"%g",ini.bigstep);
  snprintf(maxspinchangestr,sizeof(maxspinchangestr),"%g",ini.maxspinchange);

  if (EZ_initialized==0){EZ_Initialize(1,&argv,0);EZ_initialized=1;}
  if (washere==0)
     {washere=1;// at first call of this routing - setup all buttons and fields
      frame = EZ_CreateWidget(EZ_WIDGET_FRAME,NULL,
                          EZ_WIDTH, 12*TAB, EZ_HEIGHT,8*LINEWIDTH,EZ_LABEL_JUSTIFICATION,EZ_LEFT,
			  EZ_TEXT_LINE_LENGTH,TAB*12, 
                          0);
      button = EZ_CreateWidget(EZ_WIDGET_CHECK_BUTTON, frame,
                           EZ_X,0,EZ_Y,0, 
                           EZ_LABEL_STRING,      "EXIT",
			   EZ_CHECK_BUTTON_ON_VALUE, 1,
			   EZ_CHECK_BUTTON_OFF_VALUE, 0,
                           EZ_CALLBACK,ExitCallBack, &ini.exit_mcphas,EZ_CHECK_BUTTON_ON_OFF,EZ_OFF,
                           0);
      EZ_DisplayWidget(button);
      if (ini.pause_mcphas)
      {pausebutton = EZ_CreateWidget(EZ_WIDGET_CHECK_BUTTON, frame,
                           EZ_X,2*TAB,EZ_Y,0, 
                           EZ_LABEL_STRING,      "Pause",
			   EZ_CHECK_BUTTON_ON_VALUE, 1,
			   EZ_CHECK_BUTTON_OFF_VALUE, 0,
                           EZ_CALLBACK,ExitCallBack, &ini.pause_mcphas,EZ_CHECK_BUTTON_ON_OFF,EZ_ON,
                           0);
      }else{
      pausebutton = EZ_CreateWidget(EZ_WIDGET_CHECK_BUTTON, frame,
                           EZ_X,2*TAB,EZ_Y,0, 
                           EZ_LABEL_STRING,      "Pause",
			   EZ_CHECK_BUTTON_ON_VALUE, 1,
			   EZ_CHECK_BUTTON_OFF_VALUE, 0,
                           EZ_CALLBACK,ExitCallBack, &ini.pause_mcphas,EZ_CHECK_BUTTON_ON_OFF,EZ_OFF,
                           0);
      }
      EZ_DisplayWidget(pausebutton);
      if (ini.displayall)
      {displayallbutton = EZ_CreateWidget(EZ_WIDGET_CHECK_BUTTON,frame,
                           EZ_X,4*TAB,EZ_Y,0,
                           EZ_LABEL_STRING,      "displayall",
			   EZ_CHECK_BUTTON_ON_VALUE, 1,
			   EZ_CHECK_BUTTON_OFF_VALUE, 0,
                           EZ_CALLBACK,DisplayCallBack, &ini.displayall,EZ_CHECK_BUTTON_ON_OFF,EZ_ON,
                           0);
      }else
      {displayallbutton = EZ_CreateWidget(EZ_WIDGET_CHECK_BUTTON,frame,
                           EZ_X,4*TAB,EZ_Y,0,
                           EZ_LABEL_STRING,      "displayall",
			   EZ_CHECK_BUTTON_ON_VALUE, 1,
			   EZ_CHECK_BUTTON_OFF_VALUE, 0,
                           EZ_CALLBACK,DisplayCallBack, &ini.displayall,EZ_CHECK_BUTTON_ON_OFF,EZ_OFF,
                           0);
      }
      EZ_DisplayWidget(displayallbutton);
      if(ini.logfevsQ)
      {logqbutton = EZ_CreateWidget(EZ_WIDGET_CHECK_BUTTON,frame,
                           EZ_X,6*TAB,EZ_Y,0,
                           EZ_LABEL_STRING,      "Log fe vs Q [mcphas.log]",
			   EZ_CHECK_BUTTON_ON_VALUE, 1,
			   EZ_CHECK_BUTTON_OFF_VALUE, 0,
                           EZ_CALLBACK,DisplayCallBack, &ini.logfevsQ,EZ_CHECK_BUTTON_ON_OFF,EZ_ON,
                           0);
      }else
      {logqbutton = EZ_CreateWidget(EZ_WIDGET_CHECK_BUTTON,frame,
                           EZ_X,6*TAB,EZ_Y,0,
                           EZ_LABEL_STRING,      "Log fe vs Q [mcphas.log]",
			   EZ_CHECK_BUTTON_ON_VALUE, 1,
			   EZ_CHECK_BUTTON_OFF_VALUE, 0,
                           EZ_CALLBACK,DisplayCallBack, &ini.logfevsQ,EZ_CHECK_BUTTON_ON_OFF,EZ_OFF,
                           0);
      }
      EZ_DisplayWidget(logqbutton);
     if (Norm(ini.xv)>0){
      xlabel = EZ_CreateWidget (EZ_WIDGET_LABEL, frame, 
                          EZ_X,0,EZ_Y,LINEWIDTH,
                          EZ_LABEL_STRING, "x:(T Ha Hb Hc) min max step (",0);
      EZ_DisplayWidget(xlabel);
      ylabel = EZ_CreateWidget (EZ_WIDGET_LABEL, frame, 
                          EZ_X,0,EZ_Y,2*LINEWIDTH,
                          EZ_LABEL_STRING, "y:(T Ha Hb Hc) min max step (",0);
      EZ_DisplayWidget(ylabel);
     }
     else
     {
      xlabel = EZ_CreateWidget (EZ_WIDGET_LABEL, frame, 
                          EZ_X,0,EZ_Y,LINEWIDTH,
                          EZ_LABEL_STRING, "x: min max  ",0);
      EZ_DisplayWidget(xlabel);
      ylabel = EZ_CreateWidget (EZ_WIDGET_LABEL, frame, 
                          EZ_X,0,EZ_Y,2*LINEWIDTH,
                          EZ_LABEL_STRING, "y: min max  ",0);
      EZ_DisplayWidget(ylabel);
     }
     if (Norm(ini.xv)>0){
      xT = EZ_CreateWidget (EZ_WIDGET_ENTRY, frame,
		           EZ_X,5*TAB,EZ_Y,LINEWIDTH,EZ_WIDTH,TAB-10,
			   EZ_ENTRY_STRING, xTstr,
			   EZ_CALLBACK,doubleCallback,&ini.xv[0],
			   0);
      EZ_DisplayWidget(xT);  

      xHa = EZ_CreateWidget (EZ_WIDGET_ENTRY, frame,
		           EZ_X,6*TAB,EZ_Y,LINEWIDTH,EZ_WIDTH,TAB-10,
			   EZ_ENTRY_STRING, xHastr,
			   EZ_CALLBACK,doubleCallback,&ini.xv[1],
			   0);
      EZ_DisplayWidget(xHa);  

      xHb = EZ_CreateWidget (EZ_WIDGET_ENTRY, frame,
		           EZ_X,7*TAB,EZ_Y,LINEWIDTH,EZ_WIDTH,TAB-10,
			   EZ_ENTRY_STRING, xHbstr,
			   EZ_CALLBACK,doubleCallback,&ini.xv[2],
			   0);
      EZ_DisplayWidget(xHb);  
        
      xHc = EZ_CreateWidget (EZ_WIDGET_ENTRY, frame,
		           EZ_X,8*TAB,EZ_Y,LINEWIDTH,EZ_WIDTH,TAB-10,
			   EZ_ENTRY_STRING, xHcstr,
			   EZ_CALLBACK,doubleCallback,&ini.xv[3],
			   0);
      EZ_DisplayWidget(xHc);  
      
      xl1 = EZ_CreateWidget (EZ_WIDGET_LABEL, frame, 
                          EZ_X,9*TAB-20,EZ_Y,LINEWIDTH,
                          EZ_LABEL_STRING, ")",0);
      EZ_DisplayWidget(xl1);
     }
      xmin = EZ_CreateWidget (EZ_WIDGET_ENTRY, frame,
		           EZ_X,9*TAB,EZ_Y,LINEWIDTH,EZ_WIDTH,TAB-10,
			   EZ_ENTRY_STRING, xminstr,
			   EZ_CALLBACK,floatCallback,&ini.xmin,
			   0);
      EZ_DisplayWidget(xmin);  
      xmax = EZ_CreateWidget (EZ_WIDGET_ENTRY, frame,
		           EZ_X,10*TAB,EZ_Y,LINEWIDTH,EZ_WIDTH,TAB-10,
			   EZ_ENTRY_STRING, xmaxstr,
			   EZ_CALLBACK,floatCallback,&ini.xmax,
			   0);
      EZ_DisplayWidget(xmax);  
     if (Norm(ini.xv)>0){
      xstep = EZ_CreateWidget (EZ_WIDGET_ENTRY, frame,
		           EZ_X,11*TAB,EZ_Y,LINEWIDTH,EZ_WIDTH,TAB-10,
			   EZ_ENTRY_STRING, xstepstr,
			   EZ_CALLBACK,floatCallback,&ini.xstep,
			   0);
      EZ_DisplayWidget(xstep);      
      yT = EZ_CreateWidget (EZ_WIDGET_ENTRY, frame,
		           EZ_X,5*TAB,EZ_Y,2*LINEWIDTH,EZ_WIDTH,TAB-10,
			   EZ_ENTRY_STRING, yTstr,
			   EZ_CALLBACK,doubleCallback,&ini.yv[0],
			   0);
      EZ_DisplayWidget(yT);  

      yHa = EZ_CreateWidget (EZ_WIDGET_ENTRY, frame,
		           EZ_X,6*TAB,EZ_Y,2*LINEWIDTH,EZ_WIDTH,TAB-10,
			   EZ_ENTRY_STRING, yHastr,
			   EZ_CALLBACK,doubleCallback,&ini.yv[1],
			   0);
      EZ_DisplayWidget(yHa);  

      yHb = EZ_CreateWidget (EZ_WIDGET_ENTRY, frame,
		           EZ_X,7*TAB,EZ_Y,2*LINEWIDTH,EZ_WIDTH,TAB-10,
			   EZ_ENTRY_STRING, yHbstr,
			   EZ_CALLBACK,doubleCallback,&ini.yv[2],
			   0);
      EZ_DisplayWidget(yHb);  
        
      yHc = EZ_CreateWidget (EZ_WIDGET_ENTRY, frame,
		           EZ_X,8*TAB,EZ_Y,2*LINEWIDTH,EZ_WIDTH,TAB-10,
			   EZ_ENTRY_STRING, yHcstr,
			   EZ_CALLBACK,doubleCallback,&ini.yv[3],
			   0);
      EZ_DisplayWidget(yHc);
        
      yl1 = EZ_CreateWidget (EZ_WIDGET_LABEL, frame, 
                          EZ_X,9*TAB-20,EZ_Y,2*LINEWIDTH,
                          EZ_LABEL_STRING, ")",0);
      EZ_DisplayWidget(yl1);
     }
      ymin = EZ_CreateWidget (EZ_WIDGET_ENTRY, frame,
		           EZ_X,9*TAB,EZ_Y,2*LINEWIDTH,EZ_WIDTH,TAB-10,
			   EZ_ENTRY_STRING, yminstr,
			   EZ_CALLBACK,floatCallback,&ini.ymin,
			   0);
      EZ_DisplayWidget(ymin);  
      ymax = EZ_CreateWidget (EZ_WIDGET_ENTRY, frame,
		           EZ_X,10*TAB,EZ_Y,2*LINEWIDTH,EZ_WIDTH,TAB-10,
			   EZ_ENTRY_STRING, ymaxstr,
			   EZ_CALLBACK,floatCallback,&ini.ymax,
			   0);
      EZ_DisplayWidget(ymax);  
     if (Norm(ini.xv)>0){
      ystep = EZ_CreateWidget (EZ_WIDGET_ENTRY, frame,
		           EZ_X,11*TAB,EZ_Y,2*LINEWIDTH,EZ_WIDTH,TAB-10,
			   EZ_ENTRY_STRING, ystepstr,
			   EZ_CALLBACK,floatCallback,&ini.ystep,
			   0);
      EZ_DisplayWidget(ystep);      
     }
      mclabel = EZ_CreateWidget (EZ_WIDGET_LABEL, frame, 
                          EZ_X,0,EZ_Y,3*LINEWIDTH,
			  EZ_TEXT_LINE_LENGTH,TAB*12, 
                          EZ_LABEL_STRING, "No. of Monte-Carlo tries for each conf.",0);
      EZ_DisplayWidget(mclabel);
      nofrndtries = EZ_CreateWidget (EZ_WIDGET_ENTRY, frame,
		           EZ_X,7*TAB,EZ_Y,3*LINEWIDTH,EZ_WIDTH,TAB,
			   EZ_ENTRY_STRING, nofrndtriesstr,
			   EZ_CALLBACK,intCallBack,&ini.nofrndtries,
			   0);
      EZ_DisplayWidget(nofrndtries);      

      maxnofmfloopslabel = EZ_CreateWidget (EZ_WIDGET_LABEL, frame, 
                          EZ_X,0,EZ_Y,4*LINEWIDTH,
			  EZ_TEXT_LINE_LENGTH,TAB*12, 
                          EZ_LABEL_STRING, "Max.No. of MF loops  to end MF process",0);
      EZ_DisplayWidget(maxnofmfloopslabel);
      maxnofmfloops = EZ_CreateWidget (EZ_WIDGET_ENTRY, frame,
		           EZ_X,7*TAB,EZ_Y,4*LINEWIDTH,EZ_WIDTH,TAB,
			   EZ_ENTRY_STRING, maxnofmfloopsstr,
			   EZ_CALLBACK,intCallBack,&ini.maxnofmfloops,
			   0);
      EZ_DisplayWidget(maxnofmfloops);      

      maxstamflabel = EZ_CreateWidget (EZ_WIDGET_LABEL, frame, 
                          EZ_X,0,EZ_Y,5*LINEWIDTH,
			  EZ_TEXT_LINE_LENGTH,TAB*12, 
                          EZ_LABEL_STRING, "Max. stand. dev. of spinchange to end MF selfconsistency process",0);
      EZ_DisplayWidget(maxstamflabel);
      maxstamf = EZ_CreateWidget (EZ_WIDGET_ENTRY, frame,
		           EZ_X,10*TAB+10,EZ_Y,5*LINEWIDTH,EZ_WIDTH,TAB,
			   EZ_ENTRY_STRING, maxstamfstr,
			   EZ_CALLBACK,floatCallback,&ini.maxstamf,
			   0);
      EZ_DisplayWidget(maxstamf);      

      bigsteplabel = EZ_CreateWidget (EZ_WIDGET_LABEL, frame, 
                          EZ_X,0,EZ_Y,6*LINEWIDTH,
			  EZ_TEXT_LINE_LENGTH,TAB*12, 
                          EZ_LABEL_STRING, "a step [0,1] (=actual/calc.MF step)",0);
      EZ_DisplayWidget(bigsteplabel);
      bigstep = EZ_CreateWidget (EZ_WIDGET_ENTRY, frame,
		           EZ_X,6*TAB,EZ_Y,6*LINEWIDTH,EZ_WIDTH,TAB,
			   EZ_ENTRY_STRING, bigstepstr,
			   EZ_CALLBACK,floatCallback,&ini.bigstep,
			   0);
      EZ_DisplayWidget(bigstep);      

      maxspinchangelabel = EZ_CreateWidget (EZ_WIDGET_LABEL, frame, 
                          EZ_X,0,EZ_Y,7*LINEWIDTH,
			  EZ_TEXT_LINE_LENGTH,TAB*12, 
                          EZ_LABEL_STRING, "Max. total spinchange to end MF process",0);
      EZ_DisplayWidget(maxspinchangelabel);
      maxspinchange = EZ_CreateWidget (EZ_WIDGET_ENTRY, frame,
		           EZ_X,7*TAB,EZ_Y,7*LINEWIDTH,EZ_WIDTH,TAB,
			   EZ_ENTRY_STRING, maxspinchangestr,
			   EZ_CALLBACK,floatCallback,&ini.maxspinchange,
			   0);
      EZ_DisplayWidget(maxspinchange);      

    }
    else
    {EZ_ServiceEvents();
     if (ini.exit_mcphas){exit(0);} // print parameters to file and exit
    }
}

// here follow the different routines which are called when
// the user presses a button

static void ExitCallBack(EZ_Widget *button, void *data)
{ int *state = (int *)data;
   EZ_DisplayWidget(button);
 *state=  EZ_GetCheckButtonState(button,state);
  ini.print();
}


static void floatCallback(EZ_Widget *entry, void * data)
{ char *state;
  float * number= (float *)data;
  state =  EZ_GetEntryString(entry);
 (*number)=(float)strtod (state,NULL);
 ini.print();
}

static void doubleCallback(EZ_Widget *entry, void * data)
{ char *state;
  double * number= (double *)data;
 state =  EZ_GetEntryString(entry);
 (*number)=(double)strtod (state,NULL);
 ini.print();
}

static void intCallBack(EZ_Widget *entry, void *data)
{ char *state;
  int * number= (int *)data;
 state =  EZ_GetEntryString(entry);
 strtod (state,NULL);
 *number=(int)strtod (state,NULL);
 ini.print();
}

static void DisplayCallBack(EZ_Widget *button, void *data)
{ int *state = (int *)data;
   EZ_DisplayWidget(button);
 *state=  EZ_GetCheckButtonState(button,state);
 ini.print();
}

int main()
{ini.exit_mcphas=0;ini.print(); // if exit was 1 - save parameters and set exit=0
 while(1!=0)
 {sleep(1);
 checkquit();
 }
}
