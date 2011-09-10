/*
    Mac OS Portions of the chart window class

    Sergei L. Kosakovsky Pond, Spring 2000 - December 2002.
*/

#include "HYChartWindow.h"
#include "HYCanvas.h"
#include "HYUtils.h"
#include "HYPulldown.h"

#include "math.h"
#include "scrap.h"

extern   _Parameter     pi_const;


//__________________________________________________________________

void _HYChartWindow::_SetMenuBar(void)
{
    _HYWindow::_SetMenuBar();
    MenuHandle  t = GetMenuHandle (130);
    EnableMenuItem (t,4);
    EnableMenuItem (t,5);
    EnableMenuItem (t,8);
    MenuHandle chartMenu = GetMenuHandle (HY_CHART_WINDOW_MENU_ID);
    if (!chartMenu) {
        chartMenu = NewMenu(HY_CHART_WINDOW_MENU_ID,"\pChart");
        if (!chartMenu) {
            warnError (-108);
        }

        MenuHandle saveMenu  =  NewMenu(HY_CHART_WINDOW_HMENU_ID,"\pSave"),
                   printMenu =  NewMenu(HY_CHART_WINDOW_HMENU_ID+1,"\pPrint"),
                   fontMenu  =  NewMenu(HY_CHART_WINDOW_HMENU_ID+2,"\pFonts"),
                   procMenu  =  NewMenu(HY_CHART_WINDOW_HMENU_ID+3,"\pData Processing");

        InsertMenuItem  (saveMenu, "\pSave Chart.../S",10000);
        InsertMenuItem  (saveMenu, "\pSave Graphic...",10000);
        InsertMenuItem  (saveMenu, "\pSave Table...",10000);
        InsertMenuItem  (printMenu,"\pPrint Chart.../P",10000);
        InsertMenuItem  (printMenu,"\pPrint Table...",10000);
        InsertMenuItem  (fontMenu, "\pTickmark  Font...",10000);
        InsertMenuItem  (fontMenu, "\pLegend Font...",10000);
        InsertMenuItem  (fontMenu, "\pAxis Label  Font...",10000);

        InsertMenuItem  (chartMenu,"\pChart Name...",10000);
        InsertMenuItem  (chartMenu,"\pChart Options...",10000);
        InsertMenuItem  (chartMenu, "\pFonts",10000);
        InsertMenuItem  (chartMenu, "\p(-",10000);
        InsertMenuItem  (chartMenu, "\pData Processing",10000);
        InsertMenu      (saveMenu,hierMenu);
        InsertMenu      (printMenu,hierMenu);
        InsertMenu      (fontMenu,hierMenu);
        InsertMenu      (procMenu,hierMenu);
        InsertMenu      (chartMenu,132);

        SetItemCmd (chartMenu,3,hMenuCmd);
        SetItemMark(chartMenu,3,HY_CHART_WINDOW_HMENU_ID+2);
        SetItemCmd (chartMenu,5,hMenuCmd);
        SetItemMark(chartMenu,5,HY_CHART_WINDOW_HMENU_ID+3);
        if (chartProcessors.lLength == 0) {
            DisableMenuItem (chartMenu,5);
        } else {
            Str255    buffer;
            for (long k=0; k<chartProcessors.lLength; k++) {
                _String *thisItem = (_String*)chartProcessors (k),
                         chopped = thisItem->Cut (thisItem->FindBackwards (':',0,-1)+1,-1);
                StringToStr255  (chopped,buffer);
                InsertMenuItem  (procMenu, buffer,10000);
            }
        }

        t = GetMenuHandle (129);
        SetItemCmd (t,4,hMenuCmd);
        SetItemMark(t,4,HY_CHART_WINDOW_HMENU_ID);
        SetItemCmd (t,8,hMenuCmd);
        SetItemMark(t,8,HY_CHART_WINDOW_HMENU_ID+1);
    }
    InvalMenuBar();
}


//__________________________________________________________________

void _HYChartWindow::_UnsetMenuBar(void)
{
    MenuHandle chartMenu    = GetMenuHandle (HY_CHART_WINDOW_MENU_ID),
               saveMenu   = GetMenuHandle (HY_CHART_WINDOW_HMENU_ID),
               printMenu    = GetMenuHandle (HY_CHART_WINDOW_HMENU_ID+1),
               fontMenu        = GetMenuHandle (HY_CHART_WINDOW_HMENU_ID+2),
               procMenu      = GetMenuHandle (HY_CHART_WINDOW_HMENU_ID+3),
               fMenu       = GetMenuHandle (129);

    DeleteMenu (HY_CHART_WINDOW_MENU_ID);
    DeleteMenu (HY_CHART_WINDOW_HMENU_ID);
    DeleteMenu (HY_CHART_WINDOW_HMENU_ID+1);
    DeleteMenu (HY_CHART_WINDOW_HMENU_ID+2);
    DeleteMenu (HY_CHART_WINDOW_HMENU_ID+3);
    DisposeMenu (chartMenu);
    DisposeMenu (printMenu);
    DisposeMenu (saveMenu);
    DisposeMenu (fontMenu);
    DisposeMenu (procMenu);
    SetItemCmd (fMenu,4,'S');
    SetItemCmd (fMenu,8,'P');
    SetItemMark(fMenu,4,noMark);
    SetItemMark(fMenu,8,noMark);
    _HYWindow::_UnsetMenuBar();
}


//__________________________________________________________________

void        _HYChartWindow::_PrintChart(void)
{
    GrafPtr     savePort;
#ifdef      TARGET_API_MAC_CARBON
    PMRect prRect;
#else
    TPrStatus   prStatus;
    TPPrPort    printPort;
    OSErr       err;
#endif

#ifdef TARGET_API_MAC_CARBON
    OSStatus theStatus;
    Boolean isAccepted;

    PMPrintSession hyPC;
    theStatus = PMCreateSession(&hyPC);
    if (theStatus != noErr) {
        return;
    }
#endif

    if (!InitPrint(hyPC)) {
        _String errMsg ("Could not initialize printing variables.");
        WarnError (errMsg);
        terminateExecution = false;
        return;
    }

    GetPort(&savePort);

#ifdef TARGET_API_MAC_CARBON
    if (gPrintSettings != kPMNoPrintSettings) {
        theStatus = PMSessionValidatePrintSettings(hyPC,gPrintSettings, kPMDontWantBoolean);
    } else {
        theStatus = PMCreatePrintSettings(&gPrintSettings);

        if ((theStatus == noErr) && (gPrintSettings != kPMNoPrintSettings)) {
            theStatus = PMSessionDefaultPrintSettings(hyPC,gPrintSettings);
        }
    }

    if (theStatus == noErr) {
        theStatus = PMSessionPrintDialog(hyPC,gPrintSettings, gPageFormat, &isAccepted);

        if (isAccepted) {
            theStatus = PMGetAdjustedPageRect(gPageFormat, &prRect);
            if (theStatus != noErr) {
                return;
            }

            theStatus = PMSessionBeginDocument(hyPC,gPrintSettings, gPageFormat);
            if (theStatus != noErr) {
                return;
            }

            long     printW    = prRect.right-prRect.left-2,
                     printH    = prRect.bottom-prRect.top-2;

            UInt32   startPage,
                     endPage;

            PMGetFirstPage (gPrintSettings,&startPage);
            PMGetLastPage  (gPrintSettings,&endPage);
#else
    PrOpen();
    if (err=PrError()) {
        _String errMsg ("Could not print the chart. Error Code:");
        errMsg = errMsg & (long)err;
        WarnError (errMsg);
        terminateExecution = false;
        return;
    }

    if (PrJobDialog(prRecHdl)) {
        printPort = PrOpenDoc(prRecHdl, nil, nil);
        SetPort((GrafPtr)printPort);
        long     printW = (*prRecHdl)->prInfo.rPage.right-2,
                 printH = (*prRecHdl)->prInfo.rPage.bottom-2,
                 startPage = (*prRecHdl)->prJob.iFstPage,
                 endPage = (*prRecHdl)->prJob.iLstPage;
#endif
#ifdef TARGET_API_MAC_CARBON
            theStatus = PMSessionBeginPage(hyPC, gPageFormat, NULL);
            if (theStatus != noErr) {
                return;
            }
            GrafPtr ppPort;
            PMSessionGetGraphicsContext (hyPC, NULL, (void**)&ppPort);
            SetPort (ppPort);
#else
            PrOpenPage      (printPort, nil);
#endif

            _HYRect     viewRect  = ((_HYStretchCanvas*)GetObject (0))->GetCanvasSize();
            _Parameter  aspectRatio = viewRect.right/(_Parameter)viewRect.bottom;
            _HYRect     printRect = {0,0,printH, printH*aspectRatio,0};

            if (printRect.right > printW) {
                aspectRatio = printW/(_Parameter)printRect.right;
                printRect.right = printW;
                printRect.bottom *= aspectRatio;
            }

            DrawChart   (&printRect);

#ifdef TARGET_API_MAC_CARBON
            theStatus = PMSessionEndPage(hyPC);
            if (theStatus != noErr) {
                return;
            }
#else
            PrClosePage     (printPort);
#endif


#ifdef TARGET_API_MAC_CARBON
            theStatus = PMSessionEndDocument(hyPC);
            SetPort(savePort);
            if (theStatus == noErr) {
                if (gFlattenedFormat != NULL) {
                    DisposeHandle(gFlattenedFormat);
                    gFlattenedFormat = NULL;
                }

                theStatus = PMFlattenPageFormat(gPageFormat, &gFlattenedFormat);
            }

            if (theStatus == noErr) {
                if (gFlattenedSettings != NULL) {
                    DisposeHandle(gFlattenedSettings);
                    gFlattenedSettings = NULL;
                }

                theStatus = PMFlattenPrintSettings(gPrintSettings, &gFlattenedSettings);
            }

            if (gPageFormat != kPMNoPageFormat) {
                theStatus = PMRelease(gPageFormat);
                gPageFormat = kPMNoPageFormat;
            }

            if (gPrintSettings != kPMNoPrintSettings) {
                theStatus = PMRelease(gPrintSettings);
                gPrintSettings = kPMNoPrintSettings;
            }

            theStatus = PMRelease(hyPC);

#else
            PrCloseDoc(printPort);
            if (((*prRecHdl)->prJob.bJDocLoop = bSpoolLoop) && (!PrError() ) ) {
                PrPicFile(prRecHdl, nil, nil, nil, &prStatus);
            }
#endif
        }
#ifdef TARGET_API_MAC_CARBON
        else {
            theStatus = PMRelease(hyPC);
        }
#endif

#ifdef TARGET_API_MAC_CARBON
    }
#else
        PrClose();
        SetPort(savePort);
#endif


}

//__________________________________________________________________


bool        _HYChartWindow::_ProcessMenuSelection (long msel)
{
    long        menuChoice = msel&0x0000ffff;

    HiliteMenu(0);
    InvalMenuBar();

    switch (msel/0xffff) {
    case HY_CHART_WINDOW_MENU_ID: { // chart menu
        if (menuChoice==1) {
            RenameChartWindow ();
            return true;
        } else if (menuChoice==2) {
            HandleChartOptions ();
            return true;
        }
    }
    case HY_CHART_WINDOW_HMENU_ID: { // save menu
        DoSave (menuChoice-1);
        return true;
    }
    case HY_CHART_WINDOW_HMENU_ID+1: { // print menu
        DoPrint (menuChoice-1);
        return true;
    }
    case HY_CHART_WINDOW_HMENU_ID+2: { // font menu
        DoChangeFont (menuChoice-1);
        return true;
    }
    case HY_CHART_WINDOW_HMENU_ID+3: { // proc menu
        ExecuteProcessor (menuChoice-1);
        return true;
    }
    }

    return _HYTWindow::_ProcessMenuSelection(msel);
}

//__________________________________________________________________

bool _HYChartWindow::_ProcessOSEvent (Ptr vEvent)
{
    EventRecord*    theEvent = (EventRecord*)vEvent;

    if (!_HYTWindow::_ProcessOSEvent (vEvent)) {
        _HYPullDown *p1 = (_HYPullDown*)GetObject (4);

        if ((theEvent->what==mouseDown)&&(p1->GetSelection()>=8)&&(ySeries.lLength)) {
            Point localClick = theEvent->where;
            GrafPtr savedPort;
            GetPort(&savedPort);
#ifdef OPAQUE_TOOLBOX_STRUCTS
            SetPort(GetWindowPort(theWindow));
#else
            SetPort(theWindow);
#endif
            GlobalToLocal (&localClick);
            int   c = FindClickedCell(localClick.h,localClick.v),ch,cv;
            if (c<0) {
                return false;
            }
            if (c==0) { // the chart
                ch = localClick.h-componentL.lData[0];
                cv = localClick.v-componentT.lData[0];

                Point oldPt = {-1,-1};

                if (StillDown()) {
                    while (WaitMouseUp()) {
                        GetMouse(&localClick);
                        if (oldPt.h==-1) {
                            oldPt = localClick;
                        }
                        if (DeltaPoint (oldPt,localClick)) {
                            bool        redraw = false;

                            _Parameter  stepper = pi_const/180.;

                            if (abs(localClick.h-oldPt.h)>abs(localClick.v-oldPt.v)) {
                                stepper *= 1+log (fabs(localClick.h-oldPt.h))/log(2.0);
                                if (localClick.h-oldPt.h<0) {
                                    if (xyAngle>0.0) {
                                        xyAngle -= stepper;
                                        if (xyAngle<0) {
                                            xyAngle = 0;
                                        }
                                        redraw = true;
                                    }
                                } else if (xyAngle<pi_const/2) {
                                    xyAngle += stepper;
                                    if (xyAngle>pi_const/2) {
                                        xyAngle = pi_const/2;
                                    }
                                    redraw = true;
                                }
                            } else {
                                stepper *= 1+log (fabs(localClick.v-oldPt.v))/log(2.0);
                                if (localClick.v-oldPt.v>0) {
                                    if (zAngle<pi_const/2) {
                                        zAngle += stepper;
                                        if (zAngle>pi_const/2) {
                                            zAngle = pi_const/2;
                                        }
                                        redraw = true;
                                    }
                                } else if (zAngle>0.0) {
                                    zAngle -= stepper;
                                    if (zAngle<0) {
                                        zAngle = 0;
                                    }
                                    redraw = true;
                                }

                            }

                            if (redraw) {
                                ComputeProjectionSettings();
                                projectionMatrix = ComputeProjectionMatrix   ();
                                forceUpdateForScrolling = true;
                                DrawChart();
                                forceUpdateForScrolling = false;
                            }
                        }
                        oldPt = localClick;
                    }
                }
                return true;
            }
        }
        return false;
    }
    return true;
}

//__________________________________________________________________

void _HYChartWindow::_CopyChart (void)
{
#ifdef TARGET_API_MAC_CARBON
    ClearCurrentScrap();
#else
    ZeroScrap();
#endif

    _HYStretchCanvas    *sc = (_HYStretchCanvas*)GetObject (0);
    Rect                bRect;
    _HYRect             bounds;

    bounds.left         = bounds.top =
                              bRect.left            = bRect.top = 0;

    bounds.right        =   bRect.right         = sc->_HYComponent::GetMaxW();
    bounds.bottom       =   bRect.bottom        = sc->_HYComponent::GetMaxH();

    PicHandle    pic    = OpenPicture (&bRect);

    if (pic) {
        PixPatHandle         white = NewPixPat ();
        checkPointer         (white);
        RGBColor             wc = (RGBColor) {
            0xffff,0xffff,0xffff
        };
        MakeRGBPat           (white,&wc);
        FillCRect            (&bRect,white);
        DrawChart            (&bounds);
        DisposePixPat        (white);

        ClosePicture ();
        HLock   ((Handle)pic);
#ifdef TARGET_API_MAC_CARBON
        ClearCurrentScrap();
        ScrapRef         theScrapRef;
        GetCurrentScrap(&theScrapRef);
        PutScrapFlavor(theScrapRef, 'PICT', kScrapFlavorMaskNone,GetHandleSize((Handle)pic),*pic);
#else
        PutScrap (GetHandleSize((Handle)pic),'PICT',*pic);
#endif
        KillPicture (pic);
    } else {
        _String wMsg ("Failed to open a picture in _CopyChart - probably low on memory\n");
        WarnError (wMsg);
        terminateExecution = false;
    }
}

//__________________________________________________________________

void _HYDistributionChartWindow::_SetMenuBar(void)
{
    _HYChartWindow::_SetMenuBar();
    MenuHandle dMenu = GetMenuHandle (HY_CHART_WINDOW_MENU_ID+1);
    if (!dMenu) {
        dMenu = NewMenu(HY_CHART_WINDOW_MENU_ID+1,"\pCategories");
        if (!dMenu) {
            warnError (-108);
        }


        InsertMenuItem  (dMenu, "\pDefine New Variable",10000);
        InsertMenuItem  (dMenu, "\pDelete Variable",  10000);
        InsertMenuItem  (dMenu, "\pConditional Distribution",  10000);

        if (distribProcessors.lLength > 0) {
            InsertMenuItem  (dMenu, "\p(-",10000);
            Str255    buffer;
            for (long k=0; k<distribProcessors.lLength; k++) {
                _String *thisItem = (_String*)distribProcessors (k),
                         chopped = thisItem->Cut (thisItem->FindBackwards (':',0,-1)+1,-1);
                StringToStr255  (chopped,buffer);
                InsertMenuItem  (dMenu, buffer,10000);
            }
        }

        InsertMenu      (dMenu,132);
    }
    InvalMenuBar();
}

//__________________________________________________________________

void _HYDistributionChartWindow::_UnsetMenuBar(void)
{
    MenuHandle dMenu = GetMenuHandle (HY_CHART_WINDOW_MENU_ID+1);

    DeleteMenu (HY_CHART_WINDOW_MENU_ID+1);
    DisposeMenu (dMenu);

    _HYChartWindow::_UnsetMenuBar();
}

//__________________________________________________________________

bool _HYDistributionChartWindow::_ProcessMenuSelection (long msel)
{
    long        menuChoice = msel&0x0000ffff;

    HiliteMenu(0);
    InvalMenuBar();

    switch (msel/0xffff) {
    case HY_CHART_WINDOW_MENU_ID+1: { // chart menu
        if (menuChoice==1) {
            AddVariable ();
            return true;
        } else if (menuChoice==2) {
            RemoveVariable ();
            return true;
        } else if (menuChoice==3) {
            ShowMarginals ();
            return true;
        } else {
            HandleCatPostProcessor (menuChoice-5);
            return true;
        }

        break;
    }
    }

    return _HYChartWindow::_ProcessMenuSelection(msel);
}

//EOF