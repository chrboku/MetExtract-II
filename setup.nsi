!define APP_NAME "MetExtractII"
!define INSTDIR_REG_ROOT "HKLM"
!define INSTDIR_REG_KEY "Software\Microsoft\Windows\CurrentVersion\Uninstall\${APP_NAME}"
!include "AdvUninstLog.nsh"
!include "WinVer.nsh"
!include "LogicLib.nsh"

!insertmacro INTERACTIVE_UNINSTALL

Function .onInit
    !insertmacro UNINSTALL.LOG_PREPARE_INSTALL
FunctionEnd
 
Function .onInstSuccess
    !insertmacro UNINSTALL.LOG_UPDATE_INSTALL
FunctionEnd
 
Function UN.onInit
    !insertmacro UNINSTALL.LOG_BEGIN_UNINSTALL
FunctionEnd



RequestExecutionLevel user

# name the installer
OutFile "Setup.exe" 
#InstallDir "$PROGRAMFILES32\MetExtractII\"
InstallDir "$LOCALAPPDATA\MetExtractII_$$METEXTRACTVERSION$$"
ShowInstDetails show
ShowUnInstDetails show 

LicenseText "License information"
LicenseData "LICENSE.txt"


Page license 
Page directory
Page instfiles

# default section start; every NSIS script has at least one section.
Section "instfiles"
    
    # Create installation directory
    CreateDirectory "$INSTDIR"
    
    # Create uninstaller
    WriteUninstaller "$INSTDIR\Uninstall.exe"
    
    # Copy MetExtract II
    SetOutPath $INSTDIR
    
    !insertmacro UNINSTALL.LOG_OPEN_INSTALL
    
    
    File /r "PyMetExtract_$$METEXTRACTVERSION$$\"

    CreateDirectory "$INSTDIR\R"
    
    CreateDirectory "$SMPROGRAMS\MetExtractII $$METEXTRACTVERSION$$"
    CreateShortCut "$SMPROGRAMS\MetExtractII $$METEXTRACTVERSION$$\MetExtract II $$METEXTRACTVERSION$$.lnk" "$INSTDIR\MetExtractII_Main.exe"
    # CreateShortCut "$SMPROGRAMS\MetExtractII\Sample data.lnk" "$INSTDIR\sampleData"
    CreateShortCut "$SMPROGRAMS\MetExtractII $$METEXTRACTVERSION$$\Documentation.lnk" "$INSTDIR\documentation\index.html"
    
    CreateShortCut "$SMPROGRAMS\MetExtractII $$METEXTRACTVERSION$$\Uninstall.lnk" "$INSTDIR\Uninstall.exe"
    
    MessageBox MB_YESNO "Do you want to create a shortcut on the desktop?" /SD IDYES IDNO endDesktopIcon
                
        CreateShortCut "$DESKTOP\MetExtract II $$METEXTRACTVERSION$$.lnk" "$INSTDIR\MetExtractII_Main.exe"
            
    endDesktopIcon:
    
    
    Call installR
    Call configureR
    Call checkWinXPProblem
    
    !insertmacro UNINSTALL.LOG_CLOSE_INSTALL
# default section end
SectionEnd


# uninstaller section start
Section "uninstall"

    #Verify the uninstaller - last chance to back out
    MessageBox MB_OKCANCEL "Do you want to permanently remove ${APP_NAME}?" IDOK next
        Abort
        Quit
    next:
          
    RmDir /r "$INSTDIR\sampleData"
    RmDir /r "$INSTDIR\R"
    RmDir /r "$INSTDIR\documentation"
    RmDir /r "$INSTDIR\mpl-data"
    RmDir /r "$INSTDIR\Settings"
    RmDir /r "$INSTDIR\tcl"
          
    !insertmacro UNINSTALL.LOG_UNINSTALL "$INSTDIR"
    !insertmacro UNINSTALL.LOG_UNINSTALL "$APPDATA\${APP_NAME}"
    !insertmacro UNINSTALL.LOG_END_UNINSTALL
    
    DeleteRegKey /ifempty ${INSTDIR_REG_ROOT} "${INSTDIR_REG_KEY}"
  
    # second, remove the links from the start menu
    RmDir /r "$SMPROGRAMS\MetExtractII $$METEXTRACTVERSION$$"

    delete "$DESKTOP\MetExtract II $$METEXTRACTVERSION$$.lnk"
 
# uninstaller section end
SectionEnd


Function installR
    MessageBox MB_YESNO "Install R? (3.3.2; This will not affect other installation/versions of R currently installed and no registry or environment variables will be created)" /SD IDYES IDNO endinstallR
        File "resources\setupR_minimal.inf"
        #File "resources\R-2.15.2-win.exe"

        NSISdl::download "https://cran.r-project.org/bin/windows/base/old/3.3.2/R-3.3.2-win.exe" "$INSTDIR\R-3.3.2-win.exe"
        Pop $0
        ${If} $0 == "success"
            ExecWait "R-3.3.2-win.exe /SILENT /DIR=$\"$INSTDIR\R$\" /LOADINF=$\"$INSTDIR\setupR_minimal.inf$\""
        ${Else}
            MessageBox mb_iconstop "Error downloading R 2.15.2. Please download it from https://cran.r-project.org/bin/windows/base/old/3.3.2/R-3.3.2-win.exe"
        ${EndIf}

    endinstallR:
FunctionEnd

Function configureR
    MessageBox MB_YESNO "Configure R installation? This will not affect other installation/versions of R currently installed and no registry or environment variables will be created." /SD IDYES IDNO showNonConfigMessage
        CopyFiles "$INSTDIR\R\bin\i386\*.*" "$INSTDIR\R\bin"
        
        FileOpen $9 RPATH.conf w 
        FileWrite $9 "$INSTDIR\R"
        FileClose $9 
        
        Goto endConfigureR
            
    showNonConfigMessage:
        MessageBox MB_OK "Please follow instructions in the documentation or after program start to configure your R-instance."
    endConfigureR:
FunctionEnd

Function checkWinXPProblem
    ${If} ${IsWinXP}
    ${AndIfNot} ${AtLeastServicePack} 1
        MessageBox MB_YESNO "WinXP without SP1 requires the Micorosoft C++ Redistributables 2008. Do you want to download them?" /SD IDYES IDNO endWinXPTest
            ExecShell "open" "http://www.microsoft.com/en-us/download/details.aspx?id=29"
        endWinXPTest:
    ${EndIf}
FunctionEnd
