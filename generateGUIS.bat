@echo off
CALL c:\Python27\Lib\site-packages\PyQt4\pyuic4.bat -x .\mePyGuis\guis\mainwindow.ui                  -o .\mePyGuis\mainWindow.py
CALL c:\Python27\Lib\site-packages\PyQt4\pyuic4.bat -x .\mePyGuis\guis\TracerEditor.ui                -o .\mePyGuis\TracerEditor.py
CALL c:\Python27\Lib\site-packages\PyQt4\pyuic4.bat -x .\mePyGuis\guis\adductsEditor.ui               -o .\mePyGuis\adductsEditor.py
CALL c:\Python27\Lib\site-packages\PyQt4\pyuic4.bat -x .\mePyGuis\guis\groupEditor.ui                 -o .\mePyGuis\groupEditor.py
CALL c:\Python27\Lib\site-packages\PyQt4\pyuic4.bat -x .\mePyGuis\guis\heteroAtomEditor.ui            -o .\mePyGuis\heteroAtomEditor.py
CALL c:\Python27\Lib\site-packages\PyQt4\pyuic4.bat -x .\mePyGuis\guis\TSVLoaderEditor.ui             -o .\mePyGuis\TSVLoaderEditor.py
CALL c:\Python27\Lib\site-packages\PyQt4\pyuic4.bat -x .\mePyGuis\guis\ModuleSelectionWindow.ui       -o .\mePyGuis\ModuleSelectionWindow.py
CALL c:\Python27\Lib\site-packages\PyQt4\pyuic4.bat -x .\mePyGuis\guis\calcIsotopeEnrichmentDialog.ui -o .\mePyGuis\calcIsotopeEnrichmentDialog.py
CALL c:\Python27\Lib\site-packages\PyQt4\pyuic4.bat -x .\mePyGuis\guis\FE_mainWindow.ui               -o .\mePyGuis\FE_mainWindow.py
CALL c:\Python27\Lib\site-packages\PyQt4\pyuic4.bat -x .\mePyGuis\guis\FTICRwindow.ui                 -o .\mePyGuis\FTICRWindow.py
CALL c:\Python27\Lib\site-packages\PyQt4\pyuic4.bat -x .\mePyGuis\guis\combineResultsDialog.ui        -o .\mePyGuis\combineResultsDialog.py
CALL c:\Python27\Lib\site-packages\PyQt4\pyrcc4 -o resources_rc.py resources.qrc
echo "Guis created.."

::CALL c:\Python26\Lib\site-packages\PyQt4\pyuic4.bat -x guis\SettingsWizard.ui -o .\guis\pys\SettingsWizard.py