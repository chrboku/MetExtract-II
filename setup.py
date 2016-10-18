#    MetExtract II
#    Copyright (C) 2015
#
#    This program is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation; either version 2 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program; if not, write to the Free Software
#    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA



from distutils.core import setup
import py2exe
from shutil import rmtree, copy, move
import os
import sys
import zipfile

from MetExtractII_Main import MetExtractVersion

import matplotlib

import sys
sys.path.append("C:/Users/cbueschl/Documents/Dropbox/IFA Tulln/Archive/PyMetExtract/PyMetExtract")
from TableUtils import TableUtils








def replaceInFile(textToReplace, newText, filein, fileout=None):
    f = open(filein,'r')
    filedata = f.read()
    f.close()

    newdata = filedata.replace(textToReplace, newText)

    if fileout is None:
        fileout=filein
    f = open(fileout,'w')
    f.write(newdata)
    f.close()


########################################################################################################################
########################################################################################################################
########################################################################################################################



#<editor-fold desc="### check if R is installed and accessible">
def checkR():
    try:
        import rpy2.robjects as ro              # import RPy2 module
        r = ro.r                                # make R globally accessible

        v = r("R.Version()$version.string")     # if this R-commando is executed, the RPy2 connection to the
                                                # R subprocess has been established
        return True
    except:
        # The R subprocess could not be started / accessed successfully
        return False

def loadRConfFile(path):
    import os
    if os.path.isfile(path+"/RPATH.conf"):
        with open(path+"/RPATH.conf", "rb") as rconf:
            line=rconf.readline()
            os.environ["R_HOME"]=line
            return True
    else:
        return False


__RHOMEENVVAR=""
import os
from utils import get_main_dir
if "R_HOME" in os.environ.keys():
    __RHOMEENVVAR=os.environ["R_HOME"]


# try to load r configuration file (does not require any environment variables or registry keys)
if not loadRConfFile(path=get_main_dir()) or not checkR():
    os.environ["R_HOME"]=get_main_dir()+"/R"

    if checkR():
        with open("RPATH.conf", "wb") as rconf:
            rconf.write(get_main_dir()+"/R")
            tryLoad=False

            # Show a dialog box to the user that R could not be started
            from os import sys
            from PyQt4 import QtGui, QtCore

            app = QtGui.QApplication(sys.argv)

            QtGui.QMessageBox.information(None, "MetExtract",
                      "R successfully configured\nUsing MetExtract R-Installation\nPlease restart",
                      QtGui.QMessageBox.Ok)
            sys.exit(0)
    else:

        os.environ["R_HOME"]=__RHOMEENVVAR
        os.environ["R_HOME_FROM"]="RPATH environment variable"
        if not checkR():

            print "Error: R could not be loaded correctly (No RPATH.conf file or R_HOME environment variable found)\nPlease make sure it is installed and accessible"

            # Show a dialog box to the user that R could not be started
            from os import sys
            from PyQt4 import QtGui, QtCore

            app = QtGui.QApplication(sys.argv)

            if QtGui.QMessageBox.warning(None, "MetExtract",
                                      "Error: R could not be loaded\nPlease make sure it is installed and accessible\n"
                                      "The default installation path is C:\\Program Files\\R\n"
                                      "Do you want to specify the folder?",
                                      QtGui.QMessageBox.Yes | QtGui.QMessageBox.No)==QtGui.QMessageBox.Yes:
                tryLoad=True
                from utils import get_main_dir
                lastDir=get_main_dir()
                while tryLoad:
                    folder = str(QtGui.QFileDialog.getExistingDirectory(None, "Select R-directory (not bin folder)", directory=lastDir))
                    if folder=="":
                        sys.exit(1)
                    else:
                        lastDir=folder
                        os.environ["R_HOME"]=folder
                        if checkR():
                            with open("RPATH.conf", "wb") as rconf:
                                rconf.write(folder)
                                tryLoad=False

                                QtGui.QMessageBox.information(None, "MetExtract",
                                          "R successfully configured\nPlease restart",
                                          QtGui.QMessageBox.Ok)
                                sys.exit(0)
                        else:
                            if QtGui.QMessageBox.warning(None, "MetExtract",
                                          "Error: R could not be loaded from the specified location\n"
                                          "%s\n\n"
                                          "Please make sure it is installed and accessible\n"
                                          "The default installation path is C:\\Program Files\\R\n"
                                          "Do you want to specify the folder?"%folder,
                                          QtGui.QMessageBox.Yes | QtGui.QMessageBox.No)==QtGui.QMessageBox.Yes:
                                pass
                            else:
                                sys.exit(1)
            else:
                sys.exit(1)
else:
    os.environ["R_HOME_FROM"]="RPATH.conf of MetExtract II"
#</editor-fold>
#Used to locate R.dll (no idea why this is necessary)
import rpy2.robjects as ro
r = ro.r




# remove previous setup files
try:
    rmtree("./dist/")
except:
    pass
try:
    rmtree("./build/")
except:
    pass
try:
    rmtree("./distribute/")
except:
    print "Error: could not clean up py2exe environment prior to compilation\n==============================\n"

# get local files (images, R-Scripts, ...)
data_files = matplotlib.get_py2exe_datafiles()
data_files.append("./chromPeakPicking/MassSpecWaveletIdentification.r")
data_files.append("./XICAlignment.r")

err = False

import os
import shutil

def mergeDirs (root_src_dir, root_dst_dir):
    for src_dir, dirs, files in os.walk(root_src_dir):
        dst_dir = src_dir.replace(root_src_dir, root_dst_dir)
        if not os.path.exists(dst_dir):
            os.mkdir(dst_dir)
        for file_ in files:
            src_file = os.path.join(src_dir, file_)
            dst_file = os.path.join(dst_dir, file_)
            if os.path.exists(dst_file):
                os.remove(dst_file)
            shutil.move(src_file, dst_dir)

# create executable compilation
print "###################################################"
print "########## Packing MetExtractII_Main"
print "###################################################"
setup(console=[{"script": "MetExtractII_Main.py"}],
      options={"py2exe": {
                 "includes": ["sip", "matplotlib.backends.backend_tkagg"],  # use this line if above does not work
                 "dll_excludes": ["MSVCP90.dll"],
                 "excludes": ["_gtkagg", "_tkagg"]
      }})
shutil.copytree("./dist", "./dist_MetExtract_Main")
#rmtree("./build/")

print "###################################################"
print "########## Packing FragExtract"
print "###################################################"
setup(console=[{"script": "FragExtract.py"}],
      options={"py2exe": {
                 "includes": ["sip", "matplotlib.backends.backend_tkagg"],
                 "dll_excludes": ["MSVCP90.dll"],
                 "excludes": ["_gtkagg", "_tkagg"]}},
      data_files=data_files,
      requires=['matplotlib'])
shutil.copytree("./dist", "./dist_FragExtract")
#rmtree("./build/")

print "###################################################"
print "########## Packing MExtract"
print "###################################################"
setup(console=[{"script": "MExtract.py"}],
      options={"py2exe": {
                 "includes": ["sip", "matplotlib.backends.backend_tkagg"],
                 "dll_excludes": ["MSVCP90.dll"],
                 "excludes": ["_gtkagg", "_tkagg"]}},
      data_files=data_files,
      requires=['matplotlib'])
shutil.copytree("./dist", "./dist_MExtract")

print "forced waiting..."
import time
time.sleep(3)

rmtree("./dist/")
rmtree("./build/")

os.mkdir("./dist")

mergeDirs("./dist_MetExtract_Main", "./dist")
mergeDirs("./dist_MExtract", "./dist")
mergeDirs("./dist_FragExtract", "./dist")

rmtree("./dist_MetExtract_Main/")
rmtree("./dist_MExtract/")
rmtree("./dist_FragExtract/")

print "Setup finished\n==============================\n"

# copy Settings and necessary files to distribution folder
try:
    os.makedirs("./dist/Settings/")
    copy("./Settings/defaultSettings.ini", "./dist/Settings/defaultSettings.ini")
    copy("./Settings/LTQ-Orbitrap-XL__HPLC.ini", "./dist/Settings/LTQ-Orbitrap-XL__HPLC.ini")
    copy("./Settings/QExactive__HPLC.ini", "./dist/Settings/QExactive__HPLC.ini")
    copy("./Settings/QTof__HPLC.ini", "./dist/Settings/QTof__HPLC.ini")
    os.makedirs("./dist/Settings/Tracers/")
    copy("./Settings/Tracers/DON_12C13C.ini", "./dist/Settings/Tracers/DON_12C13C.ini")

    os.makedirs("./dist/chromPeakPicking/")
    copy("./chromPeakPicking/MassSpecWaveletIdentification.r", "./dist/chromPeakPicking/MassSpecWaveletIdentification.r")
    copy("./XICAlignment.r", "./dist/XICAlignment.r")
    copy("./LICENSE.txt", "./dist/LICENSE.txt")
    copy("./calculateIsotopeEnrichment.R", "./dist/calculateIsotopeEnrichment.R")
    print "Additional resources copied\n==============================\n"
except:
    print "Error: Could not copy all required files"
    err = True


# copy documentation
def copyAllFilesInFolder(src, dest):
    src_files = os.listdir(src)
    for file_name in src_files:
        full_file_name = os.path.join(src, file_name)
        if (os.path.isfile(full_file_name)):
            shutil.copy(full_file_name, dest)
try:
    import shutil
    dest="./dist/documentation"
    os.makedirs(dest)
    src="./documentation"

    copyAllFilesInFolder(src, dest)
    dest="./dist/documentation/figures/"
    os.makedirs(dest)
    src="./documentation/figures"
    copyAllFilesInFolder(src, dest)

    copy("./help/INSTALL.txt", "./dist/INSTALL.txt")

    if False:
        os.makedirs("./dist/sampleData")
        os.makedirs("./dist/sampleData/LTQ_Orbitrap_XL")
        os.makedirs("./dist/sampleData/Exactive_plus")

        copy("C:/Users/cbueschl/Desktop/implTest/LTQ_Orbitrap_XL/10_Remus_DON.mzXML", "./dist/sampleData/LTQ_Orbitrap_XL/10_Remus_DON.mzXML")
        copy("C:/Users/cbueschl/Desktop/implTest/LTQ_Orbitrap_XL/10_CM_DON.mzXML", "./dist/sampleData/LTQ_Orbitrap_XL/10_CM_DON.mzXML")
        copy("C:/Users/cbueschl/Desktop/implTest/Orbitrap.grp", "./dist/sampleData/DON_in_Wheat.grp")
        copy("C:/Users/cbueschl/Desktop/implTest/Orbitrap_LAC.grp", "./dist/sampleData/DON_in_Wheat_LAC.grp")

        replaceInFile("python \"%METEXTRACTLOCATION%/MExtract.py\"", "\"%METEXTRACTLOCATION%/MExtract.exe\"",
                      filein="C:/Users/cbueschl/Desktop/implTest/_processExperiment_LTQOrbitrapXL.bat",
                      fileout="./dist/sampleData/_processExperiment_LTQOrbitrapXL.bat")

        copy("C:/Users/cbueschl/Desktop/implTest/Exactive_plus/WheatEar_DON_posneg.mzXML", "./dist/sampleData/Exactive_plus/WheatEar_DON_posneg.mzXML")
        copy("C:/Users/cbueschl/Desktop/implTest/Exactive_posneg.grp", "./dist/sampleData/Exactive_posneg.grp")

        copy("E:/140124_pos_131_MSMS_FragExtract/Fg_Extract_ACN_CID25.mzXML", "./dist/sampleData/LTQ_Orbitrap_XL/Fg_Extract_ACN_CID25.mzXML")
        copy("E:/140124_pos_131_MSMS_FragExtract/Fg_Extract_ACN_CID35.mzXML", "./dist/sampleData/LTQ_Orbitrap_XL/Fg_Extract_ACN_CID35.mzXML")
        copy("E:/140124_pos_131_MSMS_FragExtract/FE_groups_small.grp", "./dist/sampleData/FragExtract.grp")

    print "Help files copied\n==============================\n"

except:
    print "Error: Could not copy help files"
    err = True
    import sys
    sys.exit(1)

# rename dist folder to PyMetExtract and current version of the software
meDistFolder="./dist"
try:
    from time import sleep  # sometimes, the re-naming does not work (probably some kind of lock from the OS)
    sleep(3)                # this short waiting time decreases the number of times the renaming does not work
    os.rename("./dist", "./PyMetExtract_%s"%MetExtractVersion)
    meDistFolder="./PyMetExtract_%s"%MetExtractVersion
    print "Distribution renamed\n==============================\n"

except:
    print "Error: Could not rename dist folder"
    err = True

try:
    os.makedirs('./distribute')
except:
    pass

import os

## update MetExtract Version in NSIS setup file
replaceInFile("$$METEXTRACTVERSION$$", MetExtractVersion, filein="setup.nsi", fileout="setup_curVersion.nsi")

os.system("\"c:\\Program Files (x86)\\NSIS\\makensis.exe\" setup_curVersion.nsi")
os.remove("./setup_curVersion.nsi")
os.rename("./Setup.exe", "./Setup_MetExtractII_%s.exe"%MetExtractVersion)
move("./Setup_MetExtractII_%s.exe"%MetExtractVersion, "./distribute/Setup_MetExtractII_%s.exe"%MetExtractVersion)
#shutil.copy("../../distribution/vcredist_x86.exe", "./distribute/vcredist_x86.exe")


def zipdir(path, zip):
    for root, dirs, files in os.walk(path):
        for file in files:
            zip.write(os.path.join(root, file))

# create zip archive from the executables and documentation
if not err:

    print "Zipping MetExtract (%s)" % MetExtractVersion
    zipFileName = "PyMetExtract_%s.zip" % MetExtractVersion

    zipF = zipfile.ZipFile(zipFileName, 'w')
    zipdir(meDistFolder, zipF)
    zipF.close()

    move('./%s' % zipFileName, './distribute/%s' % zipFileName)

    print "MetExtract (%s) created\n see %s\n==============================\n" % (
        MetExtractVersion, './distribute/%s' % zipFileName)

    try:
        rmtree(meDistFolder)
    except:
        print "Cleanup failed. dist and/or build directories still there\n==============================\n"


os.makedirs("./BootstrapKnitr_Template")
os.makedirs("./BootstrapKnitr_Template/dataIn")
os.makedirs("./BootstrapKnitr_Template/dataOut")
os.makedirs("./BootstrapKnitr_Template/figure")
os.makedirs("./BootstrapKnitr_Template/documentation")
os.makedirs("./BootstrapKnitr_Template/documentation/figures")
os.makedirs("./BootstrapKnitr_Template/scripts")
os.makedirs("./BootstrapKnitr_Template/scripts/templateScripts")
os.makedirs("./BootstrapKnitr_Template/scripts/extLib")
os.makedirs("./BootstrapKnitr_Template/scripts/js")

copyAllFilesInFolder("./../BootstrapKnitr_Template", "./BootstrapKnitr_Template")
copyAllFilesInFolder("./../BootstrapKnitr_Template/dataIn", "./BootstrapKnitr_Template/dataIn")
copyAllFilesInFolder("./../BootstrapKnitr_Template/dataOut", "./BootstrapKnitr_Template/dataOut")
copyAllFilesInFolder("./../BootstrapKnitr_Template/figure", "./BootstrapKnitr_Template/figure")
copyAllFilesInFolder("./../BootstrapKnitr_Template/documentation", "./BootstrapKnitr_Template/documentation")
copyAllFilesInFolder("./../BootstrapKnitr_Template/documentation/figures", "./BootstrapKnitr_Template/documentation/figures")
copyAllFilesInFolder("./../BootstrapKnitr_Template/documentation/importantPackagesWithVersion", "./BootstrapKnitr_Template/documentation/importantPackagesWithVersion")
copyAllFilesInFolder("./../BootstrapKnitr_Template/scripts/templateScripts", "./BootstrapKnitr_Template/scripts/templateScripts")
copyAllFilesInFolder("./../BootstrapKnitr_Template/scripts/extLib", "./BootstrapKnitr_Template/scripts/extLib")
copyAllFilesInFolder("./../BootstrapKnitr_Template/scripts/js", "./BootstrapKnitr_Template/scripts/js")

zipF = zipfile.ZipFile("./distribute/BootstrapKnitr_Template.zip", 'w')
zipdir("./BootstrapKnitr_Template", zipF)
zipF.close()
rmtree("./BootstrapKnitr_Template")

print "BootstrapKnitr_Template zipped"

copy("./../distribution/R-2.15.2-win.exe", "./distribute/R-2.15.2-win.exe")

print "R 2.15.2 copied"

copy("./../distribution/RStudio-0.97.551.exe", "./distribute/RStudio-0.97.551.exe")

print "RStudio copied"

copy ("./calculateIsotopeEnrichment.R", "./distribute/calcIsotopicEnrichment.R")

print "calcIsotopicEnrichment.R copied"

#rmtree("./dist/")
#rmtree("./dist/")