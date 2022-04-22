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
from termcolor import colored

from MetExtractII_Main import MetExtractVersion

import matplotlib

import sys
sys.path.append("../PyMassBankSearchTool")
from TableUtils import TableUtils





import os, stat
def on_rm_error( func, path, exc_info):
    # path contains the path of the file that couldn't be removed
    # let's just assume that it's read-only and unlink it.
    os.chmod( path, stat.S_IWRITE )
    os.unlink( path )



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

os.environ["R_USER"]=get_main_dir()+"/Ruser"
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
    print colored("Error: could not clean up py2exe environment prior to compilation\n==============================\n", "red")

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
import openpyxl


class Target:
    def __init__(self, **kw):
        self.__dict__.update(kw)

a=Target(script = "MetExtractII_Main.py")
b=Target(script = "FragExtract.py")
c=Target(script = "MExtract.py")
d=Target(script = "resultsPostProcessing/combineResults.py")
e=Target(script = "FTICRModule.py")

print "###################################################"
print "########## Packing MetExtractII_Main"
print "###################################################"
import sys
sys.setrecursionlimit(5000)
setup(console=[a,b,c,d,e],
      options={"py2exe": {
                 "includes": ["sip", "matplotlib.backends.backend_tkagg", 'scipy', 'scipy.integrate', 'scipy.special.*','scipy.linalg.*', 'scipy.sparse.csgraph._validation', 'scipy._lib.messagestream'],  # use this line if above does not work
                 "dll_excludes": ["MSVCP90.dll", "api-ms-win-core-string-l1-1-0.dll","api-ms-win-core-registry-l1-1-0.dll","api-ms-win-core-errorhandling-l1-1-0.dll","api-ms-win-core-string-l2-1-0.dll",
                                  "api-ms-win-core-profile-l1-1-0.dll","api-ms-win*.dll","api-ms-win-core-processthreads-l1-1-2.dll","api-ms-win-core-libraryloader-l1-2-1.dll","api-ms-win-core-file-l1-2-1.dll",
                                  "api-ms-win-security-base-l1-2-0.dll","api-ms-win-eventing-provider-l1-1-0.dll","api-ms-win-core-heap-l2-1-0.dll","api-ms-win-core-libraryloader-l1-2-0.dll","api-ms-win-core-localization-l1-2-1.dll",
                                  "api-ms-win-core-sysinfo-l1-1-0.dll","api-ms-win-core-synch-l1-2-0.dll","api-ms-win-core-heap-l1-2-0.dll","api-ms-win-core-handle-l1-1-0.dll","api-ms-win-core-io-l1-1-1.dll","api-ms-win-core-com-l1-1-1.dll",
                                  "api-ms-win-core-memory-l1-1-2.dll","api-ms-win-core-version-l1-1-1.dll","api-ms-win-core-version-l1-1-0.dll","api-ms-win-core-processthreads-l1-1-0.dll"],
                 "excludes": ["_gtkagg", "_tkagg", 'jinja2.asyncsupport','jinja2.asyncfilters'],
                 "packages": ["FileDialog", "openpyxl", 'reportlab','reportlab.graphics.charts','reportlab.graphics.samples','reportlab.graphics.widgets','reportlab.graphics.barcode','reportlab.graphics','reportlab.lib','reportlab.pdfbase','reportlab.pdfgen','reportlab.platypus', 'zeep', 'lxml'],
                 'dist_dir': "./dist"
      }},
      data_files=data_files,
      requires=['matplotlib'])
#rmtree("./build/")


print "forced waiting..."
import time
time.sleep(3)


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
    print "Additional resources copied\n==============================\n"
except:
    print colored("Error: Could not copy all required files", "red")
    err = True
    import sys
    sys.exit(1)

try:
    import shutil
    dest="./dist/documentation"
    src="./documentation"
    shutil.copytree(src, dest)

    print "Help files copied\n==============================\n"

except Exception as ex:
    print colored("Error: Could not copy help files: "+ex.message, "red")
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
    print colored("Error: Could not rename dist folder", "red")
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
    sleep(3)
    rmtree("./build", onerror = on_rm_error)
except:
    print colored(
        "Cleanup failed (./build) dist and/or build directories still there\n==============================\n", "red")
    import traceback

    traceback.print_exc()

try:
    sleep(3)
    rmtree(meDistFolder, onerror = on_rm_error)
except:
    print colored(
        "Cleanup failed (" + meDistFolder + ") dist and/or build directories still there\n==============================\n",
        "red")
    import traceback

    traceback.print_exc()