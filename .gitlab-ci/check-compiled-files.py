###############################################################################
##
## File: check-compiled-files.py
##
## For more information, please see: http://www.nektar.info
##
## The MIT License
##
## Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
## Department of Aeronautics, Imperial College London (UK), and Scientific
## Computing and Imaging Institute, University of Utah (USA).
##
## Permission is hereby granted, free of charge, to any person obtaining a
## copy of this software and associated documentation files (the "Software"),
## to deal in the Software without restriction, including without limitation
## the rights to use, copy, modify, merge, publish, distribute, sublicense,
## and/or sell copies of the Software, and to permit persons to whom the
## Software is furnished to do so, subject to the following conditions:
##
## The above copyright notice and this permission notice shall be included
## in all copies or substantial portions of the Software.
##
## THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
## OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
## FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
## THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
## LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
## FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
## DEALINGS IN THE SOFTWARE.
##
## Description: Check for compiled files using compile_commands.json
##
###############################################################################

import json, sys, glob, os

# Should be run from root source directory.
cwd = os.getcwd()

# Folders and extensions to check
folders = ['library', 'solvers', 'tests', 'utilities']
exts = ['cpp', 'c']

# Presently some stuff is not in the CI.
ignore_sources = [
    # CADfix API
    "library/NekMesh/CADSystem/CFI/CADSystemCFI.cpp",
    "library/NekMesh/CADSystem/CFI/CADSurfCFI.cpp",
    "library/NekMesh/CADSystem/CFI/CADCurveCFI.cpp",
    "library/NekMesh/CADSystem/CFI/CADVertCFI.cpp",
    "library/NekMesh/Module/InputModules/InputCADfix.cpp",
    "library/NekMesh/Module/OutputModules/OutputCADfix.cpp",
    # Likwid
    "solvers/CompressibleFlowSolver/Utilities/TimeRiemann.cpp",
    "solvers/CompressibleFlowSolver/Utilities/TimeRoeKernel.cpp",
    # Template for PWS
    "solvers/PulseWaveSolver/EquationSystems/TemplatePressureArea.cpp",
    # CardiacEPSolver CellMLToNektar template file
    "solvers/CardiacEPSolver/Utilities/CellMLToNektar/nektar/template/model.cpp",
]

ignore_sources = [ os.path.join(cwd, os.path.normpath(p)) for p in ignore_sources ]

with open(sys.argv[1], 'r') as f:
    compilation_data = json.load(f)
    compiled_files = [ entry['file'] for entry in compilation_data ]

    # Search for all c/cpp files in library, solvers and utilities folder.
    found_files = []
    for folder in folders:
        for ext in exts:
            found_files += glob.glob(os.path.join(cwd, folder, '**', '*.{:s}'.format(ext)), recursive=True)

    # Compare the lists of files.
    all_good = True
    for f in found_files:
        if f in ignore_sources:
            continue

        if f not in compiled_files:
            print('Uncompiled file: ' + f)
            all_good = False

    if not all_good:
        exit(1)
