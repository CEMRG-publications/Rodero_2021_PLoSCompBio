import numpy as np
import os
import vtk
from vtk.numpy_interface import dataset_adapter as dsa

def extract_UVCs(input_dir, input_file, output_dir):
    """
    Function to extract the UVCs from an ASCII vtk.
    """
    reader = vtk.vtkUnstructuredGridReader()
    reader.SetFileName(os.path.join(input_dir, input_file))
    reader.ReadAllScalarsOn()
    reader.Update()
    usg = dsa.WrapDataObject(reader.GetOutput())

    rho = usg.PointData['RHO.dat']
    phi = usg.PointData['PHI.dat']
    V = usg.PointData['V.dat']
    Z = usg.PointData['Z.dat']

    np.savetxt(os.path.join(output_dir, "COORDS_RHO.dat"), rho, fmt='%4.6f')
    np.savetxt(os.path.join(output_dir, "COORDS_PHI.dat"), phi, fmt='%4.6f')
    np.savetxt(os.path.join(output_dir, "COORDS_V.dat"), V, fmt='%4.6f')
    np.savetxt(os.path.join(output_dir, "COORDS_Z.dat"), Z, fmt='%4.6f')