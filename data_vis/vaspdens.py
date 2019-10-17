from ase.calculators.vasp import VaspChargeDensity
import numpy as np

def read_vasp_density(density_file):
    data = VaspChargeDensity(density_file)
    density = data.chg[-1]
    atoms = data.atoms[-1]
    ngridpts=np.array(density.shape)
    totalgridpts = ngridpts.prod()
    unitcell = atoms.get_cell()
    return data,density,atoms,ngridpts,totalgridpts,unitcell

class VaspDens:
    
    def __init__(self, density, atoms, grid_points, cell):
        self.density = density
        self.atoms = atoms, 
        self.grid_points = grid_points
        self.cell = cell
        
    def from_file(file_path):
        data = VaspChargeDensity(file_path)
        density = data.chg[-1]
        atoms = data.atoms[-1]
        ngridpts=np.array(density.shape)
        totalgridpts = ngridpts.prod()
        unitcell = atoms.get_cell()
        return VaspDens(density, atoms, ngridpts, unitcell)
    
    @property
    def ngridpoints(self):
        n_grid_points = self.grid_points.prod()
        return n_grid_points
    
    @property
    def cell_lengths(self):
        cell_lengths = np.sqrt(np.dot(self.cell,self.cell.transpose().diagonal()))
        return cell_lengths
    
    def _get_plane(self, lattice_vector_not_included, distance):
        if lattice_vector_not_included == 'a':
            plane, vector = [1,2], 0
        elif lattice_vector_not_included == 'b':
            plane, vector = [0,2], 1
        elif lattice_vector_not_included == 'c':
            plane, vector = [0,1], 2
        else:
            raise ValueError('lattice vector must be a, b, or c')
        index = int(round(self.grid_points[vector]*distance/self.cell_lengths[vector]))%self.grid_points[vector]
        return plane,index,vector
    
    def get_density_slice_for_plot(self, slice_in, distance):
        plane,index,vector = self._get_plane(slice_in, distance)
        if vector == 0:
            density2D = self.density[index,:,:]
        if vector == 1:
            density2D = self.density[:,index,:]
        if vector == 2:
            density2D = self.density[:,:,index]
        yontox = np.dot(self.cell[plane[0]],self.cell[plane[1]].T)/self.cell_lengths[0]
        ynormal = np.cross(self.cell[plane[0]],self.cell[plane[1]].T)/self.cell_lengths[0]
        ynormal=np.sqrt(np.dot(ynormal,ynormal.T))
        xarray=np.zeros((self.grid_points[plane[0]],self.grid_points[plane[1]]),np.float)
        yarray=np.zeros((self.grid_points[plane[0]],self.grid_points[plane[1]]),np.float)
        for i in range(self.grid_points[plane[0]]):
            for j in range(self.grid_points[plane[1]]):
                xarray[i][j]=float(i)/float(self.grid_points[plane[0]])*self.cell_lengths[plane[0]]+float(j)/float(self.grid_points[plane[1]])*yontox
                yarray[i][j]=float(j)/float(self.grid_points[plane[1]])*ynormal
        return xarray,yarray,density2D 
            
        
    
    
