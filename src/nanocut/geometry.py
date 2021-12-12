import numpy as np
from nanocut.common import EPSILON
from nanocut.output import error, printstatus

class Geometry:
    """Class for handling crystal structure, containing unit-cell-vectors,
    atom coordinates and names of atoms."""

    def __init__(self, latvecs, basis, basis_names_idx, basis_names,
                 basis_coordsys="lattice", shift=None,
                 shift_coordsys="lattice", bravais_cell=None):
        """Initializes Geometry object.
        
        Args:
            latvecs: Lattice vectors as (3, 3) array.
            basis: Coordinates of the atoms in the basis.
            basis_names_idx: Type index of every atom in the basis.
            basis_names: List of atom types.
            basis_coordsys: Coordinate system for the basis (should be "lattice"
                or "cartesian").
            bravais_cell: Specifies the Bravais superlattice.
        """
        self.latvecs = np.array(latvecs, dtype=float)
        self.basis_names_idx = basis_names_idx
        self.basis_names = basis_names
        self.basis_coordsys = basis_coordsys
        self.basis = np.array(basis, dtype=float)
        self.basis = self.coord_transform(basis, basis_coordsys)
        if shift is not None:
            shift = self.coord_transform(shift, shift_coordsys)
            self.basis += shift
        self.basis = self.mv_basis_to_prim(self.basis)
        if bravais_cell is not None:
            self.bravais_cell = np.array(bravais_cell)
        else:
            self.bravais_cell = np.eye(3, dtype=float)


    @classmethod
    def fromdict(cls, inidict):
        """Initliazes geometry from dict with type checking.
        
        Args:
            cls: Class type.
            inidict: Dictionary with settings
        """
        if "lattice_vectors" not in inidict:
            error("latvecs not defined.")
        if "basis" not in inidict:
            error("Basis not defined.")            
        try:
            latvecs = np.array(
                [ float(s) for s in inidict["lattice_vectors"].split() ])
            latvecs.shape = (3, 3)
        except ValueError:
            error("Invalid lattice vector specification.")
        if abs(np.linalg.det(latvecs)) < EPSILON:
            error("Linearly dependent lattice vectors.")
            
        basis = inidict["basis"].split()
        basis_names=[ basis.pop(idx) 
                      for idx in range(0, len(basis) * 3 // 4, 3) ]
        try:
            basis = np.array([ float(s) for s in basis ])
            basis.shape= (-1, 3)
        except ValueError:
            error("Invalid basis specification.")
        basis_names_idx = range(len(basis))
        basis_coordsys = inidict.get("basis_coordsys", "lattice")
        if basis_coordsys not in ["lattice", "cartesian"]:
            error("Invalid coordinate system specification.")
            
        shiftstr = inidict.get("shift_vector", "0.0 0.0 0.0")
        try:
            shift = np.array([ float(ss) for ss in shiftstr.split() ])
            shift.shape = (3, )
        except ValueError:
            error("Invalid shift vector for basis")
        shift_coordsys = inidict.get("shift_vector_coordsys", "lattice")
        
        bravais_cell = inidict.get("bravais_cell", "1 0 0  0 1 0  0 0 1")
        try:
            bravais_cell = np.array(
                [ float(ss) for ss in bravais_cell.split() ])
            bravais_cell.shape = (3, 3)
        except ValueError:
            error("Invalid Bravais cell definition")

        return cls(latvecs, basis, basis_names_idx, basis_names, basis_coordsys,
                   shift, shift_coordsys, bravais_cell)


    def coord_transform(self, array, array_coordsys):
        """Transforms given vector into cartesian coordinate system.
        
        Args:
            array: Array with coordinates.
            array_coordsys: Coordinate system ("lattice" or "cartesian")
        
        Returns:
            Cartesian coordinates.
        """
        if array_coordsys == "lattice":
            return np.dot(array, self.latvecs)
        elif array_coordsys == "cartesian":
            return array
        else:
            raise ValueError("Invalid coodinate system type '{:s}'".format(
                array_coordsys))


    def mv_basis_to_prim(self, basis):
        """Folds vectors into primitive cell.
        
        Args:
           basis:Coordinates.
        
        Return:
            Coordintes folded into the central cell.
        """
        invlatvecs = np.linalg.inv(self.latvecs)
        basis = np.dot(basis, invlatvecs) % 1.0
        return self.coord_transform(basis, "lattice")


    def gen_cuboid(self, cuboid, periodicity=None):
        """Generates list of lattice points containing all points in a cuboid.
        
        Args:
            cuboid: lower and upper ends of the cuboid.
            periodicity: periodicity object.
            
        Returns:
            Cartesian coordinates of the grid points.
        """

        # Get the 8 corners of the cuboid
        mesh = np.mgrid[0:2,0:2,0:2].reshape(3, -1).transpose()
        abc_corners = [ [ cuboid[mm[0],0], cuboid[mm[1],1], cuboid[mm[2],2] ]
                          for mm in mesh ]

        # Get corners of the containing parallelepiped in relative coordinates
        invlatvecs = np.linalg.inv(self.latvecs)
        nmo_corners = np.dot(abc_corners, invlatvecs)
        nmo_mininds = np.floor(nmo_corners.min(axis=0)).astype(int)
        nmo_maxinds = np.floor(nmo_corners.max(axis=0)).astype(int)
        
        # Generate mesh indexing all points inside the parallelepiped
        nallpoint = np.abs(np.prod(nmo_maxinds - nmo_mininds))
        printstatus(
            "Number of necessary grid points: {:d}".format(int(nallpoint)),
            indentlevel=1)
        nmo = np.mgrid[nmo_mininds[0]:nmo_maxinds[0]+1,
                       nmo_mininds[1]:nmo_maxinds[1]+1,
                       nmo_mininds[2]:nmo_maxinds[2]+1, ]
        nmo = nmo.reshape(3, -1).transpose()
        
        # Throw away points in the parallelepiped which are farther away from
        # the cuboid as the maximal size of the unit cell along given direction
        abc = np.dot(nmo, self.latvecs)
        buffer = np.max(np.abs(self.latvecs), axis=0)
        cond1 = np.all(abc < cuboid[0] - buffer, axis=1)
        cond2 = np.all(abc > cuboid[1] + buffer, axis=1)
        inside = np.logical_not(np.logical_or(cond1, cond2))
        return np.dot(nmo[inside], self.latvecs)
    
    def get_atom_type_names(self):
        return self.basis_names
    

    def get_name_of_atom(self, index):
        """Returns the name of an atom with given index.
        
        Args:
            index: Index of the atom.
            
        Returns:
            Name (type) of the atom.
        """
        return self.basis_names[self.basis_names_idx[index]]


    def gen_atoms(self, lattice_points):
        """Returns the coordinates and index of each atom inside the cells
        corresponding to given lattice points.
        
        Args:
           lattice_points: Cartesian coordinates of the lattice points as
               (-1,3) shaped array.
        
        Returns:
            Coordinates of all atoms in the cells with the given coordinates. 
        """
        nbasis = len(self.basis)
        nlatpoint = len(lattice_points) 
        atoms_coords = np.empty((nlatpoint * nbasis , 3), dtype=float)
        # Taking longer loop (instead over range(len(basis)) to have nicer
        # output with atoms being in the same cell having close indices.
        for ii in range(nlatpoint):
            atoms_coords[ii * nbasis : (ii + 1) * nbasis ] = (self.basis +
                lattice_points[ii])
        atoms_idx = np.resize(np.arange(nbasis), (nlatpoint * nbasis))
        return atoms_coords, atoms_idx
