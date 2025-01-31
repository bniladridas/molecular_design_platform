import logging
from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse
from pydantic import BaseModel
from typing import List, Dict

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, Draw
import numpy as np
import base64
import io

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

app = FastAPI(
    title="Molecular Design Platform",
    description="Advanced molecular generation and analysis system",
    version="0.1.0"
)

# CORS Configuration
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],  # Allows all origins
    allow_credentials=True,
    allow_methods=["*"],  # Allows all methods
    allow_headers=["*"],  # Allows all headers
)

class MoleculeGenerationRequest(BaseModel):
    base_smiles: str
    num_mutations: int = 3
    mutation_types: List[str] = ["substituent", "bond_order", "atom_swap"]

class MoleculeResponse(BaseModel):
    original_smiles: str
    generated_smiles: str
    molecular_image: str
    properties: Dict[str, float]

class MoleculeGenerator:
    @staticmethod
    def generate_random_molecule(base_smiles: str, num_mutations: int = 3) -> Chem.Mol:
        base_mol = Chem.MolFromSmiles(base_smiles)
        
        if base_mol is None:
            raise ValueError(f"Invalid SMILES: {base_smiles}")
        
        mol = Chem.MolFromSmiles(base_smiles)
        
        for _ in range(num_mutations):
            mol = MoleculeGenerator._random_mutation(mol)
        
        return mol
    
    @staticmethod
    def _random_mutation(mol):
        mutation_types = [
            MoleculeGenerator._add_random_substituent,
            MoleculeGenerator._modify_bond_order,
            MoleculeGenerator._swap_atom
        ]
        
        mutation_func = np.random.choice(mutation_types)
        return mutation_func(mol)
    
    @staticmethod
    def _add_random_substituent(mol):
        # Ensure input is a valid molecule
        if mol is None:
            return None
        
        try:
            # Possible substituents
            substituents = ['C', 'CC', 'N', 'O', 'Cl', 'Br']
            
            # Convert to SMILES to work with string representation
            mol_smiles = Chem.MolToSmiles(mol, isomericSmiles=True)
            
            # Find atoms where substituents can be added
            attachable_atoms = [
                i for i, atom in enumerate(mol.GetAtoms()) 
                if atom.GetDegree() < 4  # Limit to atoms with fewer than 4 bonds
            ]
            
            # If no attachable atoms, return original molecule
            if not attachable_atoms:
                return mol
            
            # Randomly select an atom to attach substituent
            attach_atom_idx = int(np.random.choice(attachable_atoms))
            substituent = np.random.choice(substituents)
            
            # Create modified SMILES by adding substituent
            modified_smiles = mol_smiles[:attach_atom_idx] + f'({substituent})' + mol_smiles[attach_atom_idx:]
            
            # Convert back to molecule
            mol_result = Chem.MolFromSmiles(modified_smiles)
            
            return mol_result if mol_result is not None else mol
        except Exception as e:
            print(f"Error in _add_random_substituent: {e}")
            return mol
    
    @staticmethod
    def _modify_bond_order(mol):
        # Ensure input is a valid molecule
        if mol is None:
            return None
        
        try:
            # Convert mol to a canonical SMILES 
            mol_smiles = Chem.MolToSmiles(mol, isomericSmiles=True)
            mol = Chem.MolFromSmiles(mol_smiles)
            
            # Get all bonds
            bonds = mol.GetBonds()
            if not bonds:
                return mol
            
            # Randomly select a bond
            bond_to_modify = np.random.choice(bonds)
            
            # Prepare bond order options
            bond_orders = [
                Chem.BondType.SINGLE,
                Chem.BondType.DOUBLE,
                Chem.BondType.TRIPLE
            ]
            
            # Remove the current bond order
            current_bond_type = bond_to_modify.GetBondType()
            bond_orders = [bo for bo in bond_orders if bo != current_bond_type]
            
            # Select a new bond order
            new_bond_type = np.random.choice(bond_orders)
            
            # Create an editable molecule
            rw_mol = Chem.RWMol(mol)
            
            # Modify the bond
            start_atom = int(bond_to_modify.GetBeginAtomIdx())
            end_atom = int(bond_to_modify.GetEndAtomIdx())
            bond_idx = int(bond_to_modify.GetIdx())
            
            # Remove the existing bond
            rw_mol.RemoveBond(start_atom, end_atom)
            
            # Add a new bond with the modified order
            rw_mol.AddBond(start_atom, end_atom, new_bond_type)
            
            # Convert back to immutable molecule
            mol_result = Chem.MolFromSmiles(Chem.MolToSmiles(rw_mol))
            
            return mol_result if mol_result is not None else mol
        except Exception as e:
            print(f"Error in _modify_bond_order: {e}")
            return mol

    @staticmethod
    def _swap_atom(mol):
        # Ensure input is a valid molecule
        if mol is None:
            return None
        
        try:
            # Get atom types to swap
            atom_types = ['C', 'N', 'O', 'S', 'P']
            
            # Convert to SMILES to work with string representation
            mol_smiles = Chem.MolToSmiles(mol, isomericSmiles=True)
            
            # Find atoms that can be swapped
            swappable_atoms = [
                i for i, atom in enumerate(mol.GetAtoms()) 
                if atom.GetSymbol() in atom_types
            ]
            
            # If no swappable atoms, return original molecule
            if not swappable_atoms:
                return mol
            
            # Randomly select an atom to swap
            atom_to_swap_idx = int(np.random.choice(swappable_atoms))
            current_atom = mol.GetAtomWithIdx(atom_to_swap_idx).GetSymbol()
            
            # Get possible swap atoms (excluding current atom)
            swap_options = [at for at in atom_types if at != current_atom]
            new_atom = np.random.choice(swap_options)
            
            # Create modified SMILES by replacing atom
            modified_smiles = list(mol_smiles)
            modified_smiles[atom_to_swap_idx] = new_atom
            modified_smiles = ''.join(modified_smiles)
            
            # Convert back to molecule
            mol_result = Chem.MolFromSmiles(modified_smiles)
            
            return mol_result if mol_result is not None else mol
        except Exception as e:
            print(f"Error in _swap_atom: {e}")
            return mol
    
    @staticmethod
    def calculate_properties(mol):
        return {
            'Molecular Weight': Descriptors.ExactMolWt(mol),
            'LogP': Descriptors.MolLogP(mol),
            'H-Bond Donors': Descriptors.NumHDonors(mol),
            'H-Bond Acceptors': Descriptors.NumHAcceptors(mol),
            'Topological Polar Surface Area': Descriptors.TPSA(mol)
        }
    
    @staticmethod
    def generate_molecule_image(mol):
        img = Draw.MolToImage(mol, size=(400, 400))
        buffered = io.BytesIO()
        img.save(buffered, format="PNG")
        return base64.b64encode(buffered.getvalue()).decode()

@app.post("/generate_molecule", response_model=MoleculeResponse)
async def generate_molecule(request: MoleculeGenerationRequest):
    try:
        generator = MoleculeGenerator()
        
        # Generate molecule
        generated_mol = generator.generate_random_molecule(
            request.base_smiles, 
            request.num_mutations
        )
        
        # Convert to SMILES
        generated_smiles = Chem.MolToSmiles(generated_mol)
        
        # Calculate properties
        properties = generator.calculate_properties(generated_mol)
        
        # Generate molecule image
        molecule_image = generator.generate_molecule_image(generated_mol)
        
        return MoleculeResponse(
            original_smiles=request.base_smiles,
            generated_smiles=generated_smiles,
            molecular_image=molecule_image,
            properties=properties
        )
    
    except Exception as e:
        logger.error(f"Molecule generation error: {e}")
        raise HTTPException(status_code=400, detail=str(e))

@app.get("/health")
async def health_check():
    return {"status": "healthy"}

if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8000)
