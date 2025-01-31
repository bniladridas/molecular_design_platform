from app.celery_worker import app
import logging
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, Crippen, Draw
import random
import base64
import io

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

@app.task
def add(x, y):
    return x + y

@app.task
def multiply(x, y):
    return x * y

@app.task
def generate_random_molecule(base_smiles, num_mutations=3):
    """
    Generate a new molecule by mutating a base molecule
    
    :param base_smiles: Base molecule SMILES string
    :param num_mutations: Number of random mutations to apply
    :return: Dictionary with molecule details
    """
    try:
        # Parse base molecule
        mol = Chem.MolFromSmiles(base_smiles)
        if not mol:
            raise ValueError(f"Invalid SMILES: {base_smiles}")

        # Apply mutations
        for _ in range(num_mutations):
            mol = _random_mutation(mol)

        # Generate SMILES and properties
        generated_smiles = Chem.MolToSmiles(mol)
        properties = calculate_properties(mol)
        molecular_image = generate_molecule_image(mol)

        return {
            "original_smiles": base_smiles,
            "generated_smiles": generated_smiles,
            "molecular_image": molecular_image,
            "properties": properties
        }
    except Exception as e:
        logger.error(f"Error generating molecule: {e}")
        raise

def _random_mutation(mol):
    """Apply a random mutation to the molecule"""
    mutation_strategies = [
        _add_random_substituent,
        _modify_bond_order,
        _swap_atom
    ]
    strategy = random.choice(mutation_strategies)
    return strategy(mol)

def _add_random_substituent(mol):
    """Add a random substituent to the molecule"""
    substituents = [
        'O', 'N', 'C(=O)', 'C(=N)', 'Cl', 'Br', 'F', 'S'
    ]
    substituent = random.choice(substituents)
    
    # Find a random atom to attach the substituent
    atom_indices = list(range(mol.GetNumAtoms()))
    random.shuffle(atom_indices)
    
    for idx in atom_indices:
        atom = mol.GetAtomWithIdx(idx)
        if atom.GetTotalNumHs() > 0:  # Ensure the atom can accept a substituent
            try:
                new_mol = Chem.RWMol(mol)
                new_substituent = Chem.MolFromSmarts(substituent)
                new_mol.ReplaceSubstructs(
                    Chem.MolFromSmarts(f'[#1:1]'),
                    new_substituent,
                    idx
                )
                return new_mol
            except Exception:
                continue
    
    return mol

def _modify_bond_order(mol):
    """Modify bond order of a random bond"""
    bond_indices = list(range(mol.GetNumBonds()))
    random.shuffle(bond_indices)
    
    for idx in bond_indices:
        bond = mol.GetBondWithIdx(idx)
        current_order = bond.GetBondTypeAsDouble()
        
        # Cycle through bond orders
        new_orders = [1.0, 1.5, 2.0, 3.0]
        new_orders.remove(current_order)
        new_order = random.choice(new_orders)
        
        try:
            new_mol = Chem.RWMol(mol)
            new_mol.GetBondWithIdx(idx).SetBondType(
                Chem.BondType(new_order)
            )
            return new_mol
        except Exception:
            continue
    
    return mol

def _swap_atom(mol):
    """Swap an atom with another atom type"""
    atom_types = ['C', 'N', 'O', 'S', 'P']
    atom_indices = list(range(mol.GetNumAtoms()))
    random.shuffle(atom_indices)
    
    for idx in atom_indices:
        atom = mol.GetAtomWithIdx(idx)
        if atom.GetSymbol() not in ['H']:  # Avoid swapping hydrogen
            current_type = atom.GetSymbol()
            possible_types = [t for t in atom_types if t != current_type]
            new_type = random.choice(possible_types)
            
            try:
                new_mol = Chem.RWMol(mol)
                new_mol.GetAtomWithIdx(idx).SetAtomicNum(
                    Chem.GetPeriodicTable().GetAtomicNumber(new_type)
                )
                return new_mol
            except Exception:
                continue
    
    return mol

def calculate_properties(mol):
    """Calculate molecular properties"""
    return {
        "molecular_weight": Descriptors.ExactMolWt(mol),
        "logp": Crippen.MolLogP(mol),
        "num_atoms": mol.GetNumAtoms(),
        "num_bonds": mol.GetNumBonds(),
        "num_rings": Descriptors.RingCount(mol)
    }

def generate_molecule_image(mol):
    """Generate a base64 encoded image of the molecule"""
    try:
        img = Draw.MolToImage(mol, size=(300, 300))
        buffered = io.BytesIO()
        img.save(buffered, format="PNG")
        return base64.b64encode(buffered.getvalue()).decode('utf-8')
    except Exception as e:
        logger.error(f"Error generating molecule image: {e}")
        return None

# Example usage for testing
@app.task
def test_molecule_generation():
    base_molecules = [
        "CC(=O)OC1=CC=CC=C1C(=O)O",  # Aspirin
        "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",  # Caffeine
        "CC(C)(C)NCC(O)C1=CC(=C(C=C1)O)CO"  # Salbutamol
    ]
    
    for base_smiles in base_molecules:
        result = generate_random_molecule(base_smiles)
        logger.info(f"Generated molecule from {base_smiles}: {result['generated_smiles']}")
        logger.info(f"Properties: {result['properties']}")
