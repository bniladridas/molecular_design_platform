import pytest
from rdkit import Chem
from app.main import MoleculeGenerator

def test_generate_random_molecule():
    base_smiles = 'CC(=O)OC1=CC=CC=C1C(=O)O'  # Aspirin
    generator = MoleculeGenerator()
    
    # Generate a molecule
    generated_mol = generator.generate_random_molecule(base_smiles)
    
    # Check that a valid molecule was generated
    assert generated_mol is not None
    assert isinstance(generated_mol, Chem.Mol)
    
    # Check that the generated molecule is different from the base molecule
    generated_smiles = Chem.MolToSmiles(generated_mol)
    assert generated_smiles != base_smiles

def test_molecule_properties():
    base_smiles = 'CC(=O)OC1=CC=CC=C1C(=O)O'  # Aspirin
    generator = MoleculeGenerator()
    
    # Generate a molecule
    generated_mol = generator.generate_random_molecule(base_smiles)
    
    # Calculate properties
    properties = generator.calculate_properties(generated_mol)
    
    # Check that properties are calculated
    assert 'Molecular Weight' in properties
    assert 'LogP' in properties
    assert 'H-Bond Donors' in properties
    assert 'H-Bond Acceptors' in properties
    assert 'Topological Polar Surface Area' in properties
    
    # Check that properties are numeric
    for prop_value in properties.values():
        assert isinstance(prop_value, (int, float))

def test_mutation_strategies():
    base_smiles = 'CC(=O)OC1=CC=CC=C1C(=O)O'  # Aspirin
    
    # Test multiple mutation strategies
    mutation_methods = [
        MoleculeGenerator._add_random_substituent,
        MoleculeGenerator._modify_bond_order,
        MoleculeGenerator._swap_atom
    ]
    
    for method in mutation_methods:
        mol = Chem.MolFromSmiles(base_smiles)
        mutated_mol = method(mol)
        
        assert mutated_mol is not None
        assert isinstance(mutated_mol, Chem.Mol)
        assert Chem.MolToSmiles(mutated_mol) != base_smiles
