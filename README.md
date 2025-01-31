# Molecular Design Platform 

## Overview
A comprehensive molecular generation and analysis platform that leverages RDKit, FastAPI, and React to enable advanced molecular exploration and design. Integrated from previous standalone molecular generation scripts, this platform provides robust tools for molecular mutation and analysis.

## 🚀 Features
- Interactive molecular generation
- Advanced random molecular mutation strategies:
  - Adding random substituents
  - Modifying bond orders
  - Swapping atom types
- Molecular property calculation
- Base64 encoded molecular image generation
- Visualization of generated molecules

## 🛠 Tech Stack
- **Backend**: 
  - FastAPI
  - RDKit
  - Celery
  - Python
- **Frontend**: 
  - React
  - Axios
  - Tailwind CSS
- **Infrastructure**:
  - Docker
  - PostgreSQL
  - Redis

## 📦 Prerequisites
- Docker
- Docker Compose
- Python 3.9+
- Node.js 18+

## 🔧 Installation & Setup

### Local Development
1. Clone the repository
```bash
git clone https://github.com/bniladridas/molecular_design_platform.git
cd molecular_design_platform
```

2. Start the platform
```bash
docker-compose up --build
```

### Accessing Services
- Frontend: `http://localhost:3000`
- Backend API: `http://localhost:8000/docs`

## 🧪 Usage
1. Enter a base molecule SMILES string
2. Select number of mutations
3. Generate and explore new molecular structures

## 📝 Example Base Molecules
- Aspirin: `CC(=O)OC1=CC=CC=C1C(=O)O`
- Caffeine: `CN1C=NC2=C1C(=O)N(C(=O)N2C)C`
- Salbutamol: `CC(C)(C)NCC(O)C1=CC(=C(C=C1)O)CO`

## 🔬 Mutation Strategies
### Substituent Addition
Randomly add functional groups like:
- Oxygen (`O`)
- Nitrogen (`N`)
- Carbonyl (`C(=O)`)
- Halides (`Cl`, `Br`, `F`)

### Bond Order Modification
Dynamically change bond types between:
- Single bonds
- Aromatic bonds
- Double bonds
- Triple bonds

### Atom Swapping
Replace atoms while maintaining structural integrity:
- Carbon
- Nitrogen
- Oxygen
- Sulfur
- Phosphorus

## Molecular Property Calculation
Automatically calculate key molecular properties:
- Molecular Weight
- LogP (Lipophilicity)
- Number of Atoms
- Number of Bonds
- Ring Count

## 🚧 Limitations
- Random mutations may not always produce chemically valid molecules
- Mutations are probabilistic and may require multiple generations
- Limited to simple structural modifications

## 🤝 Contributing
1. Fork the repository
2. Create your feature branch
3. Commit your changes
4. Push to the branch
5. Create a new Pull Request

## 📄 License
MIT License

## 📞 Contact
Niladri Das - bniladridas@gmail.com

## 🙏 Acknowledgments
- RDKit
- FastAPI
- React Community
- Open-source molecular generation research
