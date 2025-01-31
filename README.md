# Molecular Design Platform 🧬🔬

## Overview
A full-stack molecular generation and analysis platform leveraging RDKit, FastAPI, and React to enable advanced molecular exploration and design.

## 🚀 Features
- Interactive molecular generation
- Random molecular mutation
- Molecular property calculation
- Base molecule customization
- Visualization of generated molecules

## 🛠 Tech Stack
- **Backend**: 
  - FastAPI
  - RDKit
  - SQLAlchemy
  - Celery
- **Frontend**:
  - React
  - Tailwind CSS
- **Infrastructure**:
  - Docker
  - Docker Compose
  - Nginx
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
git clone https://github.com/bniladridas/molecular-design-platform.git
cd molecular-design-platform
```

2. Start the application
```bash
docker-compose up --build
```

3. Access the application
- Frontend: `http://localhost`
- Backend API: `http://localhost/api`

## 🧪 Usage
1. Enter a base molecule SMILES string
2. Select number of mutations
3. Generate and explore new molecular structures

## 📝 Example Base Molecules
- Aspirin: `CC(=O)OC1=CC=CC=C1C(=O)O`
- Caffeine: `CN1C=NC2=C1C(=O)N(C(=O)N2C)C`
- Salbutamol: `CC(C)(C)NCC(O)C1=CC(=C(C=C1)O)CO`

## 🔬 Mutation Strategies
- Adding random substituents
- Modifying bond orders
- Swapping atom types

## 🚧 Limitations
- Random mutations may not always produce chemically valid molecules
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
Your Name - your.email@example.com

## 🙏 Acknowledgments
- RDKit
- FastAPI
- React Community
