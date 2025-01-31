# Molecular Design Platform 

## Overview
A full-stack molecular generation and analysis platform leveraging RDKit, FastAPI, and React to enable advanced molecular exploration and design.

## Features
- Interactive molecular generation
- Random molecular mutation
- Molecular property calculation
- Base molecule customization
- Visualization of generated molecules

## Tech Stack
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

## Prerequisites
- Docker
- Docker Compose
- Python 3.9+
- Node.js 18+

## Installation & Setup

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

## Usage
1. Enter a base molecule SMILES string
2. Select number of mutations
3. Generate and explore new molecular structures

## Example Base Molecules
- Aspirin: `CC(=O)OC1=CC=CC=C1C(=O)O`
- Caffeine: `CN1C=NC2=C1C(=O)N(C(=O)N2C)C`
- Salbutamol: `CC(C)(C)NCC(O)C1=CC(=C(C=C1)O)CO`

## Mutation Strategies
- Adding random substituents
- Modifying bond orders
- Swapping atom types

## Limitations
- Random mutations may not always produce chemically valid molecules
- Limited to simple structural modifications

## User Interface

### Frontend Overview
![Molecular Design Platform Frontend](/frontend_screenshot.png)

The frontend provides an intuitive interface for molecular generation and analysis. Key features include:
- Input field for base molecule SMILES
- Mutation count selector
- Real-time molecule generation
- Visualization of generated molecules
- Property calculation display

### API Documentation
![Swagger API Documentation](/api_docs_screenshot.png)

## API Endpoints

### Molecule Generation
- **Endpoint**: `/generate_molecule`
- **Method**: POST
- **Request Body**:
  ```json
  {
    "base_smiles": "string",  // Base molecule SMILES
    "num_mutations": 3        // Number of mutations to apply
  }
  ```
- **Response**:
  ```json
  {
    "original_smiles": "string",
    "generated_smiles": "string",
    "molecular_image": "base64_encoded_image",
    "properties": {
      "molecular_weight": "float",
      "logp": "float",
      // Other molecular properties
    }
  }
  ```

### Example API Test
You can test the API using `curl`:

```bash
# Generate a molecule
curl -X POST http://localhost:8000/generate_molecule \
     -H "Content-Type: application/json" \
     -d '{
         "base_smiles": "CC(=O)OC1=CC=CC=C1C(=O)O",  # Aspirin
         "num_mutations": 3
     }'
```

## Performance and Scalability
- Built with FastAPI for high-performance async handling
- Celery worker for background task processing
- Redis for task queuing and result storage
- Dockerized microservices architecture

## Contributing
1. Fork the repository
2. Create your feature branch
3. Commit your changes
4. Push to the branch
5. Create a new Pull Request

## License
MIT License

## 📞 Contact

Contact me on [LinkedIn](https://www.linkedin.com/in/bniladridas/)

## Acknowledgments
- RDKit
- FastAPI
- React Community
