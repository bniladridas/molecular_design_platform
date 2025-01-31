# Molecular Design Platform

## Overview
The **Molecular Design Platform** is a containerized application designed for molecular simulation and mutation analysis. It utilizes a microservices architecture, including a backend API, frontend interface, database, Redis, Celery worker, and Nginx for proxying requests. The platform provides a seamless environment for molecular computations, leveraging Docker and Docker Compose for easy deployment.

## Prerequisites
Ensure you have the following installed before running the platform:
- [Docker](https://docs.docker.com/get-docker/)
- [Docker Compose](https://docs.docker.com/compose/install/)

## Running the Platform
To start the Molecular Design Platform, run the following command:

```sh
$ docker-compose up --build
```

### Expected Output
During startup, you might see a warning regarding the `version` attribute in `docker-compose.yml`. You can ignore it or remove the version attribute to avoid confusion.

Once built successfully, the following services will be up and running:

| Service         | Container Name                              | Port Mapping        |
|----------------|------------------------------------------|---------------------|
| Backend API    | `molecular_design_platform-backend-1`   | `0.0.0.0:8000->8000/tcp`  |
| Frontend       | `molecular_design_platform-frontend-1`  | `0.0.0.0:3000->3000/tcp`  |
| Database       | `molecular_design_platform-database-1`  | `5432/tcp`         |
| Redis         | `molecular_design_platform-redis-1`     | `0.0.0.0:6379->6379/tcp`  |
| Celery Worker  | `molecular_design_platform-celery_worker-1` | No exposed ports |
| Nginx Proxy    | `molecular_design_platform-nginx-1`    | `0.0.0.0:80->80/tcp`  |

## Checking Running Services
You can verify the running services with:

```sh
$ docker-compose ps
```

## Accessing the Platform
Once all services are up, access the following:

- **Backend API Documentation**: [http://localhost:8000/docs](http://localhost:8000/docs)
- **Frontend Web Interface**: [http://localhost:3000](http://localhost:3000)

## API Testing
To test the backend API, you can run:

```sh
$ curl http://localhost:8000/docs
```

## Frontend Check
Ensure the frontend is running properly by executing:

```sh
$ curl http://localhost:3000
```

If the output contains an HTML document with `<title>Molecular Design Platform</title>`, the frontend is up.

## Troubleshooting
- If any service fails to start, check logs using:
  ```sh
  $ docker-compose logs -f <service-name>
  ```
- Ensure Docker Desktop is running (if using macOS/Windows).
- Restart containers if needed:
  ```sh
  $ docker-compose down && docker-compose up --build
  ```

## Contributing
Contributions are welcome! Feel free to submit issues or pull requests to improve the Molecular Design Platform.

## License
This project is licensed under the MIT License.
